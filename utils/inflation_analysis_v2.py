# inflation_analysis_v2.py
#
# Adds basal-border support to v1:
#   - set_basal_mask()          : clips dilation to the tissue ROI (between luminal and basal borders)
#   - compute_average_thickness(): ROI area / luminal perimeter * pixel_size
#   - compute_local_widths()    : per-point width profile along the basal contour (KDTree)
#   - save_border_overlay()     : tissue image with cyan luminal + red basal contours
#   - save_width_profile_overlay(): two-panel tissue+lines / histogram figure

import os
import cv2
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from skimage.measure import regionprops, label, find_contours
from scipy.spatial import KDTree


class InflationAnalysis:
    """
    Dilation-based ring analysis with optional two-border (luminal + basal) support.

    Typical usage (two-border):
        ia = InflationAnalysis(output_folder="out", kernel_size=20, inflation_iterations=10, pixel2nm=1)
        ia.set_mask(luminal_mask_u8)
        ia.set_basal_mask(basal_mask_u8)
        ia.save_border_overlay(tissue_img, "out/border_overlay.png")
        ia.inflate_and_save()
        width_df, stats = ia.compute_local_widths(n_samples=500, min_width_um=50,
                                                   save_csv=True, iter_dir="out/inflation_num_iterations_10")
        ia.save_width_profile_overlay(tissue_img, width_df, "out/.../width_profile_overlay.png")
    """

    def __init__(
        self,
        output_folder: str,
        kernel_size: int = 40,
        inflation_iterations: int = 10,
        pixel2nm: float = 1,
        kernel_shape: str = "rect",
    ):
        self.output_folder = output_folder
        os.makedirs(self.output_folder, exist_ok=True)

        self.kernel_size       = int(kernel_size)
        self.inflation_iterations = int(inflation_iterations)
        self.pixel2nm          = float(pixel2nm)
        self.kernel_shape      = kernel_shape.lower()

        self.mask        = None   # luminal (inner) mask, uint8 0/255
        self.basal_mask  = None   # basal (outer) mask,   uint8 0/255

        # optional cell / spot assignment (kept for compatibility with v1)
        self.cell_by_gene    = None
        self.scaled_centers  = None

        # saved after inflate_and_save for use by other methods
        self._luminal_mask_orig = None
        self._iter_dir          = None

    # ------------------------------------------------------------------ masks

    @staticmethod
    def poly2mask_from_image_shape(vertices_xy, H, W, value=255):
        """
        vertices_xy : (N,2) in (x,y) = (col, row) pixel coords
        Returns uint8 mask (H,W) filled with `value` inside the polygon.
        """
        v = np.asarray(vertices_xy, dtype=float)
        if v.ndim != 2 or v.shape[1] != 2:
            raise ValueError("vertices_xy must be (N,2) array-like of (x,y)")
        x = np.clip(np.round(v[:, 0]), 0, W - 1)
        y = np.clip(np.round(v[:, 1]), 0, H - 1)
        pts = np.stack([x, y], axis=1).astype(np.int32).reshape((-1, 1, 2))
        mask = np.zeros((H, W), dtype=np.uint8)
        cv2.fillPoly(mask, [pts], color=int(value))
        return mask

    @staticmethod
    def load_polygon_csv_napari(polygon_path):
        """
        Reads Napari-exported CSV with columns 'axis-0' (row=y) and 'axis-1' (col=x).
        Returns vertices_xy as (x, y) = (axis-1, axis-0).
        """
        df = pd.read_csv(polygon_path)
        if not {"axis-0", "axis-1"}.issubset(df.columns):
            raise KeyError("Expected columns 'axis-0' and 'axis-1' in polygon CSV")
        return df[["axis-1", "axis-0"]].to_numpy()

    def set_mask(self, mask_u8):
        """Set the luminal (inner) polygon mask."""
        m = np.asarray(mask_u8)
        if m.ndim != 2:
            raise ValueError("mask must be 2D")
        self.mask = (m > 0).astype(np.uint8) * 255

    def set_basal_mask(self, mask_u8):
        """
        Set the basal (outer) polygon mask.
        Dilation at each step is clipped to this region, so rings never grow
        beyond the basal border.  Rings stop automatically when no new ROI
        pixels remain; a final catch-all ring covers any leftover pixels.
        """
        m = np.asarray(mask_u8)
        if m.ndim != 2:
            raise ValueError("basal_mask must be 2D")
        self.basal_mask = (m > 0).astype(np.uint8) * 255

    def set_cells(self, cell_by_gene_df: pd.DataFrame, centers_xy: np.ndarray, id_col="Id"):
        """Optional spot / cell assignment (same API as v1)."""
        if id_col not in cell_by_gene_df.columns:
            raise KeyError(f"Expected column '{id_col}' in cell_by_gene_df")
        self.cell_by_gene   = cell_by_gene_df.reset_index(drop=True).copy()
        self.scaled_centers = np.asarray(centers_xy, dtype=float)
        if self.scaled_centers.shape[0] != len(self.cell_by_gene):
            raise ValueError("centers_xy length must match number of rows in cell_by_gene_df")

    # ------------------------------------------------------------------ metrics

    def calculate_inflation_ring_width(self, prev_u8, nxt_u8) -> float:
        """Ring width (µm) = ring_area / prev_perimeter * pixel_size."""
        prevb = (prev_u8 > 0)
        nxtb  = (nxt_u8  > 0)
        ring  = nxtb & (~prevb)
        if not ring.any():
            return 0.0

        lbl = label(prevb, connectivity=2)
        if lbl.max() == 0:
            return 0.0

        total_ring_area = 0.0
        total_perimeter  = 0.0
        for pr in regionprops(lbl):
            comp = (lbl == pr.label)
            ra   = np.count_nonzero(ring & comp)
            if ra <= 0 or pr.perimeter <= 0:
                continue
            total_ring_area += ra
            total_perimeter  += pr.perimeter

        if total_ring_area <= 0 or total_perimeter <= 0:
            return 0.0
        return float(total_ring_area / total_perimeter) * (self.pixel2nm / 1000.0)

    def compute_average_thickness(self) -> float:
        """
        Average tissue thickness (µm) = ROI area / luminal perimeter * pixel_size.
        ROI = basal_mask & ~luminal_mask_orig.
        Requires both masks to be set and inflate_and_save() to have been called first.
        """
        if self._luminal_mask_orig is None:
            raise RuntimeError("Call inflate_and_save() before compute_average_thickness().")
        if self.basal_mask is None:
            raise RuntimeError("No basal mask set.")

        lum  = (self._luminal_mask_orig > 0)
        bas  = (self.basal_mask > 0)
        roi  = bas & (~lum)

        roi_area_px = float(np.count_nonzero(roi))

        lbl  = label(lum, connectivity=2)
        perim_px = sum(pr.perimeter for pr in regionprops(lbl))
        if perim_px <= 0:
            return 0.0

        thickness_um = (roi_area_px / perim_px) * (self.pixel2nm / 1000.0)
        return thickness_um

    # ------------------------------------------------------------------ main loop

    def inflate_and_save(
        self,
        ring_col="InflationRing",
        iter_folder_prefix="inflation",
        mask_prefix="inflated",
        save_diff_masks=True,
    ):
        """
        Runs iterative dilation.

        When basal_mask is set:
          - each dilated ring is clipped to basal_mask
          - loop stops when no new ROI pixels remain (or max iterations reached)
          - any unassigned ROI pixels become a final catch-all ring

        Saves per-iteration masks, diff masks, ring_widths.csv, and (when basal_mask
        is set) section_thickness.csv.
        """
        if self.mask is None:
            raise ValueError("No mask set. Call set_mask() first.")

        iter_dir = os.path.join(
            self.output_folder,
            f"{iter_folder_prefix}_num_iterations_{self.inflation_iterations}"
        )
        os.makedirs(iter_dir, exist_ok=True)
        self._iter_dir = iter_dir

        # save original luminal mask for thickness / width calculations
        self._luminal_mask_orig = self.mask.copy()

        if self.cell_by_gene is not None:
            self.cell_by_gene[ring_col] = 0

        prev = (self.mask > 0).astype(np.uint8) * 255

        if self.kernel_shape == "ellipse":
            kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (self.kernel_size, self.kernel_size))
        else:
            kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (self.kernel_size, self.kernel_size))

        # ROI pixels that still need to be assigned (only used in two-border mode)
        if self.basal_mask is not None:
            bas_bool = (self.basal_mask > 0)
            lum_bool = (self._luminal_mask_orig > 0)
            roi_bool = bas_bool & (~lum_bool)
            unassigned = roi_bool.copy()
        else:
            unassigned = None

        ring_widths = []
        all_coords  = []

        for i in range(self.inflation_iterations):
            print(f"Inflation iteration {i+1}")
            nxt = cv2.dilate(prev, kernel, iterations=1)
            nxt = (nxt > 0).astype(np.uint8) * 255

            # clip to basal mask if set
            if self.basal_mask is not None:
                nxt = np.where(self.basal_mask > 0, nxt, 0).astype(np.uint8)

            diff_bool = (nxt > 0) & ~(prev > 0)

            if not diff_bool.any():
                print("  → no further growth, stopping.")
                break

            # remove assigned pixels from unassigned pool
            if unassigned is not None:
                unassigned &= ~diff_bool

            w = self.calculate_inflation_ring_width(prev, nxt)
            ring_widths.append({"Ring Width (um)": w})

            ring_number = i + 1
            coords = []
            if self.scaled_centers is not None and self.cell_by_gene is not None:
                for idx, (x, y) in enumerate(self.scaled_centers):
                    xi, yi = int(round(x)), int(round(y))
                    if 0 <= yi < diff_bool.shape[0] and 0 <= xi < diff_bool.shape[1]:
                        if diff_bool[yi, xi]:
                            cid = self.cell_by_gene.iloc[idx]["Id"]
                            coords.append((cid, x, y))
                ring_ids = [cid for cid, _, _ in coords]
                self.cell_by_gene.loc[self.cell_by_gene["Id"].isin(ring_ids), ring_col] = ring_number
            all_coords.append(coords)

            Image.fromarray(nxt).save(os.path.join(iter_dir, f"{mask_prefix}_mask_{ring_number}.png"))
            if save_diff_masks:
                Image.fromarray(diff_bool.astype(np.uint8) * 255).save(
                    os.path.join(iter_dir, f"{mask_prefix}_diff_{ring_number}.png"))
            pd.DataFrame(coords, columns=["cell_id", "X", "Y"]).to_csv(
                os.path.join(iter_dir, f"{mask_prefix}_coords_{ring_number}.csv"), index=False)

            prev = nxt
            self.mask = nxt

            # in two-border mode: stop if all ROI pixels are covered
            if unassigned is not None and not unassigned.any():
                print("  → all ROI pixels assigned, stopping.")
                break

        # catch-all final ring for any leftover ROI pixels (two-border mode only)
        if unassigned is not None and unassigned.any():
            ring_number = len(ring_widths) + 1
            print(f"  Final catch-all ring {ring_number}: {np.count_nonzero(unassigned)} leftover pixels.")
            ring_widths.append({"Ring Width (um)": 0.0})   # width undefined for catch-all

            if save_diff_masks:
                Image.fromarray(unassigned.astype(np.uint8) * 255).save(
                    os.path.join(iter_dir, f"{mask_prefix}_diff_{ring_number}.png"))
            Image.fromarray(((prev.astype(bool) | unassigned).astype(np.uint8) * 255)).save(
                os.path.join(iter_dir, f"{mask_prefix}_mask_{ring_number}.png"))
            pd.DataFrame([], columns=["cell_id", "X", "Y"]).to_csv(
                os.path.join(iter_dir, f"{mask_prefix}_coords_{ring_number}.csv"), index=False)
            all_coords.append([])

        # ---- write ring_widths.csv ----
        pd.DataFrame(ring_widths).to_csv(os.path.join(iter_dir, "ring_widths.csv"), index=False)

        if self.cell_by_gene is not None:
            self.cell_by_gene.to_csv(os.path.join(iter_dir, "cell_by_gene_with_ring.csv"), index=False)

        # ---- section_thickness.csv (two-border mode) ----
        if self.basal_mask is not None:
            avg_thickness = self.compute_average_thickness()
            n_rings = len(ring_widths)
            pd.DataFrame([{
                "n_rings":          n_rings,
                "avg_thickness_um": avg_thickness,
                "kernel_size_px":   self.kernel_size,
                "pixel2nm":         self.pixel2nm,
            }]).to_csv(os.path.join(iter_dir, "section_thickness.csv"), index=False)
            print(f"\nAverage section thickness: {avg_thickness:.2f} µm  ({n_rings} rings)")

        # ---- final overlay ----
        from matplotlib.colors import hsv_to_rgb
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(self.mask, cmap="gray", origin="lower")

        if all_coords:
            hsv  = [(k / max(1, len(all_coords)), 1, 1) for k in range(len(all_coords))]
            cols = hsv_to_rgb(hsv)
            for k, coords in enumerate(all_coords):
                if not coords:
                    continue
                _, xs, ys = zip(*coords)
                ax.scatter(xs, ys, c=[cols[k]], s=10, label=f"Ring {k+1}")
            ax.legend()
        else:
            ax.set_title("No rings inflated; showing final mask")

        ax.axis("off")
        plt.tight_layout()
        plt.savefig(os.path.join(iter_dir, "final_rings.png"))
        plt.close()

    # ------------------------------------------------------------------ width profile

    def compute_local_widths(
        self,
        n_samples: int = 500,
        min_width_um: float = None,
        save_csv: bool = True,
        iter_dir: str = None,
    ):
        """
        Per-point local width along the tissue section.

        Samples n_samples points along the basal border contour and finds the
        nearest point on the luminal border contour using a KDTree.

        Parameters
        ----------
        n_samples    : number of points to sample along the basal contour
        min_width_um : exclude points with width < this from stats (cut-tissue filter).
                       Set to None to include all points.
        save_csv     : if True, saves width_profile.csv and width_stats.csv to iter_dir
        iter_dir     : folder for CSV output (defaults to self._iter_dir)

        Returns
        -------
        df    : DataFrame with columns basal_x, basal_y, lum_x, lum_y,
                width_px, width_um, included
        stats : dict with mean_um, std_um, median_um, n_total, n_included
        """
        if self._luminal_mask_orig is None or self.basal_mask is None:
            raise RuntimeError("Both luminal and basal masks must be set and inflate_and_save() called first.")

        iter_dir = iter_dir or self._iter_dir
        if iter_dir:
            os.makedirs(iter_dir, exist_ok=True)

        # extract contour points
        lum_contours  = find_contours((self._luminal_mask_orig > 0).astype(float), 0.5)
        bas_contours  = find_contours((self.basal_mask > 0).astype(float), 0.5)

        if not lum_contours or not bas_contours:
            raise ValueError("Could not find contours in one or both masks.")

        # use the longest contour for each border
        lum_pts = max(lum_contours, key=len)   # (N, 2) in (row, col) = (y, x)
        bas_pts = max(bas_contours, key=len)

        # sample n_samples evenly spaced points along basal contour
        idx = np.round(np.linspace(0, len(bas_pts) - 1, n_samples)).astype(int)
        bas_sampled = bas_pts[idx]   # (n_samples, 2) as (row, col)

        # build KDTree on luminal contour (row, col)
        tree = KDTree(lum_pts)
        dists, near_idx = tree.query(bas_sampled)

        lum_nearest = lum_pts[near_idx]   # (n_samples, 2)

        width_px = dists
        width_um = dists * (self.pixel2nm / 1000.0)

        included = np.ones(n_samples, dtype=bool)
        if min_width_um is not None:
            # when pixel2nm≤1 no real calibration is set — treat min_width_um as pixels
            if self.pixel2nm <= 1.0:
                included = width_px >= min_width_um
            else:
                included = width_um >= min_width_um

        df = pd.DataFrame({
            "basal_y":  bas_sampled[:, 0],
            "basal_x":  bas_sampled[:, 1],
            "lum_y":    lum_nearest[:, 0],
            "lum_x":    lum_nearest[:, 1],
            "width_px": width_px,
            "width_um": width_um,
            "included": included,
        })

        w_inc = width_um[included]
        stats = {
            "mean_um":   float(w_inc.mean())   if w_inc.size else 0.0,
            "std_um":    float(w_inc.std())    if w_inc.size else 0.0,
            "median_um": float(np.median(w_inc)) if w_inc.size else 0.0,
            "n_total":   int(n_samples),
            "n_included": int(included.sum()),
        }

        if save_csv and iter_dir:
            df.to_csv(os.path.join(iter_dir, "width_profile.csv"), index=False)
            pd.DataFrame([stats]).to_csv(os.path.join(iter_dir, "width_stats.csv"), index=False)

        return df, stats

    # ------------------------------------------------------------------ overlays

    def save_border_overlay(self, tissue_img, output_path: str):
        """
        Save the tissue image with clearly visible borders:
          - bright yellow : luminal (inner) border
          - bright cyan   : basal   (outer) border
        Line thickness scales with image size.
        Works whether called before or after inflate_and_save().
        """
        img = np.asarray(tissue_img).copy()
        if img.ndim == 2:
            img = np.stack([img] * 3, axis=-1)
        elif img.shape[2] == 4:
            img = img[:, :, :3]

        H, W = img.shape[:2]
        thickness = max(4, H // 200)   # scales: ~5 px on 1000px image, ~10 on 2000px

        # use self.mask as luminal source (valid before inflate_and_save too)
        lum_mask = self._luminal_mask_orig if self._luminal_mask_orig is not None else self.mask

        fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
        ax.imshow(img)

        def plot_contour(mask, color, label):
            cnts, _ = cv2.findContours((mask > 0).astype(np.uint8),
                                       cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            first = True
            for c in cnts:
                c = c.squeeze()
                if c.ndim < 2:
                    continue
                xs = np.append(c[:, 0], c[0, 0])
                ys = np.append(c[:, 1], c[0, 1])
                ax.plot(xs, ys, color=color, linewidth=thickness * 0.4,
                        label=label if first else "_nolegend_", solid_capstyle="round")
                first = False

        if lum_mask is not None:
            plot_contour(lum_mask, "#ffe600", "Luminal border")   # bright yellow
        if self.basal_mask is not None:
            plot_contour(self.basal_mask, "#00e5ff", "Basal border")  # bright cyan

        ax.legend(loc="upper right", fontsize=14,
                  facecolor="black", labelcolor="white", framealpha=0.7)
        ax.axis("off")
        plt.tight_layout(pad=0)
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Border overlay saved: {output_path}")

    def save_width_profile_overlay(
        self,
        tissue_img,
        width_df: pd.DataFrame,
        output_path: str,
        n_lines: int = 300,
        cmap: str = "plasma",
    ):
        """
        Two-panel figure:
          Left  : tissue image with luminal/basal borders and width lines
                  (lines colored by width_um; excluded points drawn in grey)
          Right : histogram of included widths

        Parameters
        ----------
        tissue_img  : ndarray (H,W) or (H,W,3)
        width_df    : DataFrame returned by compute_local_widths()
        output_path : save path for the PNG
        n_lines     : max number of lines to draw (subsample if more rows)
        cmap        : colormap name for width coloring
        """
        img = np.asarray(tissue_img)
        if img.ndim == 2:
            img = np.stack([img] * 3, axis=-1)
        elif img.shape[2] == 4:
            img = img[:, :, :3]
        if img.max() > 1.5:
            img = img.astype(np.float32) / 255.0

        included = width_df["included"].to_numpy()
        width_px_vals = width_df["width_px"].to_numpy()
        width_um_vals = width_df["width_um"].to_numpy()

        # use pixel widths for display when no real calibration set
        use_px = self.pixel2nm <= 1.0
        display_widths = width_px_vals if use_px else width_um_vals
        unit_label = "px" if use_px else "µm"

        w_inc = display_widths[included]
        vmin = w_inc.min() if included.any() else display_widths.min()
        vmax = w_inc.max() if included.any() else display_widths.max()
        if vmin == vmax:
            vmin, vmax = vmin * 0.9, vmax * 1.1 + 1
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap_obj = cm.get_cmap(cmap)

        # subsample lines
        idx = np.round(np.linspace(0, len(width_df) - 1, min(n_lines, len(width_df)))).astype(int)
        df_sub = width_df.iloc[idx]
        disp_sub = display_widths[idx]
        inc_sub  = included[idx]

        # line thickness scales with image size
        H_img = img.shape[0]
        lw = max(1.0, H_img / 800)

        fig, axes = plt.subplots(1, 2, figsize=(16, 7))

        # --- Left panel ---
        axes[0].imshow(img, origin="upper")

        for i, (_, row) in enumerate(df_sub.iterrows()):
            xs = [row["basal_x"], row["lum_x"]]
            ys = [row["basal_y"], row["lum_y"]]
            if inc_sub[i]:
                color = cmap_obj(norm(disp_sub[i]))
                axes[0].plot(xs, ys, color=color, linewidth=lw, alpha=0.85, solid_capstyle="round")
            else:
                axes[0].plot(xs, ys, color=(0.3, 0.3, 0.3), linewidth=lw * 0.6,
                             alpha=0.4, solid_capstyle="round")

        # draw borders
        if self._luminal_mask_orig is not None:
            cnts, _ = cv2.findContours((self._luminal_mask_orig > 0).astype(np.uint8),
                                       cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            for c in cnts:
                c = c.squeeze()
                if c.ndim == 2:
                    axes[0].plot(c[:, 0], c[:, 1], color="cyan", linewidth=1.5, label="luminal")
        if self.basal_mask is not None:
            cnts, _ = cv2.findContours((self.basal_mask > 0).astype(np.uint8),
                                       cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            for c in cnts:
                c = c.squeeze()
                if c.ndim == 2:
                    axes[0].plot(c[:, 0], c[:, 1], color="red", linewidth=1.5, label="basal")

        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm)
        sm.set_array([])
        plt.colorbar(sm, ax=axes[0], label=f"Width ({unit_label})", fraction=0.03, pad=0.02)
        axes[0].axis("off")
        axes[0].set_title("Local tissue width", fontsize=12)

        # --- Right panel: histogram ---
        # fall back to all points if none included
        plot_widths = w_inc if included.any() else display_widths
        plot_label  = "included" if included.any() else "all points (none passed filter)"

        axes[1].hist(plot_widths, bins=30, color="#2196a8", edgecolor="white", linewidth=0.5,
                     label=plot_label)
        if plot_widths.size > 0:
            axes[1].axvline(plot_widths.mean(), color="black", lw=1.4, ls="--",
                            label=f"Mean={plot_widths.mean():.1f} {unit_label}")
            axes[1].axvline(np.median(plot_widths), color="#e63946", lw=1.4, ls="-.",
                            label=f"Median={np.median(plot_widths):.1f} {unit_label}")
            cv_str = f"{plot_widths.std()/plot_widths.mean()*100:.1f}%" if plot_widths.mean() > 0 else "n/a"
        else:
            cv_str = "n/a"
        axes[1].set_xlabel(f"Width ({unit_label})", fontsize=11)
        axes[1].set_ylabel("Count", fontsize=11)
        axes[1].set_title(
            f"Width distribution\n"
            f"n={int(included.sum())} included / {len(included)} total  |  CV={cv_str}",
            fontsize=11,
        )
        axes[1].legend(fontsize=9)
        axes[1].spines["top"].set_visible(False)
        axes[1].spines["right"].set_visible(False)

        plt.tight_layout()
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Width profile overlay saved: {output_path}")
