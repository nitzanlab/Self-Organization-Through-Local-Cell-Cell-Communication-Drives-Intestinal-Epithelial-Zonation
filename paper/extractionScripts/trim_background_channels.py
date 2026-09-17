"""Drop unread channels from the bundle's background TIFFs.

crop_monolayer.py reads channel 3 of pasadena_roi1; transplant_neighborhoods.py
reads channels 2 (GFP) and 3 (DAPI) of the two nov23 files. Everything else in
these 4-channel stacks is dead weight in the archive. The full acquisitions
remain on Dropbox under shared_yael/sg/.

Reads and writes one plane at a time (~3.9 GB peak for the largest file), then
verifies by comparing a SHA-256 of every kept plane against the original.
"""
import sys, os, time, gc, hashlib
import numpy as np
import tifffile

BUNDLE = "/Users/yaelheyman/Documents/zonation_data_bundle/backgrounds"
KEEP = {                      # folder -> source channels to keep, in output order
    "pasadena_roi1":   [3],
    "nov23_72hr_roi1": [2, 3],
    "nov23_12hr_roi2": [2, 3],
}

def plane_hash(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16]

def trim(folder, keep):
    src = os.path.join(BUNDLE, folder, "hyb_background_aligned.tiff")
    dst = src + ".trimmed"
    src_hashes = []
    with tifffile.TiffFile(src) as tf:
        C, H, W = tf.series[0].shape
        assert tf.series[0].dtype == np.uint16 and len(tf.pages) == C
        print(f"[{folder}] {(C,H,W)} uint16 -> keep {keep}", flush=True)
        t0 = time.time()
        with tifffile.TiffWriter(dst, bigtiff=True) as tw:
            for oi, ch in enumerate(keep):
                plane = tf.pages[ch].asarray()
                src_hashes.append(plane_hash(plane))
                # contiguous=True appends into ONE series, so tifffile.imread()
                # returns (n_kept, H, W) -- which is what the plotting code indexes.
                tw.write(plane, photometric="minisblack",
                         rowsperstrip=1024, contiguous=True)
                print(f"  ch{ch} -> out[{oi}]  sha={src_hashes[-1]}  "
                      f"({time.time()-t0:.0f}s)", flush=True)
                del plane; gc.collect()
    print(f"[{folder}] {os.path.getsize(dst)/2**30:.1f} GiB "
          f"(was {os.path.getsize(src)/2**30:.1f} GiB)", flush=True)
    return src, dst, src_hashes

def verify(src, dst, keep, src_hashes):
    with tifffile.TiffFile(src) as a, tifffile.TiffFile(dst) as b:
        C, H, W = a.series[0].shape
        exp = (len(keep), H, W) if len(keep) > 1 else (H, W)
        # what the plotting code actually sees
        assert len(b.series) == 1, f"{len(b.series)} series, expected 1"
        assert b.series[0].shape == exp, f"shape {b.series[0].shape} != {exp}"
        assert b.series[0].dtype == np.uint16, b.series[0].dtype
        assert len(b.pages) == len(keep), f"{len(b.pages)} pages != {len(keep)}"
        for oi, ch in enumerate(keep):
            plane = b.pages[oi].asarray()
            h = plane_hash(plane)
            del plane; gc.collect()
            assert h == src_hashes[oi], f"ch{ch}: {h} != {src_hashes[oi]}"
            print(f"  verified out[{oi}] == source ch{ch}  sha={h}", flush=True)
    return True

if __name__ == "__main__":
    for folder in sys.argv[1:] or list(KEEP):
        src, dst, hs = trim(folder, KEEP[folder])
        verify(src, dst, KEEP[folder], hs)
        print(f"[{folder}] VERIFIED\n", flush=True)
