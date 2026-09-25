# Self-Organization Through Local Cell-Cell Communication Drives Intestinal Epithelial Zonation

Code to reproduce the analyses and figure panels in the paper.

## 1. Clone the repository

```bash
git clone https://github.com/nitzanlab/Self-Organization-Through-Local-Cell-Cell-Communication-Drives-Intestinal-Epithelial-Zonation.git
cd Self-Organization-Through-Local-Cell-Cell-Communication-Drives-Intestinal-Epithelial-Zonation
```

## 2. Install dependencies

```bash
pip install -r utils/requirements.txt
```

## 3. Download the data

Create one directory to hold everything — its path is what you set in step 4.

### (a) This paper's data — Zenodo

**DOI: [10.5281/zenodo.17956269](https://doi.org/10.5281/zenodo.17956269)**

Download every file from that record into one directory, then follow the
`RECONSTRUCT.md` included there. It explains how to concatenate the files distributed
as `.part` pieces (with an MD5 for each reassembled file), how to unzip the archives in
place, and how to rename the three background images, which are distributed under
distinct filenames because they share a filename in the directory layout.

### (b) Moor et al. 2018 — villus zonation reference

From Moor, A. E., Harnik, Y., Ben-Moshe, S., Massasa, E. E., Rozenberg, M., Eilam, R.,
... & Itzkovitz, S. (2018). *Spatial reconstruction of single enterocytes uncovers broad
zonation along the intestinal villus axis.* Cell, 175(4), 1156-1167:

- `table_A_LCM_TPM_values.tsv`
- `table_D_zonation_reconstruction.tsv`

Save both into a directory named `in_vivo_villus_data`.

### (c) Mouse Visium HD — GEO accession GSE303705

Save into a directory named `mouse_visium`, giving
`mouse_visium/GSE303705_RAW/` and `mouse_visium/rep2/day0/` etc.

### Resulting layout

```
<data directory>/
├── raw/                     # monolayer transcripts, segmentation, cell-by-gene tables
├── backgrounds/
│   ├── pasadena_roi1/hyb_background_aligned.tiff
│   ├── nov23_72hr_roi1/hyb_background_aligned.tiff
│   └── nov23_12hr_roi2/hyb_background_aligned.tiff
├── sprinkling_nov_23/       # cell transplantation
├── perturbations/           # pharmacological perturbations
├── unperturbed/
│   └── monolayer_erosion/
├── sprinkled/
│   ├── 12hr/{roi1,roi2,roi3}
│   └── 72hr/{roi1,roi2,roi3,roi4}
├── in_vivo_villus_data/     # from (b), downloaded separately
└── mouse_visium/            # from (c), downloaded separately
```

The background images in the Zenodo record keep only the channels the figures read
(acquisition channel 3 for `pasadena_roi1`; channels 2 and 3 for the two nov23 images).
The plotting code indexes them at their stored positions. The full multi-channel
acquisitions are archived separately.

## 4. Point the code at the data

Either set an environment variable:

```bash
export ZONATION_DATA_DIR=/path/to/your/data/directory
```

or edit `HOME_DIR` in `utils/constant.py`. Every other data path derives from it. If the
directory is missing, importing `utils.constant` fails immediately with a message saying
so, rather than failing later inside a plotting function.

## 5. Generate the figures

```python
from paper.plotScripts.plotALL import plot_all_figures
plot_all_figures()
```

Each figure runs independently, so one that fails does not prevent the others from
being produced; the run ends with a summary of which succeeded and how many files each
wrote.

Output goes to `paper/graphs/`:

| directory | figure |
|---|---|
| `autonomous_zonation_figure_plots/` | Figure 1 — autonomous zonation |
| `scale_invariance_plots/` | Figure 2 — scale invariance |
| `zonation_plasticity_plots/` | Figure 3 — cell transplantation |
| `neighborhood_zone_adoption_plots/` | Figure 4 — zone confusion |
| `pharmacological_perturbations/` | Figure 5 — pharmacological perturbations |
| `continuous_regenerative_figure_plots/` | Figure 6 — continuous regenerative response |
