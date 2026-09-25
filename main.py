"""Entry point: regenerate the paper figures.

Set the data location first — either export ZONATION_DATA_DIR, or edit HOME_DIR in
utils/constant.py. See README.md.
"""

from paper.plotScripts.plotALL import *

if __name__ == '__main__':
    set_style()

    # Every figure. Each one runs independently, so a failure in one does not stop
    # the rest; a summary at the end reports what succeeded and how many files it wrote.
    plot_all_figures()

    # Or run a single figure by uncommenting one of these instead:
    # plot_all_autonomous_figure_plots(saved_datasets=True)   # Figure 1  autonomous zonation
    # plot_scale_invariance_figures()                         # Figure 2  scale invariance
    # plot_zonation_plasticity_plots(calculate=False)         # Figure 3  cell transplantation
    # plot_all_neighborhood_zone_adoption_plots()             # Figure 4  zone confusion
    # plot_all_pharmacological_perturbation_plots()           # Figure 5  perturbations
    # plot_all_continuous_regenerative_response_plots()       # Figure 6  regenerative response
