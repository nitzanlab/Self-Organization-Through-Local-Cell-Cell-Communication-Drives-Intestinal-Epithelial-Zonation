
from paper.plotScripts.autonomous_zonation_plots import *
from paper.plotScripts.zonation_plasticity import *
from paper.plotScripts.neighborhood_zone_adoption import *
from paper.plotScripts.continuous_regenerative_response import *
from paper.plotScripts.scale_invariance_profiles import plot_scale_invariance_figures
from paper.plotScripts.pharmacological_perturbations import plot_all_pharmacological_perturbation_plots

def plot_all_figures():
    plot_all_autonomous_figure_plots(saved_datasets=False) #first time running this function, saved_datasets
    #must be False because they have not been calculated yet. After the first run, saved datasets
    #can be changed to True and thus use the caclulated results to reproduce plots
    plot_zonation_plasticity_plots(calculate=True) #first time running it, the function must have calculate=True,
    #but for the next runs, can run with calculate=False, since the necessary calculations have been done and saved
    #and thus, can be accessed and used instead of recalculated
    plot_all_neighborhood_zone_adoption_plots()
    plot_all_continuous_regenerative_response_plots()
    plot_scale_invariance_figures()  # Figure 2: scale invariance (Visium, rep2 sub-ROI)
    plot_all_pharmacological_perturbation_plots()  # Figure 4: pharmacological perturbations
