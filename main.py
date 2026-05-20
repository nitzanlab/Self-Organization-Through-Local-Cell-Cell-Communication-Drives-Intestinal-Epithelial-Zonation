
from paper.plotScripts.plotALL import *

from paper.plotScripts.continuous_regenerative_response import *
from paper.plotScripts.continuous_regenerative_response import *
from paper.plotScripts.erosion_calculations import compare_profiles_16_61_um

if __name__ == '__main__':
   set_style()
   #run all analyses and plot:
   #plot_all_figures()
   # plot_zonation_plasticity_plots(calculate=False)
   # plot_all_autonomous_figure_plots(saved_datasets=True)
   # calculate_eroded_transcription_densities(plot_erosion=True, save_components=True)
   # plot_all_neighborhood_zone_adoption_plots()
   # plot_all_continuous_regenerative_response_plots()
plot_all_autonomous_figure_plots(saved_datasets=True)
# plot_zonation_plasticity_plots(calculate=False)
# plot_all_pharmacological_perturbation_plots()
# plot_all_neighborhood_zone_adoption_plots()
# plot_all_continuous_regenerative_response_plots()
# plot_scale_invariance_figures()
# unweighted (your current request)
   # compare_profiles_16_61_um(genes=ORGANOID_GENE_NAMES, do_weighted=False)
#    compare_profiles_plotlogic_vs_invivo(
#     genes=ORGANOID_GENE_NAMES,  # or ORGANOID_GENE_NAMES_NOGFP
#     lo_um=16, hi_um=61,
#     normalize=True, smooth=True, do_weighted=False
# )
# 1) Build correlations & overlays (already done earlier)
# compare_monolayer_profiles_to_invivo(...)  # make sure do_weighted=True if you want weighted scoring
# compare_monolayer_profiles_to_invivo(
#     result_dict=load_unperturbed_monolayer_gene_densities(),
#     genes=ORGANOID_GENE_NAMES,   # or ORGANOID_GENE_NAMES_NOGFP
#     LO=16, HI=61,
#     normalize=True,
#     smooth=True,
#     do_weighted=True,            # <— enable weighted so it’s available
#     out_dir=None,
#     save_csv=True
# )

# 2) Summarize enterocytes at threshold 0.5
# summarize_enterocyte_correlations(
#     corr_csv=None,                              # will pick the default path
#     enterocyte_csv=os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, "sorted_indexing.csv"),
#     threshold=0.5,
#     out_dir=None,
#     include_celltypes=("enterocyte",),          # adjust or extend if needed
#     corr_col_priority=("pearson_weighted", "pearson", "spearman")
# )


# # when you’re ready to consider SD in the metric:
#    compare_profiles_16_61_um(genes=ORGANOID_GENE_NAMES, do_weighted=True)
