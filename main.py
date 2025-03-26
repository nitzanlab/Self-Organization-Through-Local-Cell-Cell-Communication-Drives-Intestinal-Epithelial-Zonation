
from paper.plotScripts.plotALL import *

from paper.plotScripts.continuous_regenerative_response import *


if __name__ == '__main__':
   set_style()
   #run all analyses and plot:
   #plot_all_figures()
   plot_zonation_plasticity_plots(calculate=False)
   #plot_all_autonomous_figure_plots(saved_datasets=True)
   #plot_all_neighborhood_zone_adoption_plots()