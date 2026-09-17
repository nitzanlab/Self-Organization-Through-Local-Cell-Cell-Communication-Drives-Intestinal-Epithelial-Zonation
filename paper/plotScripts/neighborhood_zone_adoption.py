import matplotlib.pyplot as plt
import pandas as pd

from utils.imports import *
from utils.constant import *
from paper.extractedData.load_csvs import *

def plot_all_neighborhood_zone_adoption_plots():
    ###panel a: schematic diagram created in  https://BioRender.com
    ###panel b: expected cell zones in unperturbed monolayer
    #plot_expected_cell_zones()

    ###panel c: expected zones in inserted cells, subregion of inserted monolayers
    # plot_inserted_cells_expected_zones_spatially()

    ###panel d: correlation in expected zone to neighbors
    #plot_expected_zone_correlations(ORGANOID_GENE_NAMES_NOGFP, num_neigh=5)

    #plot_expected_zone_correlations(ZONE_MAPPING_GENES, num_neigh=5)

    ###panel e: expected zone distribution
    # plot_inserted_cells_expected_zone_distribution()
    #
    # ###panel f: inserted cells zone confusion
    plot_inserted_cells_zone_confusion_distribution('12hr')
    plot_inserted_cells_zone_confusion_distribution('72hr')
    # #
    # # ###panel g: gene contribution to zone confusion
    # plot_zone_confusion_gene_contribution('72hr', genes=ZONE_MAPPING_GENES)


def plot_zone_confusion_gene_contribution(tmpt, genes=ZONE_MAPPING_GENES):
    """
    This function plots the contribution of each of the zone mapping genes to the zone probability of cells
    in the same zone confusion range.
    :param tmpt: the post insertion time when SEQfish was used to measure the gene expression of the monolayer.
    :param genes: the genes used to map the cells to zones.
    :return:
    """
    all_confusion_heatmaps = []
    for entropy_range in ENTROPY_RANGES:
        gene_contribution = calculate_gene_zone_confusion_contribution(tmpt, genes, 3, 14, True,
                                                                       entropy_range)
        all_confusion_heatmaps.append(gene_contribution)
    vmin = min(d.min() for d in all_confusion_heatmaps)
    vmax = max(d.max() for d in all_confusion_heatmaps)

    quadmesh = None
    # Create figure and subplots
    fig, axes = plt.subplots(1, 4, figsize=(8, 3), sharex=True, gridspec_kw={'wspace': 0.05})
    fig.suptitle("Gene Zone Confusion Contribution", fontweight="bold")
    for i, ax in enumerate(axes): # Only add colorbar to the last subplot
        hm = sns.heatmap(all_confusion_heatmaps[i], ax=ax, vmin=vmin, xticklabels=genes,vmax=vmax, cmap="plasma", cbar=False)
        quadmesh = hm
        ax.set_title(f"{ENTROPY_RANGES[i][0]}-{ENTROPY_RANGES[i][1]}")  # Set title for each heatmap
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontweight="bold")
        ax.set_xlabel("Genes")  # Label only x-axis (it will be shared)
        ax.set_ylabel("")  # Remove y-axis labels for aesthetics
        ax.set_yticks([])

        if i==0:
            ax.set_ylabel('Cells')
        else:
            ax.set_ylabel("")
        ax.set_xticks# Remove y-axis labels for other heatmaps
    # Create a separate colorbar on the right side
    cbar_ax = fig.add_axes([0.92, 0.35, 0.02, 0.5])  # [x-position, y-position, width, height]
    cbar = fig.colorbar(quadmesh.get_children()[0], cax=cbar_ax)  # Use get_children()[0] to access the heatmap colors
    cbar.ax.set_ylabel("Zone Gene Contribution", fontweight="bold")  # Add colorbar title

    # Adjust layout to prevent cutting off xticks
    plt.subplots_adjust(left=0.05, right=0.9, top=0.85, bottom=0.35)  # Adjust right space for colorbar
    os.makedirs(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, 'gene_zone_contribution.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()
    plt.close()


def calculate_gene_zone_confusion_contribution(tmpt, genes, begin, end, is_sprinkled, entropy_range, num_cells=10):
    """
    This function calculates the gene confusion contribution given cells in the same zone confusion range.
    This allows viewing the impact of different cells in determine cell zone
    :param tmpt: time following insertion in which the gene expression across the monolayer was measured
    :param genes: genes used for mapping cell zone
    :param begin: erosion iteration from which expression profiles were measured
    :param end: erosion iteration until which expression profiles were measured
    :param is_sprinkled: if the cells for which we view the gene zone confusion contribution are inserted or not
    :param entropy_range: the zone confusion range
    :param num_cells: number of cells to present in the heatmap
    :return: an array of shape num cells X genes showing how much each cell expresses each gene such in a way that
    """
    adata = get_tmpt_all_rois_adata(tmpt)
    adata = adata[adata[:, ZONE_MAPPING_GENES].X.sum(axis=1) > 0.1]
    adata = map_monolayer_to_transcript_density_profiles(adata, genes, begin, end)
    adata = get_transcript_density_entropy(adata)
    title = ''
    if is_sprinkled:
        adata = adata[adata.obs['spc']]
        title = 'GFP'
    else:
        adata = adata[~adata.obs['spc']]
        title = 'non GFP'
    adata_ent_range =  adata[(adata.obs['zone_entropy']>=entropy_range[0]) & (adata.obs['zone_entropy']<=entropy_range[1])]
    random_cells_idx = np.random.choice(adata_ent_range.shape[0], size=num_cells, replace=False)
    adata_ent_range_subset = adata_ent_range[random_cells_idx,:]
    selected_data = adata_ent_range_subset.obsm['transcript_zone_dist']
    sorted_indices = np.argsort(selected_data[:, 0]) #probability for zone 0

    adata_ent_range_subset = adata_ent_range[random_cells_idx, :]

    gene_contribution = adata_ent_range_subset[:, genes].X[sorted_indices,:]
    normalized_gene_contirbution = gene_contribution / (gene_contribution.sum(axis=1)[:, np.newaxis]+1e-3)
    return normalized_gene_contirbution


def plot_expected_cell_zones():
    """
    This function plots the cell in the unperturbed monolayer based on their expected zone -
    based on how they match the zones measured by eroding the monolayer from its edge inwards.
    """
    plot_expected_position_mapping_unperturbed_cosine_sim(ZONE_MAPPING_GENES, 'zone mapping genes', to_plot=True,
                                                          to_save=True)


def plot_expected_position_mapping_unperturbed_cosine_sim(genes, gene_title='',zoned=False, x_range=None, y_range=None, to_plot=True, to_save=False, save_name=''):
    """
    This function measures the expected zone of the cells in  the unperturbed monolayer and plots the cells spatially,
    colored based on their expected zone.
    to the transcript density profiles measured.
    :param genes:  Genes used for mapping
    :param gene_title: the name of the gene group
    :param zoned: if the plot is of a subregion of the monolayer
    :param x_range: if zoned, then the region plotted will be cells in the given x range
    :param y_range: if zoned, then the region plotted will be cells in the given y range
    :param to_plot: if to plot the expected zone of the cells spatially
    :param to_save: if to save the cells with their expected zones and x,y coordinates as a pandas dataframe
    :param save_name:  name of the dataframe
    """

    #load the unperturbed monolayer
    adata = get_unperturbed_monolayer_adata()

    #map the cells in the unperturbed monolayer to the transcript density profiles based on the given genes
    adata = map_monolayer_to_transcript_density_profiles(adata, genes)
    s=6

    if zoned: #subset of cells that are in the given zone
        aoi = (adata.obsm[COORDINATES][X_COORDINATES] > x_range[0]) & (
                adata.obsm[COORDINATES][X_COORDINATES] < x_range[1]) & (
                          adata.obsm[COORDINATES][Y_COORDINATES] > y_range[0]) & (
                          adata.obsm[COORDINATES][Y_COORDINATES] < y_range[1])

        adata = adata[aoi]
        s=10

    vmin = adata.obs[EXPECTED_ZONE_TB].min()
    vmax = adata.obs[EXPECTED_ZONE_TB].max()
    sc1 = plt.scatter(adata.obsm[COORDINATES][X_COORDINATES], -adata.obsm[COORDINATES][Y_COORDINATES], cmap='viridis',
                      c=adata.obs[EXPECTED_ZONE_TB], vmin=vmin, vmax=vmax, s=s)
    cbar = plt.colorbar(sc1, label='Expected Zone')
    plt.title(f'Unperturbed Monolayer {gene_title}\n Transcript-Based Zone Mapping ')
    plt.xticks([])
    plt.yticks([])
    if to_plot:
        plt.show()
    if to_save:
        save_spatial_signal(adata, EXPECTED_ZONE_TB, f'{gene_title}_expected_zones_unperturbed_monolayer')


def plot_inserted_cells_expected_zones_spatially():
    """
    This function plots the panels showcasing the expected zones of inserted cells in comparison to their
    neighboring non inserted cells' expected zones. We show an example of a subregion of the monolayer measured 12 hours following
    cell insertion and another example of the monolayer measured 72 hours following cell insertion
    """
    adata_GFP_12hr_roi1 = load_sprinkled_adata_hr_tmpt('12hr', 'roi1')
    adata_GFP_12hr_roi1 = map_monolayer_to_transcript_density_profiles(adata_GFP_12hr_roi1, ZONE_MAPPING_GENES)
    plot_GFP_adata_signal_spatially(adata_GFP_12hr_roi1, 'transcript_exp_pos', '12hr roi1 GFP', zoned=True, x_range=GFP_12HR_ROI1_X,
                                    y_range=GFP_12HR_ROI1_Y, title='12hr_roi1_GFP_expected_zone')

    adata_GFP_72hr_roi2 = load_sprinkled_adata_hr_tmpt('72hr', 'roi2')
    adata_GFP_72hr_roi2 = map_monolayer_to_transcript_density_profiles(adata_GFP_72hr_roi2, ZONE_MAPPING_GENES)
    plot_GFP_adata_signal_spatially(adata_GFP_72hr_roi2, 'transcript_exp_pos', '72hr roi2 GFP', zoned=True,
                                    x_range=GFP_72HR_ROI2_X,
                                    y_range=GFP_72HR_ROI2_Y, title='72hr_roi2_GFP_expected_zone')




def plot_expected_zone_correlations(genes, num_neigh):
    """
    This function plots a bar plot showing the correlation in zone of each cell and the mean zone of its nearest neighbors.
    Each bar shows the mean +- 1 stadnard deviation of the correlation for each cell group  - unperturbed monolayer
    cells, inserted cells to their neighbors (for both timepoints), non-inserted cells to their neighbors (for both timepoints)
    :param genes: genes used to calculate cell zones
    :param num_neigh: number of neighbors takes into account for measuring the mean zone of the cell's neighbors
    """
    #plot_all_expected_position_by_transcript_corr_with_GFP_neighbors(genes):
    wt_env_zone_corr = get_env_transcript_density_zonation_correlation_in_wt_monolayer(genes, num_neigh)
    corrs_GFP = {'72hr': [], '12hr': []}
    corrs_non_GFP = {'72hr': [], '12hr': []}
    mean_GFP_corrs_tmpts = []
    std_GFP_corrs_tmpts = []
    mean_nonGFP_corrs_tmpts = []
    std_nonGFP_corrs_tmpts = []
    for tmpt in TMPT_TO_ROIS_DICT.keys():
        for roi in TMPT_TO_ROIS_DICT[tmpt]:
            GFP_corr, nonGFP_corr = get_neighbor_corrs(tmpt, roi, genes, num_neigh)
            if GFP_corr is None:
                continue
            corrs_GFP[tmpt].append(GFP_corr)
            corrs_non_GFP[tmpt].append(nonGFP_corr)
        mean_nonGFP_corrs_tmpts.append(np.mean(corrs_non_GFP[tmpt]))
        mean_GFP_corrs_tmpts.append(np.mean(corrs_GFP[tmpt]))
        std_nonGFP_corrs_tmpts.append(np.std(corrs_non_GFP[tmpt]))
        std_GFP_corrs_tmpts.append(np.std(corrs_GFP[tmpt]))
    x_pos = np.arange(5)
    means = np.concatenate((mean_GFP_corrs_tmpts, mean_nonGFP_corrs_tmpts, [wt_env_zone_corr]))
    errors = np.concatenate((std_GFP_corrs_tmpts, std_nonGFP_corrs_tmpts, [0]))
    plt.figure(figsize=(8, 6))
    bar_colors = '#9B30FF'  # Varying shades of purple

    # Plotting
    bars = plt.bar(x_pos, means, yerr=errors, color=bar_colors, capsize=5, width=0.6, edgecolor='black',
                   linewidth=0.7)

    # Labels and formatting
    plt.xticks(x_pos, ['12hr\n GFP', '72hr\n GFP', '12hr\nnon-GFP', '72hr\nnon-GFP', 'Unperturbed'],
               fontsize=16)
    plt.xlabel('Cell Group Type', fontsize=20)
    plt.ylabel('Correlation', fontsize=20)
    plt.title('Expected Zone Correlation', fontsize=22)
    plt.tight_layout(pad=4)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    os.makedirs(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, 'zone_neighborhood_correlation_bar_plot.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()
    plt.close()

def get_neighbor_corrs(tmpt, roi, genes, num_neigh, radius=500):
    adata = load_sprinkled_adata_hr_tmpt(tmpt, roi)
    adata = map_monolayer_to_transcript_density_profiles(adata, genes)
    #find k nearest neighbors in radius for all cells
    tree = BallTree(adata.obsm['spatial'].values)

    indices, distances = tree.query_radius(adata.obsm['spatial'].values, r=radius, return_distance=True)
    enough_indices_idxs = [i for i in range(adata.shape[0]) if len(indices[i]) >= num_neigh]
    neighbors_of_enough_cells = [indices[i][:num_neigh] for i in range(adata.shape[0]) if len(indices[i])>=num_neigh]
    adata_neighs = adata[enough_indices_idxs,:]
    adata_neighs.obs['neighbor_mean_transcript_zone'] = np.array(adata.obs['transcript_exp_pos'])[neighbors_of_enough_cells].mean(axis=1)
    adata_neigh_spc = adata_neighs[adata_neighs.obs['spc']]
    adata_neigh_non_spc = adata_neighs[~adata_neighs.obs['spc']]
    GFP_zone_corr, GFP_pval = pearsonr(adata_neigh_spc.obs['transcript_exp_pos'],  adata_neigh_spc.obs['neighbor_mean_transcript_zone'])
    non_GFP_zone_corr, non_GFP_pval = pearsonr(adata_neigh_non_spc.obs['transcript_exp_pos'], adata_neigh_non_spc.obs['neighbor_mean_transcript_zone'])
    print( GFP_zone_corr, GFP_pval ,f'GFP {tmpt} {roi}')
    print( non_GFP_zone_corr, non_GFP_pval , f'non GFP {tmpt} {roi}')
    return GFP_zone_corr, non_GFP_zone_corr

def get_env_transcript_density_zonation_correlation_in_wt_monolayer(genes, num_neigh):
    wt_adata = get_unperturbed_monolayer_adata(num_neigh)
    wt_adata = map_monolayer_to_transcript_density_profiles(wt_adata, genes, 3, 14)
    wt_neighbors = np.array(wt_adata.obs['transcript_exp_pos'])[np.array(wt_adata.obsm['neighbors_idx'])]
    correlation_in_exp_pos, pval = stats.pearsonr(wt_adata.obs['transcript_exp_pos'],
                                                  wt_neighbors.mean(axis=1))
    print(f'corr: {correlation_in_exp_pos} pval:{pval:.3f}')
    return correlation_in_exp_pos

def plot_inserted_cells_expected_zone_distribution():
    """
    This function plots the expected zone distribution as a histogram for each timepoitn
    """
    for tmpt in TMPT_TO_ROIS_DICT.keys():
        expected_zone_distribution(tmpt, ZONE_MAPPING_GENES,3, 14)

def expected_zone_distribution(tmpt, genes, begin, end):
    """
    This function plots the distirbution of the expected zones of the cells comparing non-inserted cells to inserted cells
    in all of the regions of interest in the same monolayer
    :param tmpt: time in hours following cell insertion
    :param genes: genes used to map cells probabilistically to zones
    :param begin: erosion step from which the transcript density profiles were used for mapping cell zone
    :param end:  erosion step until which the transcipt density profiles were used for mapping cell zone
    :return:
    """
    adata = get_tmpt_all_rois_adata(tmpt)
    adata = map_monolayer_to_transcript_density_profiles(adata, genes, begin, end, binned=True)
    adata = get_transcript_density_position_std(adata)
    bins = np.linspace(0, 5, 20)
    adata_non_spc = adata[~adata.obs['spc']]
    adata_spc = adata[adata.obs['spc']]
    plt.hist(adata_non_spc.obs['transcript_exp_pos'], density=True, bins=bins, label='non GFP', alpha=0.7, color='grey')
    plt.hist(adata_spc.obs['transcript_exp_pos'], density=True, bins=bins, label='GFP', alpha=0.7, color='green')
    plt.title(f'Expected Zone Distribution {tmpt}')
    stat, p_value = ks_2samp(adata_non_spc.obs['transcript_exp_pos'], adata_spc.obs['transcript_exp_pos'])
    print(f"KS Statistic: {stat}, p-value: {p_value}")
    plt.xlabel('Expected Zone')
    plt.ylabel('Relative Density')
    plt.ylim(0,1.2)
    plt.legend(loc='upper left')
    plt.tight_layout()
    os.makedirs(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, f'expected_zone_distribution_{tmpt}.pdf')
    plt.savefig(file_name, format='pdf',bbox_inches='tight')
    plt.show()
    plt.close()

def get_transcript_density_entropy(adata):
    """
    This function adds an obs to the anndata given, adata, of the zone confusion of each cell,
    which is measured as the entropy of the zone probability distribution.
    """
    zone_entropy = np.apply_along_axis(lambda x: entropy(x, base=2), axis=1, arr=adata.obsm['transcript_zone_dist'])
    adata.obs['zone_entropy'] = zone_entropy
    return adata
def get_transcript_density_position_std(adata):
    """

    :param adata:
    :return:
    """
    exp_pos = np.array(adata.obs['transcript_exp_pos'])
    pos_dist = adata.obsm['transcript_zone_dist']

    squared_diff = (np.tile(np.arange(pos_dist.shape[1]), (len(exp_pos), 1)) - exp_pos[:,
                                                                               np.newaxis]) ** 2  # Shape: (n_samples, n_positions)
    weighted_squared_diff = pos_dist * squared_diff
    variance = np.sum(weighted_squared_diff, axis=1)
    pos_std = np.sqrt(variance)
    adata.obs['transcript_pos_std'] = pos_std
    return adata

def plot_inserted_cells_zone_confusion_distribution(tmpt):
    """
    This function plots the cell zone confusion distribution for all cells measured in the same monolayer, tmpt hours after
    cell insertion. The histogram compares the zone confusion distribution of inserted cells vs. non inserted cells
    """
    adata = get_tmpt_all_rois_adata(tmpt)
    adata = adata[adata[:, ZONE_MAPPING_GENES].X.sum(axis=1) > 0.1]
    adata = map_monolayer_to_transcript_density_profiles(adata, ZONE_MAPPING_GENES, 3, 14)
    adata = get_transcript_density_entropy(adata)

    adata_gfp = adata[adata.obs['spc']]
    adata_non_gfp = adata[~adata.obs['spc']]
    bins = np.linspace(1.6, 2.6, 20)
    plt.figure(figsize=(8, 6))      # own figure: never draw onto a leaked one
    plt.hist(adata_non_gfp.obs['zone_entropy'], alpha=0.7, label='non GFP', density=True, color='grey', bins=bins)
    plt.hist(adata_gfp.obs['zone_entropy'], alpha=0.7, label='GFP', density=True, color='green', bins=bins)
    plt.xlabel('Zone Entropy')
    plt.ylabel('Density')
    plt.ylim(0, 6)
    plt.title(f'{tmpt} Zone Entropy')
    plt.legend(loc='upper left')
    os.makedirs(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH, f'zone_confusion_distribution_{tmpt}.pdf')
    plt.savefig(file_name, format='pdf', bbox_inches='tight')
    plt.tight_layout()
    plt.show()
    plt.close()
        # ---- Summary statistics: mean confusion for GFP vs non-GFP ----
    mean_entropy_gfp = float(adata_gfp.obs['zone_entropy'].mean())
    std_entropy_gfp = float(adata_gfp.obs['zone_entropy'].std())
    n_gfp = int(adata_gfp.shape[0])

    mean_entropy_non = float(adata_non_gfp.obs['zone_entropy'].mean())
    std_entropy_non = float(adata_non_gfp.obs['zone_entropy'].std())
    n_non = int(adata_non_gfp.shape[0])
    stat, p_value = ks_2samp(adata_gfp.obs['transcript_exp_pos'], adata_non_gfp.obs['transcript_exp_pos'])
    print(f"[{tmpt}] GFP cells: mean zone entropy = {mean_entropy_gfp:.3f} ± {std_entropy_gfp:.3f} (n={n_gfp})")
    print(f"[{tmpt}] non-GFP cells: mean zone entropy = {mean_entropy_non:.3f} ± {std_entropy_non:.3f} (n={n_non})")
    print(f"KS Statistic: {stat}, p-value: {p_value} {tmpt}")

def plot_GFP_adata_signal_spatially(adata, signal_name, adata_type, zoned=False, x_range=None, y_range=None, to_plot=True, title=''):
    """
    This function plots a signal over the cells spatially in monolayers with inserted cells, labeling the inserted
    vs not inserted cells.
    :param adata: the anndata containing the cell information
    :param signal_name: the signal to be plotted spatially , which is saved in adata.obs
    :param adata_type: the meta data of the adata - tmpt , roi,..
    :param zoned: boolean: if to plot a subregion of the monolyaer
    :param x_range: if plotting a subregion of the monolyaer, this is the x range of the plot
    :param y_range: if plotting a subregion of the monolayer, this is the y range of the plot
    :param to_plot: boolean: if to show the plot
    :param title: the title of the plot
    """
    adata_spc = adata[adata.obs['spc']]
    adata_non_spc = adata[~adata.obs['spc']]
    plt.figure(figsize=(8, 6))
    s = 30
    if zoned:
        plt.figure(figsize=(6, 4))
        aoi_spc = (adata_spc.obsm[COORDINATES][X_COORDINATES] > x_range[0]) & (
                adata_spc.obsm[COORDINATES][X_COORDINATES] < x_range[1]) & (
                          adata_spc.obsm[COORDINATES][Y_COORDINATES] > y_range[0]) & (
                          adata_spc.obsm[COORDINATES][Y_COORDINATES] < y_range[1])
        aoi_non_spc = (adata_non_spc.obsm[COORDINATES][X_COORDINATES] > x_range[0]) & (
                adata_non_spc.obsm[COORDINATES][X_COORDINATES] < x_range[1]) & (
                          adata_non_spc.obsm[COORDINATES][Y_COORDINATES] > y_range[0]) & (
                          adata_non_spc.obsm[COORDINATES][Y_COORDINATES] < y_range[1])
        adata_spc = adata_spc[aoi_spc]
        adata_non_spc = adata_non_spc[aoi_non_spc]

        aoi = (adata.obsm[COORDINATES][X_COORDINATES] > x_range[0]) & (
                adata.obsm[COORDINATES][X_COORDINATES] < x_range[1]) & (
                          adata.obsm[COORDINATES][Y_COORDINATES] > y_range[0]) & (
                          adata.obsm[COORDINATES][Y_COORDINATES] < y_range[1])
        adata = adata[aoi]
        save_spatial_signal(adata, signal_name, f'zoom_in_{title}', sprinkled=True)


    vmin = min(adata_spc.obs[signal_name].min(), adata_non_spc.obs[signal_name].min())
    vmax = max(adata_spc.obs[signal_name].max(), adata_non_spc.obs[signal_name].max())

    sc1 = plt.scatter(adata_non_spc.obsm[COORDINATES][X_COORDINATES], adata_non_spc.obsm[COORDINATES][Y_COORDINATES], cmap='plasma',
                      c=adata_non_spc.obs[signal_name], vmin=vmin, vmax=vmax, s=s)
    plt.scatter(adata_spc.obsm[COORDINATES][X_COORDINATES], adata_spc.obsm[COORDINATES][Y_COORDINATES], label='GFP',
                edgecolors='lime', cmap='plasma', c=adata_spc.obs[signal_name], vmin=0, vmax=4, s=s)

    cbar = plt.colorbar(sc1, label='Expected Zone')
    plt.legend(loc='upper left', fontsize='large')
    plt.title(f'{adata_type} monolayer {signal_name}\n ')
    plt.xticks([])
    plt.yticks([])
    plt.show()
    plt.close()


##mapping based on the in vivo expression profiles
def map_monolayer_to_invivo_expression_profiles(adata, genes):
    """
    This function maps the cell in the monolayer probabilistically to zones defined by the in vivo villus data collected
    in Moor et al 2018 Cell , 'Spatial Reconstruction of Single Enterocytes Uncovers Broad Zonation along the Intestinal
    Villus axis
    :param adata: anndata of the monolayer to be mapped
    :param genes: the genes used for mapping the cells to zones
    :return:
    """
    invivo_expression = load_TPM_LCM_intestine_atlas()
    invivo_exp_profile, gene_exp_list = preprocess_LCM_atlas_only_core_reference_genes(invivo_expression, genes)
    genes_paneled_from_gene_list = [gene for gene in gene_exp_list if gene in adata.var_names]
    adata_subset = adata[:, genes_paneled_from_gene_list].copy()
    cosine_sim_per_position = cosine_similarity(adata_subset.X, invivo_exp_profile.T)
    row_sums = np.sum(cosine_sim_per_position, axis=1)
    zero_sum_rows = (row_sums==0)
    cosine_sim_per_position[~zero_sum_rows] = cosine_sim_per_position[~zero_sum_rows] / row_sums[~zero_sum_rows, np.newaxis]
    cosine_sim_per_position[zero_sum_rows] = 1/cosine_sim_per_position.shape[1]
    max_pos = np.argmax(cosine_sim_per_position, axis=1)
    adata.obs['invivo_max_pos'] = max_pos
    adata.obs['invivo_exp_pos'] = cosine_sim_per_position @ np.arange(
        cosine_sim_per_position.shape[1])
    adata.obsm['invivo_position_dist'] = cosine_sim_per_position
    #adata.obs['invivo_position_dist_std'] = np.std(adata.obsm['invivo_position_dist'], axis=1)
    return adata