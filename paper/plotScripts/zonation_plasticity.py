
from utils.imports import *
from utils.constant import *
from paper.extractedData.load_csvs import *
def plot_zonation_plasticity_plots(calculate = True):
    ###necessary calculations needed to perform, and save for analyses and plots
    if calculate:
        calculate_unperturbed_monolayer_neighborhood_gene_expression_correlations(save=True)
        calculate_unperturbed_monolayer_morans_i(save=True)

    ###panel a: schematic diagram created in BioRender.com
    ###panel b and c: raw image example #TODO Yael
    ###panel d and e : ### correlation comparison #TODO Yael

    ###panel f: gene expression neighborhood correlation in inserted cells
    #plot_gene_expression_neighborhood_correlation_in_inserted_cells()
    plot_gene_correlation_histograms_GFP_to_wt(len(ORGANOID_GENE_NAMES_NOGFP),'all_paneled_genes')
    plot_gene_correlation_histograms_GFP_to_wt(20,
                                            'top_20_moransi_genes')

    ###panel g: inserted cell neighborhood correlation across timepoints  #TODO Yael


def calculate_unperturbed_monolayer_morans_i(save=True):
    pass

def calculate_unperturbed_monolayer_neighborhood_gene_expression_correlations(save=True):
    pass

def plot_gene_correlation_histograms_GFP_to_wt(num_genes, save_title):
    """
    This function plots the distribution of correlation in gene expression between a cell and its neighboring cells in the following
    three settings. 1) in the unperturbed monolayer 2) between inserted cells for 12hr and their non inserted neighboring cells , averaged across all rois,
    3) like 2) but for cells that have been inserted for 72hr hours at time of measuring expression
    :param num_genes: the top n genes in regards to the gene's moran's I value in the unperturbed monolayer
    :return:
    """
    morans_i_wt = load_unperturbed_monolayer_genes_morans_i()

    #get the name of then genes with the highest moran's I in the unperturbed monolayer
    top_n_moran_i_genes = np.array(morans_i_wt.sort_values(by='morans_i', ascending=False).head(num_genes).index)

    # get top n morans i genes index in from all non GFP gnes : ORGANOID_GENE_NAMES_NO_GFP
    moran_i_genes_idx = [np.where(ORGANOID_GENE_NAMES_NOGFP == item)[0][0] for item in top_n_moran_i_genes]

    #load unperturbed monoalyer gene expression correlations:
    wt_non_GFP_gene_corr = pd.read_csv(os.path.join(WT_MONOLAYER_DIR, 'wt_gene_correlation.csv'), index_col=0)
    wt_non_GFP_gene_corr = wt_non_GFP_gene_corr.reindex(ORGANOID_GENE_NAMES_NOGFP)
    wt_non_GFP_gene_corr = wt_non_GFP_gene_corr['corr'].loc[top_n_moran_i_genes]

    #get correlation in expression between inserted cells and their neighboring cells in each time point
    GFP_72hr_gene_corr = get_all_GFP_to_env_one_tmpt_gene_exp_corr('72hr').iloc[moran_i_genes_idx]
    GFP_12hr_gene_corr = get_all_GFP_to_env_one_tmpt_gene_exp_corr('12hr').iloc[moran_i_genes_idx]

    bins = np.linspace(-0.25, 0.8, 20)
    plt.hist(wt_non_GFP_gene_corr, bins=bins, color='purple', alpha=0.8)
    plt.hist(GFP_72hr_gene_corr, bins=bins, color='orange', alpha=0.8)
    plt.hist(GFP_12hr_gene_corr, bins=bins, color='green', alpha=0.8)

    plt.axvline(np.mean(GFP_72hr_gene_corr), label='mean inserted 72hr ', color='orange', linestyle='--', linewidth=2)
    plt.axvline(np.mean(GFP_12hr_gene_corr), label='mean inserted 12hr', color='green', linestyle='--', linewidth=2)
    plt.axvline(np.mean(wt_non_GFP_gene_corr), label='mean unperturbed ', color='purple', linestyle='--', linewidth=2)
    if num_genes == len(ORGANOID_GENE_NAMES_NOGFP):
        plt.legend()
        plt.xlabel('Correlation')
        plt.ylabel('Frequency')
        plt.title(f'Neighboring Cells Gene Expression Correlation')
    plt.tight_layout()
    os.makedirs(ZONATION_PLASTICITY_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(ZONATION_PLASTICITY_PLOTS_FOLDER_PATH,
                             f'{save_title}_neighborhood_expression_correlations_all_genes.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()

def get_all_GFP_to_env_one_tmpt_gene_exp_corr(tmpt):
    """
    This function measures the correlation in expression between inserted cells and their non-inserted neighboring cells per gene, averaged across
    the rois
    :param tmpt: the timepoint - the amount fo hours the cells were inserted for before measuring gene expression
    """
    gene_corrs = {}
    for roi in TMPT_TO_ROIS_DICT[tmpt]:
        gene_corrs[roi] = get_GFP_to_env_gene_exp_corr_one_roi(tmpt, roi)
    mean_corrs_per_tmpt = pd.DataFrame(gene_corrs).mean(axis=1)
    return mean_corrs_per_tmpt


def get_GFP_to_env_gene_exp_corr_one_roi(tmpt, roi, radius=500):
    """
    This function measures the correlation in expression between inserted cells and their non-inserted neighboring cells per gene
    in one roi.
    :param tmpt: he timepoint - the amount fo hours the cells were inserted for before measuring gene expression
    :param roi: the region of interest in the given timepoint
    :param radius: the radius used for finding the neighboring cells of the inserted cells
    :return:
    """
    all_gene_corrs =[]
    adata = load_sprinkled_adata_hr_tmpt(tmpt, roi)
    non_GFP_indices_cond = ~adata.obs['spc']  # Boolean mask for cells not of the excluded type
    non_GFP_coordinates = adata.obsm['spatial'][non_GFP_indices_cond].values  # Coordinates of valid cells
    # Map the filtered valid indices back to the original DataFrame
    non_GFP_idx_map = np.where(non_GFP_indices_cond)[0]
    adata_spc = adata[adata.obs['spc']]
    tree = BallTree(non_GFP_coordinates)  # the coordinate space
    # query around the sprinkled cells
    indices, distances = tree.query_radius(adata_spc.obsm['spatial'].values, r=radius, return_distance=True)
    # Now map back the indices to the original dataset, the indcies are in order of the non_GFP_coordinates, not the whol adata index space
    filtered_neighbors = []
    for neighbors in indices:
        filtered_neighbors.append(non_GFP_idx_map[neighbors[:4]])

    enough_non_GFP_neighbors_idxs = [i for i in range(len(filtered_neighbors)) if len(filtered_neighbors[i]) >= 4]
    non_GFP_neighbors = [neigh for i, neigh in enumerate(filtered_neighbors) if i in enough_non_GFP_neighbors_idxs]
    for gene in ORGANOID_GENE_NAMES_NOGFP:
        non_GFP_neighbors_exp_pos = np.array(adata[:, gene].X.flatten())[np.array(non_GFP_neighbors)]
        GFP_corr, pval = stats.pearsonr(adata_spc[:, gene].X[enough_non_GFP_neighbors_idxs].flatten(),
                                        non_GFP_neighbors_exp_pos[:, :3].mean(axis=1).flatten())
        all_gene_corrs.append(GFP_corr)

    return np.array(all_gene_corrs)