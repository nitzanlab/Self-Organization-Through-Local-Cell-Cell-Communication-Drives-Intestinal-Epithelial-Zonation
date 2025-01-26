
from utils.imports import *
from utils.constant import *
from paper.extractedData.load_csvs import *
def plot_zonation_plasticity_plots():
    ###panel a: schematic diagram created in BioRender.com
    ###panel b and c: raw image example #TODO Yael
    ###panel d and e : ### correlation comparison #TODO Yael

    ###panel f: gene expression neighborhood correlation in inserted cells
    #plot_gene_expression_neighborhood_correlation_in_inserted_cells()
    plot_gene_correlation_histograms_GFP_to_wt()

    ###panel g: inserted cell neighborhood correlation across timepoints  #TODO Yael

    pass

def plot_gene_correlation_histograms_GFP_to_wt():
    morans_i_wt = pd.read_csv(os.path.join(WT_MONOLAYER_DIR, 'morans_i_wt_monolayer.csv'), index_col=0)
    num_top_moran_genes = len(morans_i_wt)
    top_20_moran_i_genes = np.array(morans_i_wt.sort_values(by='morans_i', ascending=False).head(num_top_moran_genes).index)
    #get top 20 morani genes index in ORGANOID_GENE_NAMES_NO_GFP
    moran_i_genes_idx = [np.where(ORGANOID_GENE_NAMES_NOGFP == item)[0][0] for item in top_20_moran_i_genes]
    #indices = [ORGANOID_GENE_NAMES_NOGFP.index(item) for item in top_20_moran_i_genes]
    wt_non_GFP_gene_corr = pd.read_csv(os.path.join(WT_MONOLAYER_DIR, 'wt_gene_correlation.csv'), index_col=0)
    wt_non_GFP_gene_corr = wt_non_GFP_gene_corr.reindex(ORGANOID_GENE_NAMES_NOGFP)
    wt_non_GFP_gene_corr = wt_non_GFP_gene_corr['corr'].loc[top_20_moran_i_genes]
    GFP_72hr_gene_corr = get_all_GFP_to_env_one_tmpt_gene_exp_corr('72hr').iloc[moran_i_genes_idx]
    GFP_12hr_gene_corr = get_all_GFP_to_env_one_tmpt_gene_exp_corr('12hr').iloc[moran_i_genes_idx]

    bins = np.linspace(-0.25,0.8,20)
    # plt.hist(GFP_72hr_gene_corr, bins=bins, color= 'orange')
    # plt.axvline(np.mean(GFP_72hr_gene_corr), label='mean inserted 72hr ', color='orange', linestyle='--', linewidth=2)
    # plt.axvline(np.mean(GFP_12hr_gene_corr), label='mean inserted 12hr', color='green',linestyle='--', linewidth=2)
    # plt.axvline(np.mean(wt_non_GFP_gene_corr), label='mean wt ', color='purple',linestyle='--', linewidth=2)
    # plt.legend()
    # plt.title('72hr GFP gene corr')
    # plt.show()
    #
    # plt.hist(wt_non_GFP_gene_corr, bins=bins, color='purple')
    # plt.axvline(np.mean(GFP_12hr_gene_corr), label='mean inserted 12hr', color='green',linestyle='--', linewidth=2)
    # plt.axvline(np.mean(GFP_72hr_gene_corr), label='mean inserted 72hr ', color='orange',linestyle='--', linewidth=2)
    # plt.axvline(np.mean(wt_non_GFP_gene_corr), label='mean wt ', color='purple',linestyle='--', linewidth=2)
    # plt.legend()
    # plt.title('wt gene corr')
    # plt.show()
    #
    # plt.hist(GFP_12hr_gene_corr, bins=bins, color= 'green')
    # plt.axvline(np.mean(GFP_12hr_gene_corr), label='mean inserted 12hr', color='green', linestyle='--', linewidth=2)
    # plt.axvline(np.mean(GFP_72hr_gene_corr), label='mean inserted 72hr', color='orange',linestyle='--', linewidth=2)
    # plt.axvline(np.mean(wt_non_GFP_gene_corr), label='mean wt ',color='purple',linestyle='--', linewidth=2)
    # plt.legend()
    # plt.title('12hr GFP gene corr')
    # plt.show()

    plt.hist(wt_non_GFP_gene_corr, bins=bins, color='purple', alpha=0.8)
    plt.hist(GFP_72hr_gene_corr, bins=bins, color='orange', alpha=0.8)
    plt.hist(GFP_12hr_gene_corr, bins=bins, color='green',alpha=0.8)


    plt.axvline(np.mean(GFP_72hr_gene_corr), label='mean inserted 72hr ', color='orange', linestyle='--', linewidth=2)
    plt.axvline(np.mean(GFP_12hr_gene_corr), label='mean inserted 12hr', color='green', linestyle='--', linewidth=2)
    plt.axvline(np.mean(wt_non_GFP_gene_corr), label='mean unperturbed ', color='purple', linestyle='--', linewidth=2)
    #plt.legend(fontsize=20)
    #plt.xlabel('gene expression correlation')
    #plt.ylabel('frequency')
    #plt.title(f'gene exp env corr\n - top {str(num_top_moran_genes)} morans i genes')
    plt.tight_layout()
    plt.show()

def get_all_GFP_to_env_one_tmpt_gene_exp_corr(tmpt):
    #gene_corrs = {f'{roi}':[] for roi in TMPT_TO_ROIS_DICT[tmpt]}
    gene_corrs = {}
    for roi in TMPT_TO_ROIS_DICT[tmpt]:
        gene_corrs[roi] = get_GFP_to_env_gene_exp_corr_one_roi(tmpt, roi)
    mean_corrs_per_tmpt = pd.DataFrame(gene_corrs).mean(axis=1)
    return mean_corrs_per_tmpt


def get_GFP_to_env_gene_exp_corr_one_roi(tmpt, roi, radius=500):
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
        filtered_neighbors.append(non_GFP_idx_map[neighbors[:6]])

    enough_non_GFP_neighbors_idxs = [i for i in range(len(filtered_neighbors)) if len(filtered_neighbors[i]) >= 6]
    non_GFP_neighbors = [neigh for i, neigh in enumerate(filtered_neighbors) if i in enough_non_GFP_neighbors_idxs]
    for gene in ORGANOID_GENE_NAMES_NOGFP:
        non_GFP_neighbors_exp_pos = np.array(adata[:, gene].X.flatten())[np.array(non_GFP_neighbors)]
        GFP_corr, pval = stats.pearsonr(adata_spc[:, gene].X[enough_non_GFP_neighbors_idxs].flatten(),
                                        non_GFP_neighbors_exp_pos[:, :5].mean(axis=1).flatten())
        all_gene_corrs.append(GFP_corr)

    return np.array(all_gene_corrs)