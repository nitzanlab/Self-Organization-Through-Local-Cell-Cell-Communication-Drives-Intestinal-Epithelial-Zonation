from utils.imports import *
from utils.constant import *
from paper.extractedData.load_csvs import *
def plot_all_neighborhood_zone_adoption_plots():
    ###panel a: schematic diagram created in  https://BioRender.com
    ###panel b: expected cell zones in unperturbed monolayer
    plot_expected_cell_zones()

    ###panel c: expected zones in inserted cells, subregion of inserted monolayers
    #plot_inserted_cells_expected_zones_spatially()

    ###panel d: correlation in expected zone to neighbors
    #plot_expected_zone_correlations(ORGANOID_GENE_NAMES_NOGFP)
    #plot_expected_zone_correlations(ZONE_MAPPING_GENES)

    ###panel e: expected zone distribution
    #plot_inserted_cells_expected_zone_distribution()

    ###panel f: inserted cells zone confusion
    #plot_inserted_cells_zone_confusion('12hr')
    #plot_inserted_cells_zone_confusion('72hr')

    ###panel g: gene contribution to zone confusion
    #plot_gene_contribution_to_zone_confusion()


def plot_expected_cell_zones():
    smooth_exp_genes = get_invivo_smooth_exp(ORGANOID_GENE_NAMES_NOGFP).index
    plot_expected_position_mapping_unperturbed_cosine_sim(smooth_exp_genes, 'smooth exp genes', to_plot=True)
    pass


def map_monolayer_to_invivo_expression_profiles(adata, genes):
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

def plot_expected_position_mapping_unperturbed_cosine_sim(genes, gene_title='',zoned=False, x_range=None, y_range=None, to_plot=True, to_save=False, save_name=''):
    adata = get_unperturbed_monolayer_adata()
    #adata = map_monolayer_to_invivo_expression_profiles(adata, ZONE_MAPPING_GENES)
    adata = map_monolayer_to_transcript_density_profiles(adata,ZONE_MAPPING_GENES)
    s=6
    if zoned:
        aoi = (adata.obsm['spatial']['center_x'] > x_range[0]) & (
                adata.obsm['spatial']['center_x'] < x_range[1]) & (
                          adata.obsm['spatial']['center_y'] > y_range[0]) & (
                          adata.obsm['spatial']['center_y'] < y_range[1])

        adata = adata[aoi]
        s=10

    vmin = adata.obs['transcript_exp_pos'].min()
    vmax = adata.obs['transcript_exp_pos'].max()
    sc1 = plt.scatter(adata.obsm['spatial']['center_x'], -adata.obsm['spatial']['center_y'], cmap='viridis',
                      c=adata.obs['transcript_exp_pos'], vmin=vmin, vmax=vmax, s=s)


    cbar = plt.colorbar(sc1, label='expected position')
    plt.title(f'unperturbed monolayer {gene_title}\n in vivo expression-based zone mapping ')
    plt.xticks([])
    plt.yticks([])
    if to_plot:
        plt.show()
    if to_save:
        sprinkled_zonation_dir = os.path.join(WT_MONOLAYER_DIR, 'sprinkled_zonation')
        zonation_dir = os.path.join(sprinkled_zonation_dir, f'unperturbed_{save_name}')
        os.makedirs(zonation_dir, exist_ok=True)
        save_path = os.path.join(zonation_dir, f"zones_plot.png")
        plt.savefig(save_path)

def plot_inserted_cells_expected_zones_spatially():
    pass

def plot_expected_zone_correlations(genes):
    #plot_all_expected_position_by_transcript_corr_with_GFP_neighbors(genes):
    wt_env_zone_corr = get_env_transcript_density_zonation_correlation_in_wt_monolayer(genes)
    corrs_GFP = {'72hr': [], '12hr': []}
    corrs_non_GFP = {'72hr': [], '12hr': []}
    mean_GFP_corrs_tmpts = []
    std_GFP_corrs_tmpts = []
    mean_nonGFP_corrs_tmpts = []
    std_nonGFP_corrs_tmpts = []
    for tmpt in TMPT_TO_ROIS_DICT.keys():
        for roi in TMPT_TO_ROIS_DICT[tmpt]:
            GFP_corr, nonGFP_corr = get_neighbor_corrs(tmpt, roi, genes)
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
    plt.ylabel('Expected Zone Correlation', fontsize=20)
    plt.title('Expected Zone Correlation - closest neighbors', fontsize=22)
    plt.tight_layout(pad=4)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

def get_neighbor_corrs(tmpt, roi, genes, radius=500):
    adata = load_sprinkled_adata_hr_tmpt(tmpt, roi)
    adata = map_monolayer_to_transcript_density_profiles(adata, genes)
    #find k nearest neighbors in radius for all cells
    tree = BallTree(adata.obsm['spatial'].values)

    indices, distances = tree.query_radius(adata.obsm['spatial'].values, r=radius, return_distance=True)
    enough_indices_idxs = [i for i in range(adata.shape[0]) if len(indices[i]) >= 3]
    neighbors_of_enough_cells = [indices[i][:3] for i in range(adata.shape[0]) if len(indices[i])>=3]
    adata_neighs = adata[enough_indices_idxs,:]
    adata_neighs.obs['neighbor_mean_transcript_zone'] = np.array(adata.obs['transcript_exp_pos'])[neighbors_of_enough_cells].mean(axis=1)
    adata_neigh_spc = adata_neighs[adata_neighs.obs['spc']]
    adata_neigh_non_spc = adata_neighs[~adata_neighs.obs['spc']]
    GFP_zone_corr, GFP_pval = pearsonr(adata_neigh_spc.obs['transcript_exp_pos'],  adata_neigh_spc.obs['neighbor_mean_transcript_zone'])
    non_GFP_zone_corr, non_GFP_pval = pearsonr(adata_neigh_non_spc.obs['transcript_exp_pos'], adata_neigh_non_spc.obs['neighbor_mean_transcript_zone'])
    print( GFP_zone_corr, GFP_pval ,f'GFP {tmpt} {roi}')
    print( non_GFP_zone_corr, non_GFP_pval , f'non GFP {tmpt} {roi}')
    return GFP_zone_corr, non_GFP_zone_corr

def get_env_transcript_density_zonation_correlation_in_wt_monolayer(genes):
    wt_adata = get_unperturbed_monolayer_adata()
    wt_adata = map_monolayer_to_transcript_density_profiles(wt_adata, genes, 3, 14)
    wt_neighbors = np.array(wt_adata.obs['transcript_exp_pos'])[np.array(wt_adata.obsm['neighbors_idx'])]
    correlation_in_exp_pos, pval = stats.pearsonr(wt_adata.obs['transcript_exp_pos'],
                                                  wt_neighbors.mean(axis=1))
    print(f'corr: {correlation_in_exp_pos} pval:{pval:.3f}')
    return correlation_in_exp_pos

def plot_inserted_cells_expected_zone_distribution():
    expected_zone_distribution('72hr',ZONE_MAPPING_GENES,3, 14)
    expected_zone_distribution('12hr', ZONE_MAPPING_GENES, 3, 14)

def expected_zone_distribution(tmpt, genes, begin, end):
    adata = get_tmpt_all_rois_adata(tmpt)
    adata = map_monolayer_to_transcript_density_profiles(adata, genes, begin, end, binned=True)
    adata = get_transcript_density_position_std(adata)
    bins = np.linspace(0, 5, 20)
    adata_non_spc = adata[~adata.obs['spc']]
    adata_spc = adata[adata.obs['spc']]
    plt.hist(adata_non_spc.obs['transcript_exp_pos'], density=True, bins=bins, label='non GFP', alpha=0.7)
    plt.hist(adata_spc.obs['transcript_exp_pos'], density=True, bins=bins, label='GFP', alpha=0.7)
    plt.title(f'expected position distribution {tmpt}')
    stat, p_value = ks_2samp(adata_non_spc.obs['transcript_exp_pos'], adata_spc.obs['transcript_exp_pos'])
    print(f"KS Statistic: {stat}, p-value: {p_value}")
    plt.legend()
    plt.show()

def get_transcript_density_entropy(adata):
    zone_entropy = np.apply_along_axis(lambda x: entropy(x, base=2), axis=1, arr=adata.obsm['transcript_zone_dist'])
    adata.obs['zone_entropy'] = zone_entropy
    return adata
def get_transcript_density_position_std(adata):
    exp_pos = np.array(adata.obs['transcript_exp_pos'])
    pos_dist = adata.obsm['transcript_zone_dist']

    squared_diff = (np.tile(np.arange(pos_dist.shape[1]), (len(exp_pos), 1)) - exp_pos[:,
                                                                               np.newaxis]) ** 2  # Shape: (n_samples, n_positions)
    weighted_squared_diff = pos_dist * squared_diff
    variance = np.sum(weighted_squared_diff, axis=1)
    pos_std = np.sqrt(variance)
    adata.obs['transcript_pos_std'] = pos_std
    return adata

def plot_inserted_cells_zone_confusion(tmpt):
    adata = get_tmpt_all_rois_adata(tmpt)
    adata = adata[adata[:, ZONE_MAPPING_GENES].X.sum(axis=1) > 0.1]
    adata = map_monolayer_to_transcript_density_profiles(adata, ZONE_MAPPING_GENES, 3, 14)
    adata = get_transcript_density_entropy(adata)

    adata_gfp = adata[adata.obs['spc']]
    adata_non_gfp = adata[~adata.obs['spc']]
    bins = np.linspace(1.6, 2.6, 20)
    plt.hist(adata_non_gfp.obs['zone_entropy'], alpha=0.7, label='non GFP', density=True, color='purple', bins=bins)
    plt.hist(adata_gfp.obs['zone_entropy'], alpha=0.7, label='GFP', density=True, color='grey', bins=bins)
    plt.xlabel('zone entropy')
    plt.ylabel('density')
    plt.ylim(0, 6)
    plt.title(f'{tmpt} Zone Entropy')
    plt.legend()
    plt.show()
    stat, p_value = ks_2samp(adata_gfp.obs['transcript_exp_pos'], adata_non_gfp.obs['transcript_exp_pos'])
    print(f"KS Statistic: {stat}, p-value: {p_value} {tmpt}")




def gene_confusion_contribution(tmpt, genes,begin, end, is_sprinkled= True):
    #instead of mapping cells to positions using the transcript density profiles
    #get transcript density:
    transcript_df = get_transcript_density_profile_for_monolayer_mapping(begin, end, genes)
    binned_df = transcript_df.groupby(np.arange(len(transcript_df)) // 2).mean()
    transcript_df = (binned_df-binned_df.min())/(binned_df.max()-binned_df.min())
    # sns.heatmap(transcript_df.T)
    # plt.title('trasncript df')
    # plt.show()
    #get gene confusion contribution per cell
    adata = load_sprinkled_adata_hr_tmpt(tmpt, 'roi2')

    adata = map_monolayer_to_transcript_density_profiles(adata, genes, begin, end, binned=True)
    adata = get_transcript_density_position_std(adata)
    if is_sprinkled:
        adata = adata[adata.obs['spc']]
        title='GFP'
    else:
        adata = adata[~adata.obs['spc']]
        title = 'non-GFP'
    full_title = f'{title} {tmpt}'
    bins = np.linspace(0,140,20)
    # for gene in ZONE_MAPPING_GENES:
    #     plt.hist(adata[:,gene].X, density=True, bins=bins)
    #     plt.title(f'{gene} expression in {full_title}')
    #     plt.ylim(0,0.15)
    #     plt.show()
    plot_confusion_contribution_cell_groups(adata, genes, transcript_df, 2, title=full_title)


def plot_confusion_contribution_cell_groups(adata, genes, transcript_df, num_cells, title=''):
    #plot examples of top and bottom cells
    #can also do this per zone
    obs_df = adata.obs['transcript_pos_std']
    top_indices = adata.obs.nlargest(num_cells, 'transcript_pos_std').index
    bottom_indices = adata.obs.nsmallest(num_cells, 'transcript_pos_std').index
    low_pos_std_exp = []
    for j, ind_top in enumerate(bottom_indices):
        cell = adata[ind_top, genes]
        dot_prod = np.multiply(cell.X, transcript_df)
        norms = np.maximum(np.linalg.norm(cell.X) * np.linalg.norm(transcript_df, axis=1),
                           np.ones(transcript_df.shape[0]))
        cell_confusion_contribution = dot_prod / norms[:, np.newaxis]
        sns.heatmap(expit(cell_confusion_contribution.T))
        plt.title(f"{j} lowest pos std cell gene confusion contribution\n {title}\n {float(cell.obs['transcript_pos_std']):.2f}")
        plt.show()
    #     low_pos_std_exp.append(cell.X[0])
    # low_pos_std_df = pd.DataFrame(np.array(low_pos_std_exp), columns=genes).T
    # low_pos_std_df = (low_pos_std_df - low_pos_std_df.min()) / (low_pos_std_df.max()-low_pos_std_df.min())
    # sns.heatmap(low_pos_std_df)
    # plt.title(f'cell mapping gene expression, low pos std\n {title}')
    # plt.show()

    high_pos_std_exp = []
    for i,ind in enumerate(top_indices):
        cell = adata[ind, genes]
        dot_prod = np.multiply(cell.X, transcript_df)
        norms = np.maximum(np.linalg.norm(cell.X) * np.linalg.norm(transcript_df, axis=1),
                           np.ones(transcript_df.shape[0]))
        cell_confusion_contribution = dot_prod / norms[:, np.newaxis]
        sns.heatmap(expit(cell_confusion_contribution.T))
        plt.title(f"{i} highest pos std cell gene confusion contribution\n {title}\n {float(cell.obs['transcript_pos_std']):.2f}")
        plt.show()
        #high_pos_std_exp.append(cell.X[0])
    # high_pos_std_exp_df = pd.DataFrame(np.array(high_pos_std_exp), columns=genes).T
    # high_pos_std_exp_df = (high_pos_std_exp_df - high_pos_std_exp_df.min())/ high_pos_std_exp_df - high_pos_std_exp_df.min()
    # sns.heatmap(high_pos_std_exp_df)
    # plt.title(f'cell mapping gene expression, high pos std\n {title}')
    # plt.show()

def plot_pos_std_signal_spatially_all_monolayers(genes):
    # adata_up = get_adata_with_zonation_signals_up(genes)
    # plot_adata_signal_spatially(adata_up, 'pos_std', 'unperturbed')
    for tmpt in ['12hr']:#TMPT_TO_ROIS_DICT:
        for roi in ['roi1']:# TMPT_TO_ROIS_DICT[tmpt]:
            adata_GFP = get_adata_with_zonation_GFP_monolayer(tmpt, roi, genes)
            plot_GFP_adata_signal_spatially(adata_GFP, 'exp_pos', f'{tmpt} {roi} GFP',zoned=True, x_range=GFP_12HR_ROI1_X, y_range=GFP_12HR_ROI1_Y)



def get_adata_with_zonation_GFP_monolayer(tmpt, roi, genes, num_neighs=5):
    adata = load_sprinkled_adata_hr_tmpt(tmpt, roi, num_neighs)
    #adata_exp_pos = map_monolayer_to_invivo_expression_profiles(adata, genes)
    adata_pos_std = get_position_std(adata, genes)
    return adata_pos_std

def get_position_std(adata, genes):
    adata = map_monolayer_to_invivo_expression_profiles(adata, genes)
    exp_pos = np.array(adata.obs['invivo_exp_pos'])
    pos_dist = adata.obsm['invivo_position_dist']

    squared_diff = (np.tile(np.arange(pos_dist.shape[1]), (len(exp_pos),1)) - exp_pos[:, np.newaxis]) ** 2  # Shape: (n_samples, n_positions)
    weighted_squared_diff = pos_dist * squared_diff
    variance = np.sum(weighted_squared_diff, axis=1)
    pos_std = np.sqrt(variance)
    adata.obs['pos_std'] = pos_std
    return adata


def plot_GFP_adata_signal_spatially(adata, signal_name, adata_type, zoned=False, x_range=None, y_range=None, to_plot=True):
    adata_spc = adata[adata.obs['spc']]
    adata_ent = adata[adata.obsm['cell_type_dist']['enterocyte'] > 0.3]
    plt.figure(figsize=(8, 6))
    s = 30
    if zoned:
        plt.figure(figsize=(6, 4))
        aoi_spc = (adata_spc.obsm['spatial']['center_x'] > x_range[0]) & (
                adata_spc.obsm['spatial']['center_x'] < x_range[1]) & (
                          adata_spc.obsm['spatial']['center_y'] > y_range[0]) & (
                          adata_spc.obsm['spatial']['center_y'] < y_range[1])
        aoi_ent = (adata_ent.obsm['spatial']['center_x'] > x_range[0]) & (
                adata_ent.obsm['spatial']['center_x'] < x_range[1]) & (
                          adata_ent.obsm['spatial']['center_y'] > y_range[0]) & (
                          adata_ent.obsm['spatial']['center_y'] < y_range[1])
        adata_spc = adata_spc[aoi_spc]
        adata_ent = adata_ent[aoi_ent]


    vmin = min(adata_spc.obs[signal_name].min(), adata_ent.obs[signal_name].min())
    vmax = max(adata_spc.obs[signal_name].max(), adata_ent.obs[signal_name].max())

    sc1 = plt.scatter(adata_ent.obsm['spatial']['center_x'], adata_ent.obsm['spatial']['center_y'], cmap='plasma',
                      c=adata_ent.obs[signal_name], vmin=vmin, vmax=vmax, s=s)
    plt.scatter(adata_spc.obsm['spatial']['center_x'], adata_spc.obsm['spatial']['center_y'], label='GFP',
                edgecolors='lime', cmap='plasma', c=adata_spc.obs[signal_name], vmin=vmin, vmax=vmax, s=s)

    cbar = plt.colorbar(sc1, label='Color Attribute')
    plt.title(f'{adata_type} monolayer {signal_name}\n ')
    plt.xticks([])
    plt.yticks([])
    plt.show()