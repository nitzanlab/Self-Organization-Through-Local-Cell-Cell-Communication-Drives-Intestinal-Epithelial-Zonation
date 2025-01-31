"""
This python file includes functions that allow loading the csv data into anndata for analyses
"""
from utils.imports import *
from utils.constant import *

####Unperturbed Monolayer Data Loader #####

def load_unperturbed_intestinal_organoid_cell_by_gene_mat():
    cell_by_gene_file_path = os.path.join(DATA_DIR, 'cell_by_gene_cluster_annotations.csv')
    #cell_by_gene_file_path = os.path.join(DATA_DIR,'cell_by_gene_mat.csv')

    # Load CSV file into a Pandas DataFrame
    cell_by_gene_data_organoid = pd.read_csv(cell_by_gene_file_path)

    cell_by_gene_data_organoid = cell_by_gene_data_organoid.sort_values(by='object_id', ascending=True)

    #cell_by_gene_data_organoid = cell_by_gene_data_organoid.sort_values(by='Unnamed: 0', ascending=True)
    cell_by_gene_data_organoid.drop('Unnamed: 0', axis=1, inplace=True)
    return cell_by_gene_data_organoid

def load_unperturbed_cell_coords():
    coordinates_file_path = os.path.join(DATA_DIR,'segmentation_cells.csv')
    coordinates_data = pd.read_csv(coordinates_file_path)
    coordinates_data = coordinates_data.sort_values(by='label', ascending=True)
    return coordinates_data

def get_unperturbed_monolayer_adata(num_neigh=5):
    cell_by_gene = load_unperturbed_intestinal_organoid_cell_by_gene_mat()[ORGANOID_GENE_NAMES_NOGFP]
    adata = ad.AnnData(X=cell_by_gene)
    cell_coords = load_unperturbed_cell_coords()
    sc.pp.normalize_total(adata)
    adata.var_names = cell_by_gene.columns
    cell_coords_reindexed = cell_coords.rename(index=dict(zip(cell_coords.index, adata.obs_names)))
    adata.obsm[COORDINATES] = cell_coords_reindexed
    num_neighs_to_use = num_neigh + 1
    nbrs = NearestNeighbors(n_neighbors=num_neighs_to_use, algorithm='auto').fit(adata.obsm[COORDINATES][[X_COORDINATES,Y_COORDINATES]])
    distances, neigh_idxs = nbrs.kneighbors(adata.obsm[COORDINATES][[X_COORDINATES,Y_COORDINATES]])
    adata.obsm['neighbors_idx'] = np.array(neigh_idxs[:, 1:])
    return adata

def load_unperturbed_monolayer_genes_morans_i():
    morans_i_unperturbed_df = pd.read_csv(os.path.join(WT_MONOLAYER_DIR, 'morans_i_wt_monolayer.csv'), index_col=0)
    return morans_i_unperturbed_df
def load_unperturbed_monolayer_transcripts():
    data = pd.read_csv(os.path.join(WT_MONOLAYER_DIR,'transcrips_20240925.csv'))
    return data

def load_unperturbed_monolayer_gene_densities():
    with open(
            os.path.join(WT_MONOLAYER_DIR,'wt_monolayer_gene_densities_new.pkl'),
            'rb') as f:
        result_dict = pickle.load(f)
    return result_dict


def map_monolayer_to_transcript_density_profiles(adata, genes, begin=3, end=14, binned=True):
    genes_in_monolayer = [gene for gene in genes if gene in adata.var_names]
    transcript_df = get_transcript_density_profile_for_monolayer_mapping(begin, end, genes)[genes_in_monolayer]
    if binned:
        binned_df = transcript_df.groupby(np.arange(len(transcript_df)) // 2).mean()
        transcript_df = (binned_df - binned_df.min()) / (binned_df.max() - binned_df.min())
    # sns.heatmap(transcript_df.T)
    # plt.title('transcript')
    # plt.show()
    #adata = adata[adata[:, genes].X.sum(axis=1) > 1e-1]
    #adata_subset = adata[:, genes_in_monolayer].copy()
    adata_X = adata.X.copy()

    # Perform min-max scaling for each gene (column)
    min_vals = adata_X.min(axis=0)
    max_vals = adata_X.max(axis=0)
    # Avoid division by zero for constant columns
    range_vals = np.maximum(max_vals - min_vals, 1e-9)

    # Normalize each gene to the range [0, 1]
    adata.X = (adata_X - min_vals) / range_vals
    cosine_sim_per_position = cosine_similarity(adata[:,genes_in_monolayer].X, transcript_df)
    row_sums = np.sum(cosine_sim_per_position, axis=1)
    #cosine_sim_per_position_normalized = cosine_sim_per_position/row_sums[:,np.newaxis]
    zero_sum_rows = (row_sums == 0)
    cosine_sim_per_position[~zero_sum_rows] = cosine_sim_per_position[~zero_sum_rows] / row_sums[
         ~zero_sum_rows, np.newaxis]
    cosine_sim_per_position[zero_sum_rows] = 1 / cosine_sim_per_position.shape[1]
    adata.obsm['transcript_zone_dist'] = cosine_sim_per_position
    adata.obs['transcript_exp_pos'] = cosine_sim_per_position@np.arange(cosine_sim_per_position.shape[1])
    max_pos = np.argmax(cosine_sim_per_position, axis=1)
    adata.obs['transcript_density_max_pos'] = max_pos
    return adata

def get_transcript_density_profile_for_monolayer_mapping(begin=3, end=13, genes=ORGANOID_GENE_NAMES_NOGFP):
    transcript_density = load_unperturbed_monolayer_gene_densities()
    transcript_df = pd.DataFrame({gene: values['densities'] for gene, values in transcript_density.items()})[
                        genes].iloc[:20]  # like rings
    transcript_df = transcript_df.iloc[begin:end]
    transcript_df_normalized = (transcript_df - transcript_df.min())/(transcript_df.max() - transcript_df.min())

    #TODO additional normalization?
    return transcript_df_normalized

def save_to_pickle_monolayer_masking_components(xedges, yedges, binary_mask_cleaned, extent):
    save_one_monolayer_masking_component_to_pickle(xedges, 'xedges')
    save_one_monolayer_masking_component_to_pickle(yedges, 'yedges')
    save_one_monolayer_masking_component_to_pickle(binary_mask_cleaned, 'binary_mask_cleaned')
    save_one_monolayer_masking_component_to_pickle(extent, 'extent')


def save_one_monolayer_masking_component_to_pickle(component, component_name):
    save_path = os.path.join(WT_MONOLAYER_DIR, f'{component_name}.pkl')
    with open(fr'{save_path}','wb') as f:
        pickle.dump(component, f)

def load_iteration_average_width():
    iteration_widths = load_one_monolayer_masking_component_from_pickle('iteration_widths')
    return np.mean(iteration_widths)

def load_erosion_components():
    xedges = load_one_monolayer_masking_component_from_pickle('xedges')
    yedges = load_one_monolayer_masking_component_from_pickle('yedges')
    binary_mask_cleaned = load_one_monolayer_masking_component_from_pickle('binary_mask_cleaned')
    extent = load_one_monolayer_masking_component_from_pickle('extent')
    ring_masks = load_one_monolayer_masking_component_from_pickle('monolayer_ring_masks')
    return xedges, yedges, binary_mask_cleaned, ring_masks, extent


def load_one_monolayer_masking_component_from_pickle(component_name):
    load_path = os.path.join(WT_MONOLAYER_DIR, f'{component_name}.pkl')
    with open(
            fr'{load_path}','rb') as f:
        component = pickle.load(f)
    return component

def load_transcript_densities_unperturbed_monolayer():
    with open(
            os.path.join(WT_MONOLAYER_DIR, 'wt_monolayer_gene_densities_new.pkl'),
            'rb') as f:
        result_dict = pickle.load(f)
    return result_dict
def save_transcript_densities_unperturbed_monolayer(result_dict):
    with open( os.path.join(WT_MONOLAYER_DIR, 'wt_monolayer_gene_densities_new.pkl'),
            'wb') as f:
        pickle.dump(result_dict, f)

### load sprinkled data by timepoint and roi
def load_cell_by_gene_by_hr_and_roi(hr:str ,roi:str):
    tmpt_dir = os.path.join(MULT_ROIS_DIR, hr)
    roi__dir = os.path.join(tmpt_dir, roi)
    cell_by_gene_path = os.path.join(roi__dir, 'cell_by_gene.csv')
    cell_by_gene = pd.read_csv(cell_by_gene_path, index_col=0)
    return cell_by_gene

def load_cell_coordinates_by_hr_and_roi(hr:str, roi:str):
    tmpt_dir = os.path.join(MULT_ROIS_DIR, hr)
    roi__dir = os.path.join(tmpt_dir, roi)
    cell_coords_path = os.path.join(roi__dir, 'cell_attributes.csv')
    cell_coords = pd.read_csv(cell_coords_path, index_col=0)
    return cell_coords

def load_sprinkled_adata_hr_tmpt(tmpt, roi, num_neigh=5):
    cell_by_gene = load_cell_by_gene_by_hr_and_roi(tmpt, roi)

    cell_coords = load_cell_coordinates_by_hr_and_roi(tmpt, roi)[['center_x', 'center_y']]
    #cell_annotation = load_tacco_annotation_by_hr_roi(tmpt, roi)
    adata = ad.AnnData(X=cell_by_gene)
    adata.var_names = cell_by_gene.columns

    cell_coords_reindexed = cell_coords.rename(index=dict(zip(cell_coords.index, adata.obs_names)))
    #cell_annotation_reindexed = cell_annotation.rename(index=dict(zip(cell_annotation.index, adata.obs_names)))


    adata.obsm['spatial'] = cell_coords_reindexed
    #adata.obsm['cell_type_dist'] = cell_annotation_reindexed
    #adata.obs['cell_type_confusion'] = np.apply_along_axis(entropy, axis=1,
    #                                                                arr=adata.obsm['cell_type_dist'])
    spc_cond = cell_by_gene['GFP'] > SPC_GFP_THRESH
    adata.obs['spc'] = np.array(spc_cond)
    num_neighs_to_use = num_neigh + 1
    sc.pp.filter_cells(adata, min_counts=150)
    # sc.pp.filter_genes(adata, min_counts=10)
    sc.pp.normalize_total(adata, target_sum=1e3)
    # sc.pp.log1p(adata)
    nbrs = NearestNeighbors(n_neighbors=num_neighs_to_use, algorithm='auto').fit(adata.obsm['spatial'])
    distances, neigh_idxs = nbrs.kneighbors(adata.obsm['spatial'])
    adata.obsm['neighbors_idx'] = np.array(neigh_idxs[:, 1:])
    adata = add_GFP_adjacent_obs_to_adata(adata)
    return adata

def add_GFP_adjacent_obs_to_adata(adata):
    spc_neigh_idx = np.unique(np.concatenate(adata[adata.obs['spc']].obsm['neighbors_idx']))
    adata.obs['GFP_adj'] = False
    adata.obs.iloc[spc_neigh_idx, adata.obs.columns.get_loc('GFP_adj')] = True
    return adata

def get_tmpt_all_rois_adata(tmpt, num_neigh=5):
    adata_all_rois = load_sprinkled_adata_hr_tmpt(tmpt, 'roi1', num_neigh)
    for i,roi in enumerate(TMPT_TO_ROIS_DICT[tmpt][1:]):
        adata_one_roi = load_sprinkled_adata_hr_tmpt(tmpt, roi, num_neigh)
        adata_all_rois = adata_all_rois.concatenate(adata_one_roi, batch_key='rois',batch_categories=TMPT_TO_ROIS_DICT[tmpt][:i+2])
    return adata_all_rois

#### in vivo data
def load_TPM_LCM_intestine_atlas():
    LCM_TPM_data = pd.read_csv(os.path.join(SHALEV_DATA_DIR, 'table_A_LCM_TPM_values.tsv'),
                delimiter='\t')
    return LCM_TPM_data

def get_LCM_atlas_gene_subset(atlas, genes):
    LCM_atlas_top_lndrmk, ordered_gene_list = preprocess_LCM_atlas_only_core_reference_genes(atlas, genes)
    return LCM_atlas_top_lndrmk,ordered_gene_list

def load_invivo_reconstruction():
    reconstruction = pd.read_csv(os.path.join(SHALEV_DATA_DIR, 'table_D_zonation_reconstruction.tsv'),
                               delimiter='\t', index_col=0)
    return reconstruction
def preprocess_LCM_atlas_only_core_reference_genes(LCM_atlas ,genes_list):
    preproceed_LCM_atlas = {}
    filted_LCM_atlas = LCM_atlas[LCM_atlas['external_gene_name'].isin(genes_list)]
    filtered_gene_names = LCM_atlas['external_gene_name'][LCM_atlas['external_gene_name'].isin(genes_list)]
    ordered_filtered_genes_list = [gene for gene in genes_list if gene in np.array(filtered_gene_names)]

    for i in range(NUM_POSITIONS_LCM_ATLAS):
        villus_columns = [col for col in filted_LCM_atlas.columns if col.startswith(f'Villus_{i+1}')]
        preproceed_LCM_atlas[f'Villus_{i+1}'] = filted_LCM_atlas[villus_columns].mean(axis=1)
    LCM_atlas_meaned = pd.DataFrame(preproceed_LCM_atlas)
    min_valus = LCM_atlas_meaned.min(axis=1)
    max_valus = LCM_atlas_meaned.max(axis=1)

    LCM_atlas_meaned_normalized = (LCM_atlas_meaned.sub(min_valus,axis=0)).div(max_valus-min_valus,axis=0)
    #set same index as the filtered LCM atlas with the multiple mice , before meaned and normalized
    LCM_atlas_meaned_normalized.set_index(filted_LCM_atlas['external_gene_name'], inplace=True)
    LCM_rearranged = LCM_atlas_meaned_normalized.loc[ordered_filtered_genes_list]
    return LCM_rearranged, ordered_filtered_genes_list

def get_invivo_smooth_exp(gene_names):
    invivo_exp_raw = load_TPM_LCM_intestine_atlas()
    invivo_exp, genes = get_LCM_atlas_gene_subset(invivo_exp_raw, gene_names)
    diffs = np.diff(invivo_exp, axis=1)
    # Identify rows that are monotonically increasing or decreasing
    monotonic_increasing = np.all(diffs >= 0, axis=1)
    monotonic_decreasing = np.all(diffs <= 0, axis=1)
    # Count the number of sign changes (inflection points)
    sign_changes = np.sum(np.diff(np.sign(diffs), axis=1) != 0, axis=1)
    # Identify rows that have exactly one sign change (peaking rows)
    peaking_rows = (sign_changes == 1)
    # Combine conditions: Monotonic or peaking rows
    valid_rows = monotonic_increasing | monotonic_decreasing | peaking_rows
    # Get the indices of rows that satisfy these conditions
    valid_row_indices = np.where(valid_rows)[0]
    invivo_smooth_exp = invivo_exp.iloc[valid_row_indices]
    #sns.heatmap(invivo_smooth_exp, cmap='Reds')
    #plt.show()
    return invivo_smooth_exp


def save_spatial_signal(adata, signal_name, save_name, sprinkled=False):
    signal_df = pd.DataFrame(adata.X, columns=adata.var_names)
    signal_df['x'] = np.array(adata.obsm[COORDINATES][X_COORDINATES])
    signal_df['y'] = np.array(adata.obsm[COORDINATES][Y_COORDINATES])
    signal_df['signal'] = np.array(adata.obs[signal_name])
    if sprinkled:
        signal_df['sprinkled'] = np.array(adata.obs['spc'])
    signal_df.to_csv(os.path.join(WT_MONOLAYER_DIR, f'{save_name}.csv'))