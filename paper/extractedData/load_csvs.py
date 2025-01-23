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
    adata.obsm['spatial'] = cell_coords_reindexed
    num_neighs_to_use = num_neigh + 1
    nbrs = NearestNeighbors(n_neighbors=num_neighs_to_use, algorithm='auto').fit(adata.obsm['spatial'][['center_x','center_y']])
    distances, neigh_idxs = nbrs.kneighbors(adata.obsm['spatial'][['center_x','center_y']])
    adata.obsm['neighbors_idx'] = np.array(neigh_idxs[:, 1:])
    return adata

def load_unperturbed_monolayer_transcripts():
    data = pd.read_csv(os.path.join(WT_MONOLAYER_DIR,'transcrips_20240925.csv'))
    return data

def load_unperturbed_monolayer_gene_densities():
    with open(
            os.path.join(WT_MONOLAYER_DIR,'wt_monolayer_gene_densities_new.pkl'),
            'rb') as f:
        result_dict = pickle.load(f)
    return result_dict

#TODO load erosion plots
def save_to_pickle_monolayer_masking_components(xedges, yedges, binary_mask_cleaned, extent):
    save_one_monolayer_masking_component_to_pickle(xedges, 'xedges')
    save_one_monolayer_masking_component_to_pickle(yedges, 'yedges')
    save_one_monolayer_masking_component_to_pickle(binary_mask_cleaned, 'binary_mask_cleaned')
    save_one_monolayer_masking_component_to_pickle(extent, 'extent')


def save_one_monolayer_masking_component_to_pickle(component, component_name):
    save_path = os.path.join(WT_MONOLAYER_DIR, f'{component_name}.pkl')
    with open(fr'{save_path}','wb') as f:
        pickle.dump(component, f)

def load_erosion_components():
    xedges = load_one_monolayer_masking_component_from_pickle('xedges')
    yedges = load_one_monolayer_masking_component_from_pickle('yedges')
    binary_mask_cleaned = load_one_monolayer_masking_component_from_pickle('binary_mask_cleaned')
    extent = load_one_monolayer_masking_component_from_pickle('extent')
    ring_masks = load_one_monolayer_masking_component_from_pickle('ring_masks')
    return xedges, yedges, binary_mask_cleaned, ring_masks, extent


def load_one_monolayer_masking_component_from_pickle(component_name):
    load_path = os.path.join(WT_MONOLAYER_DIR, f'{component_name}.pkl')
    with open(
            fr'{load_path}','rb') as f:
        component = pickle.load(f)
    return component

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
