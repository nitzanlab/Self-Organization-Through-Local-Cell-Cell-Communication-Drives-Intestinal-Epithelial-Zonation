from utils.imports import *

def set_style():
 plt.rcParams.update(plt.rcParamsDefault)
 plt.rcParams.update({'figure.titlesize': 8, 'figure.titleweight': 'bold', 'axes.titlesize': 8,
                      'axes.titleweight': "bold", 'axes.labelsize': 8, 'axes.labelweight': 'bold',
                      "ytick.labelsize": 8, "xtick.labelsize": 8, 'legend.fontsize': 8,
                      'font.family': 'DejaVu Sans'})
 plt.rcParams.update({'figure.figsize': (4, 3)})
 plt.rcParams.update({'savefig.dpi': 300})


#cell type constants
CELL_TYPES_TO_MARKER_GENES = {'enterocyte':['Alpi', 'Aldob','Sis','Apoa1'], 'goblet_cells':['Muc2'],'EEC':['Chga'],'Tuft_cells':['Dclk1'],'paneth_cells':['Lyz1'], 'stem_cells':['Lgr5','Olfm4'], 'regenerative':['Msln','Ahnak']}
CELL_TYPES_TO_MARKER_GENES_REG_RESPONSE = {'enterocyte':['Aldob','Ada','Apoa4','Apoa1','Alpi','Sis','Apob'], 'goblet_cells':['Muc2'],'EEC':['Chga'],'Tuft_cells':['Dclk1'],'paneth_cells':['Lyz1'], 'stem_cells':['Lgr5','Olfm4'], 'regenerative':['Mki67','Clu' ,'Yap1','Ly6a','Msln','Ahnak']}

UNPERTURBED_ZOOMED_IN_X = [10000, 17000]
UNPERTURBED_ZOOMED_IN_Y = [14000, 21000]

UNPERTURBED_ZOOMED_IN_X_SEC_CELL = [15000, 16500]
UNPERTURBED_ZOOMED_IN_Y_SEC_CELL = [14000, 21000]

# UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2 = [12000, 16000]
# UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2 = [17500, 21000]

UNPERTURBED_ZOOMED_IN_PROGENITOR_X= [13000, 15000]
UNPERTURBED_ZOOMED_IN_PROGENITOR_Y= [17800, 18800]

UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2 = [13500, 15500]
UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2 = [17750, 19500]


#EROSION constants
XY_SPACING = 10
EROSION_STEP = 5
NUM_ITERATIONS = 30

AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH = os.path.join(os.getcwd(), 'paper','graphs','autonomous_zonation_figure_plots')





CELL_TYPE_ANNOT_ORDER = ['regenerative', 'enterocyte','EEC','goblet_cells','stem_cells','paneth_cells','Tuft_cells']
ALL_CELL_TYPES_PLOT_ORDER =['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells','enterocyte','regenerative']
DOMINANT_CELL_TYPES = ['regenerative', 'enterocyte']
CELl_TYPE_WO_ENT_REG_PLOT_ORDER = ['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells']
CELl_TYPEs_WO_REGENERATIVE = ['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells','enterocyte']
##gene constants

DATA_DIR= 'C:/Users/micha/thesis/code/data/intestinal_organoid/non_sprinkled_july23_pasadena/'
ORGANOID_GENE_NAMES = np.load(os.path.join(DATA_DIR, 'organoid_gene_names.npy'))
ORGANOID_GENE_NAMES_NOGFP = ORGANOID_GENE_NAMES[ORGANOID_GENE_NAMES !='GFP']
ORGANOID_NO_CLU_GENE_NAMES = ORGANOID_GENE_NAMES_NOGFP[ORGANOID_GENE_NAMES_NOGFP!='Clu']
SPATIAL_CORRELATION_GENE_PANELS = ['Aldob', 'Chgb']
ENTEROCYTE_GENES =  ['Alpi','Sis','Apoa1','Aldob', 'Apoa4', 'Ada']

ALL_CELL_TYPE_GENES=  ['Ada', 'Apoa4', 'Apoa1', 'Alpi', 'Sis','Aldob','Clu','Msln','Ahnak','Muc2','Chga','Dclk1','Lyz1','Lgr5','Olfm4']

ENT_REG_GENES = ['Aldob', 'Alpi', 'Pigr', 'Apoa4', 'Ly6d', 'Ccna2', 'Mki67', 'Chga', 'Muc2', 'Lyz1', 'Dclk1',
                     'Msln', 'Clu', 'Olfm4', 'Lgr5']
ENT_GEG_GENES_2 = ['Ada', 'Sis', 'Apob', 'Ahnak', 'Apoa1']
ENT_REG_ALL_GENES = np.concatenate((ENT_REG_GENES, ENT_GEG_GENES_2))
KEY_ENT_AND_REG_GENES = ['Ada', 'Apoa4', 'Apoa1', 'Alpi', 'Sis','Clu','Msln','Ahnak']
ZONATION_COMPARISON_GENES =  ['Ada', 'Apoa4', 'Lgals3', 'Slc28a2', 'Aldob', 'Ccl25', 'Fabp1', 'Lypd8', 'Reg3b', 'Reg3g', 'Sis',
                    'Spink4']
INVIVO_ZONATION_GENES = ['Neat1','Malat1','Reg3g','Reg3b','Reg3a','Nlrp6','Lypd8','Il18','Ccl25', 'Apobec1','Apob',
                         'Apoa4','Apoa1','NPC1L1','Slc15a1','Slc5a1','Slc2a5','Slc2a2','Slc7a9','Slc7a8','Slc7a7', 'Alpi','Aldob','Sis','Ada','Nt5e']

genes_in_order_density_measure = ['Anxa5','Ada', 'Apoa4', 'Apoa1', 'Alpi', 'Aldob', 'Sis']
genes_in_order_density_measure_reversed = ['Sis', 'Aldob', 'Alpi', 'Apoa1', 'Apoa4', 'Anxa5']

REGENERATIVE_GENES = ['Ahnak','Msln','Yap1','Ly6a','Clu']
REG_SUB_GROUP = ['Ahnak','Msln','Clu']

UNPERTURBED_CELL_TYPE_CLUSTERS = {0:'regenerative1',1:'regenerative2',2:'regenerative3',3:'enterocyte',4:'secretory', 5:'progenitor',6:'regenerative',7:'unknown'}

FUNCTIONAL_GENES = ['Apobec1', 'Apob', 'Apoa4', 'Apoa1', 'Nlpc1l1', 'Slc15a1', 'Slc5a1', 'Slc2a5', 'Slc2a2', 'Slc7a9',
                    'Slc7a8', 'Slc7a7']
CHOLESTEROL_GENES = ['Apobec1', 'Apob', 'Apoa4', 'Apoa1', 'Nlpc1l1']
PEPTIDE_GENE = ['Slc15a1']
CARBOHYDRATE_GENES = ['Slc5a1', 'Slc2a5', 'Slc2a2']
AMINO_ACID_GENES = ['Slc7a9', 'Slc7a8', 'Slc7a7']
reg_genes = ['Reg3g', 'Reg3b', 'Reg3a', 'Nlrp6', 'Lypd8', 'Ill18', 'Ccl25']

TRANSCRIPT_DENSITY_GENES = ['Anxa5','Apoa4','Apoa1','Alpi','Aldob','Sis']

ZONE_MAPPING_GENES_WO_ALDOB =  ['Ada', 'Apoa4', 'Apoa1', 'Alpi', 'Sis']
ZONE_MAPPING_GENES_WO_ALPI =  ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Sis']
ZONE_MAPPING_GENES_WO_SIS = ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Alpi']
ZONE_MAPPING_GENES_WO_ADA = ['Apoa4', 'Apoa1', 'Aldob', 'Alpi', 'Sis']
ZONE_MAPPING_GENES = ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Alpi', 'Sis']
ZONE_MAPPING_GENES_CHECK = ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Alpi']

ZONE_MAPPING_WITH_REG = ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Alpi', 'Sis','Msln','Clu','Ahnak']
##invivo

BHB_HIGH_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\spatial_substance_perturb\20241015_150428_062__WellB07_ChannelCY5,A594,CY3,YFP,DAPI,Brightfield_Seq0000 - Stitched_crop.nd2'
PERTURB_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\spatial_substance_perturb'
PERTURB_CHANNELS = ['cy3 ', 'cy5 ', 'yfp ', 'a594_']
PERTURBED_DATASETS = ['bmp 4 high.csv','bmp 2 high 2.csv','iwp-2 - wnt inhibitor high.csv','ldn - bmp inhibitor - low.csv']
CHANNELS_TO_GENES = {'yfp ': ['Olfm4','Lgr5'], 'cy3 ' : 'Sis', 'a594_' : 'Msln', 'cy5 ' : 'Apoa4'}
PERTURB_ENT_GENES = ['Sis','Apoa4']
PERTURB_SIS_THRESH = 6000
PERTURB_APOA4_THRESH = 7100
PERTURB_MAX_EXP = {'bmp 2 high 2':0.75, 'bmp 2 high':0.75,'bmp 2 low 2':0.75,'bmp 4 high 2':0.6,'bmp 4 high':0.5,'bmp 4 low 2':0.8,
                   'bmp2 low':0.8,'iwp-2 - wnt inhibitor high':0.65,'iwp-2 high':0.6, 'iwp-2 low 2':0.65, 'ldn - bmp inhibitor - low':0.65,
                   'ldn low 2':0.65, 'bmp 4 low':0.8, 'iwp-2 high 2':0.6,'iwp-2 wnt inhib low':0.65
                   }


SHALEV_DATA_DIR = 'C:/Users/micha/thesis/code/data/intestinal_organoid/shalev_data/'
NUM_POSITIONS_LCM_ATLAS = 5
WT_MONOLAYER_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena'
GENE_TO_ZONE_DIR =  r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\reconstruction_new\Ada_0.5'
VILLUS_LOCATION_COLUMNS = [ 'V1_mean', 'V2_mean', 'V3_mean', 'V4_mean', 'V5_mean', 'V6_mean']
INITIAL_RES_DIR =  r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\minimal_reference_atlas_new'
CHECK_RES_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\enterocyte_recon_with_atlas'
NOVOSPARC_MAPPING_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\autonomous_zonation_dir\ada_sis'
NOVOSPARC_MAPPING_W_ZONATED_GENES_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\autonomous_zonation_dir\using_only_zonated_genes_for_reconstruction\ada_sis'
NOVOSPARC_MAPPING_DIR_GENERAL = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\autonomous_zonation_dir'
INITIAL_MAPPING = r'C:\Users\micha\thesis\code\data\intestinal_organoid\non_sprinkled_july23_pasadena\minimal_reference_atlas_new\ada\alpha_linear_0_5'
ZONE_CONFUSION_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\sprinkled_cells\zone_confusion'

#novosparc mapping - tissue zonation
TISSUE_SDGE = 'sdge_8662_cells_5_locations.txt'#'sdge_7749_cells_5_locations.txt'
ENT_TISSUE_SDGE = 'sdge_4953_cells_5_locations.txt'
TISSUE_GW = 'cell_to_positions_tp_Ada_Sis.txt'
TISSUE_GW_GENERAL = 'cell_to_positions_tp.txt'
REFERENCE_GENES = ['Ada', 'Apoa4', 'Aldob', 'Sis']
ALL_BOTTOM_LNDMRK_GENES_SHALEV = ['Fabp1', 'Plac8', 'Lgals4', 'Rpl41', 'Cox4i1', 'Cox6c', 'Atp5e', 'Atp5j2', 'Uqcrq',
                                  'Cox7a2',
                                  'Cox6b1', 'Rps14', 'Reg3b', 'Uqcrh', 'Cox5b', 'Uqcr11','Atpif1', 'Txn1', 'Atp5g3',
                                  'Atp5l',
                                  'Uqcr10','Rplp1', 'Cox7b', 'Rps29', 'Rps27', 'Rps18', 'Atp5j', 'Cox5a', 'Ckmt1',
                                  'Rps8',
                                  'Ndufa4', '2010107E04Rik', 'Atp5o', 'Reg3g', 'Lypd8', 'Rpl35a', 'Atp5a1', 'Tm4sf5',
                                  'Chchd10', 'Rps2', 'Cox7a1', 'Atp5h', 'Tma7', 'Rpl18', '2210407C18Rik', 'Sis',
                                  'Rpl38',
                                  'Rps12', 'Rpl39', 'Rps28', 'Rps27l', 'Ndufc1', 'Minos1', 'Ndufa1', 'Ccl25',
                                  'Ndufb8',
                                  'Ndufa5', 'Gsta1', 'Gpx2', 'Uba52', 'Ndufb6', 'Spink4']

ALL_TOP_LNDMRK_GENES_SHALEV = ['Apoa4', 'Apoc3', 'Krt20', 'Lgals3', 'Ada', 'Tmsb4x', 'S100a10', 'Cldn7',
                               '2010109I03Rik',
                               'Clca4a', 'Serpinb1a', 'Pmp22', 'Ifrd1', 'Slc28a2', 'S100a6', 'Psma7', 'Slc25a22',
                               'Chchd1',
                               'Tbk1', 'Fgd4','Ap1g1', 'Mrpl48', 'Gm10680', 'Lrrc41', 'Gcn1l1', 'Slc17a5', 'Sprr2a2',
                               'Cfap20', 'Myo7a', 'Pam', 'Zfp280d', 'Ythdc2', 'Cep57', 'Acad9', 'Chek2']





#wt monolayer mapping to invivo
NUM_NEIGHBORS_S = 5
NUM_NEIGHBORS_T = 2

LOW_ENT_RANGE = [1.6,1.8]
MEDIUM_ENT_RANGE = [1.8,2.0]
HIGH_ENT_RANGE = [2.1,2.3]
VERY_HIGH_ENT_RANGE = [2.4,2.5]
UNIFORM_DIST_RANGE = [2.5,2.6]


##Tacco
TACCO_RESULTS_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\sprinkled_cells\tacco_res'

#GFP cells
SPC_GFP_THRESH = 2

CLU_THRESH = 1



#10x gene set by clusters
GENE_CLUSTER_DIR = r'C:\Users\micha\thesis\code\data\intestinal_organoid\10x_cluster_genes'
GENE_CLUSTER_RES = 'merfish_gene_list_with_all_cluster_umi_counts.csv'


##ROIS constants
TMPT_TO_ROIS_DICT = {'12hr': ['roi1', 'roi2', 'roi3'], '72hr':['roi1', 'roi2', 'roi3', 'roi4']}
MULT_ROIS_DIR = 'C:/Users/micha/thesis/code/data/intestinal_organoid/sprinkled_cells_multiple_rois/'

#GFP cells
SPRINKLED_CELL_DIR = 'C:/Users/micha/thesis/code/data/intestinal_organoid/sprinkled_cells/'

GFP_DIR = 'C:/Users/micha/thesis/code/data/intestinal_organoid/sprinkled_cells/GFP_exp'