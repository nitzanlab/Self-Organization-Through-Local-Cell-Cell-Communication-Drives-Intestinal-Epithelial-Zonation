from utils.imports import *

"""
in order to load the data properly and be able to conduct the analyses and reproduce the figure results:
1. download the data from this paper from :
2. download from 'Moor, A. E., Harnik, Y., Ben-Moshe, S., Massasa, E. E., Rozenberg, M., Eilam, R., ... & Itzkovitz, S.
 (2018). Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis. Cell,
  175(4), 1156-1167.' 
  (a) table_A_LCM_TPM_values.tsv
  (b) table_D_zonation_reconstruction.tsv 
  and save them in an additional director 'in_vivo_villus_data'
3. at this point the directories should be in the following structure:
  HOME_DIR
  ---in_vivo_villus_data 
    --table_A_LCM_TPM_values.tsv
    --table_D_zonation_reconstruction.tsv
  ---- sprinkled 
    -- 12hr
        --roi1
        --roi2
        --roi3
    -- 72hr
       --roi1
       --roi2
       --roi3
       --roi4
 ---- unperturbed
        --- monolayer_erosion
4. change the directory path for HOME_DIR to where you have saved these directories in the format explained above 
"""

##CHANGE the directory path to where you have saved the directories above
HOME_DIR = 'C:/Users/micha/thesis/code/data/intestinal_organoid/repo_dir/'


#can change
SPRINKLED_DIR_NAME = 'sprinkled'
IN_VIVO_DIR_NAME = 'in_vivo_villus_data'
UNPERTURBED_DIR_NAME = 'unperturbed'
EROSION_DIR_NAME = 'monolayer_erosion'

#fits the structure explained at the top of the file. we suggest not to change
SPRINKLED_DIR = os.path.join(HOME_DIR,SPRINKLED_DIR_NAME)
IN_VIVO_VILLUS_DIR = os.path.join(HOME_DIR, IN_VIVO_DIR_NAME)
UNPERTURBED_DIR = os.path.join(HOME_DIR, UNPERTURBED_DIR_NAME)
EROSION_DIR = os.path.join(UNPERTURBED_DIR,EROSION_DIR_NAME)


##ROIS constants
TMPT_TO_ROIS_DICT = {'12hr': ['roi1', 'roi2', 'roi3'], '72hr':['roi1', 'roi2', 'roi3', 'roi4']}



NUM_POSITIONS_LCM_ATLAS = 5

VILLUS_LOCATION_COLUMNS = [ 'V1_mean', 'V2_mean', 'V3_mean', 'V4_mean', 'V5_mean', 'V6_mean']

AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH = os.path.join(os.getcwd(), 'paper','graphs','autonomous_zonation_figure_plots')
CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH = os.path.join(os.getcwd(), 'paper','graphs','continuous_regenerative_figure_plots')
ZONATION_PLASTICITY_PLOTS_FOLDER_PATH = os.path.join(os.getcwd(), 'paper','graphs','zonation_plasticity_plots')
NEIGHBORHOOD_ZONE_ADOPTION_FOLDER_PATH = os.path.join(os.getcwd(), 'paper','graphs','neighborhood_zone_adoption_plots')

#cell type constants
CELL_TYPES_TO_MARKER_GENES = {'enterocyte':['Alpi', 'Aldob','Sis','Apoa1'], 'goblet_cells':['Muc2'],'EEC':['Chga'],'Tuft_cells':['Dclk1'],'paneth_cells':['Lyz1'], 'stem_cells':['Lgr5','Olfm4'], 'regenerative':['Msln','Ahnak']}
CELL_TYPES_TO_MARKER_GENES_REG_RESPONSE = {'enterocyte':['Aldob','Ada','Apoa4','Apoa1','Alpi','Sis','Apob'], 'goblet_cells':['Muc2'],'EEC':['Chga'],'Tuft_cells':['Dclk1'],'paneth_cells':['Lyz1'], 'stem_cells':['Lgr5','Olfm4'], 'regenerative':['Mki67','Clu' ,'Yap1','Ly6a','Msln','Ahnak']}

# UNPERTURBED_ZOOMED_IN_X = [10000, 17000]
# UNPERTURBED_ZOOMED_IN_Y = [14000, 21000]

UNPERTURBED_ZOOMED_IN_X = [10000, 17000]
UNPERTURBED_ZOOMED_IN_Y = [15000, 22000]


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
PIXEL2NM = 107.11
EROSION_RINGS_ZOOM_IN = [1500, 2000]



MID_X_RANGE = [9219, 17780]
MID_Y_RANGE = [14269, 20986]


ZONATION_EXAMPLE_ZONE_X = [10000,17000]
ZONATION_EXAMPLE_ZONE_Y = [14000,21000]
TOP_GENES = ['Ada','Apoa4','Apoa1']
BOTTOM_GENES = ['Sis','Alpi']

COORDINATES = 'spatial'
X_COORDINATES = 'center_x'
Y_COORDINATES = 'center_y'
EXPECTED_ZONE_TB = 'transcript_exp_pos' #the expected zone of the cell along the transcript profile axis


GFP_12HR_ROI1_X = [18000, 22000]
GFP_12HR_ROI1_Y = [22000, 26000]


GFP_72HR_ROI2_X = [4500, 8000]
GFP_72HR_ROI2_Y = [22000, 28000]



CELL_TYPE_ANNOT_ORDER = ['regenerative', 'enterocyte','EEC','goblet_cells','stem_cells','paneth_cells','Tuft_cells']
ALL_CELL_TYPES_PLOT_ORDER =['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells','enterocyte','regenerative']
DOMINANT_CELL_TYPES = ['regenerative', 'enterocyte']
CELl_TYPE_WO_ENT_REG_PLOT_ORDER = ['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells']
CELl_TYPEs_WO_REGENERATIVE = ['EEC','Tuft_cells','paneth_cells','goblet_cells','stem_cells','enterocyte']
##gene constants


ORGANOID_GENE_NAMES = np.array([
    "Nupr1", "Ahnak", "Lypd8", "Ier3", "Pmepa1", "Slc12a2", "Apob",
    "Junb", "Lyz1", "Pigr", "Anxa13", "Mki67", "Olfm4", "S100g",
    "Anxa10", "Sis", "Prxl2a", "Basp1", "Sqstm1", "Txndc5", "Sprr1a",
    "Ccdc71l", "Ppp1r1b", "Rhoc", "Atf5", "Serpinb9b", "Ccna2", "Cdca7",
    "Il1rn", "Kcne3", "Apoc2", "Cryab", "Clps", "Plat", "Tac1", "Nlrp6",
    "Sult6b2", "Slc5a1", "Gkn3", "Egfr", "Selenom", "Klf4", "Sptssb",
    "Apoa4", "Il18", "Rnase1", "Lgr5", "Reg3g", "Yap1", "Hepacam2",
    "Slc2a2", "Insm1", "Neurog3", "Neurod1", "Vim", "Reg3b", "Pcsk1",
    "Rab3c", "Pclo", "Slc28a2", "Cck", "Msln", "Fos", "Tm4sf20", "Anxa5",
    "Clca3b", "Plaur", "Hook1", "Pycard", "Pclaf", "Ccn2", "Aldh1b1",
    "Gstm3", "Zfp36l2", "Fcgbp", "Jaml", "Cps1", "Smim24", "Ndrg1",
    "Ccl25", "Chga", "Smoc2", "Btc", "Tuba1a", "Apoa1", "Adh6a", "Alpi",
    "S100a7a", "Clca1", "Cavin3", "Rgcc", "Mal", "Add3", "Muc2", "Axin2",
    "Myb", "Sct", "Cyp2c29", "Ly6d", "Lor", "Slc7a7", "Mmp7", "Tph1",
    "H4c9", "Gstt1", "Guca2a", "Cpe", "Apobec1", "Ada", "Slc7a9", "Scg2",
    "Npc1l1", "Slc2a5", "Nt5e", "Reg3a", "Slc7a8", "Dclk1", "Slc15a1",
    "Ang4", "Gstm1", "Ly6a", "Rps23", "Itln1", "Reg1", "Agr2", "Tff3",
    "Prap1", "Spink4", "Rbp2", "Zg16", "Fabp1", "Reg4", "Tm4sf4",
    "Anxa1", "Lgals3", "F3", "Clu", "Chgb", "Aldob", "GFP"])




ORGANOID_GENE_NAMES_NOGFP = ORGANOID_GENE_NAMES[ORGANOID_GENE_NAMES !='GFP']

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

TRANSCRIPT_DENSITY_GENES = ['Anxa5','Apoa4','Apoa1','Alpi','Aldob','Sis']


ZONE_MAPPING_GENES = ['Ada', 'Apoa4', 'Apoa1', 'Aldob', 'Alpi', 'Sis']





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


LOW_ENT_RANGE = [1.6,1.8]
MEDIUM_ENT_RANGE = [1.8,2.0]
HIGH_ENT_RANGE = [2.1,2.3]
VERY_HIGH_ENT_RANGE = [2.4,2.5]
UNIFORM_DIST_RANGE = [2.5,2.6]
ENTROPY_RANGES = [LOW_ENT_RANGE, HIGH_ENT_RANGE, VERY_HIGH_ENT_RANGE, UNIFORM_DIST_RANGE]


#GFP cells
SPC_GFP_THRESH = 2

def set_style():
 plt.rcParams.update(plt.rcParamsDefault)  # Reset to default

 plt.rcParams.update({
     'figure.titlesize': 8, 'figure.titleweight': 'bold',
     'axes.titlesize': 8, 'axes.titleweight': "bold",
     'axes.labelsize': 8, 'axes.labelweight': 'bold',
     "ytick.labelsize": 8, "xtick.labelsize": 8,
     'legend.fontsize': 8, 'font.family': 'Arial',
     'figure.figsize': (4, 3), 'savefig.dpi': 300,
     "pdf.fonttype": 42,  # Keep text editable in PDFs
     "ps.fonttype": 42,  # Keep text editable in PS files
     "text.usetex" : False
 })

