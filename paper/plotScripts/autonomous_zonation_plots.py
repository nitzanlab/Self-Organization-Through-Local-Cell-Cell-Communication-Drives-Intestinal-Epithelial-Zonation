"""
This python file contains functions relevant for the figure which depicts that
2D intestinal organoids (monolayer) show autonomous zonated gene expression patterns, that
match those observed in the canonical zonation genes in in-vivo intestine
"""
from paper.extractedData.load_csvs import *
from paper.plotScripts.erosion_calculations import *
def plot_all_autonomous_figure_plots():
    ### panel a is schematic
    ###panel b : #TODO:Yael


    ###panel c: #Todo:Yael

    ###panel d: ##Todo:Yael

    ###panel e: bottom/(bottom+top) villus gene expression pattern
    #plot_gene_groups_expression_on_wt_monolayer(['Ada','Apoa4','Apoa1'], ['Sis','Alpi'], 'top', 'bottom')

    ###panel f: erosion rings
    #result_dict = plot_erosion_rings()

    ###panel g: erosion expression profiles
    #result_dict = load_transcript_densities_unperturbed_monolayer()
    #plot_density_profiles(result_dict, genes_in_order_density_measure, spread_plots=True)

    #panel h: heatmap comparison in vivo reconstruction to monolayer, erosion measured expression from edge to interior
    #plot_wt_monolayer_gene_density_to_invivo_comparisons()
    expression_profile_heatmap_comparison_invivo_gene_density(genes_in_order_density_measure)



def plot_gene_groups_expression_on_wt_monolayer(gene_group1, gene_group2, group1_name, group2_name):
    adata = get_unperturbed_monolayer_adata()
    plt.figure(figsize=(8, 6))
    gene_exp = adata[:, gene_group1].X.mean(axis=1) / (
                adata[:, gene_group1].X.mean(axis=1) + adata[:, gene_group2].X.mean(axis=1))
    signal_df = pd.DataFrame(adata.X, columns=adata.var.index,  index=adata.obs_names)

    signal_df['x'] = adata.obsm['spatial']['center_x']
    signal_df['y'] = adata.obsm['spatial']['center_y']
    signal_df['signal'] = gene_exp
    signal_df.to_csv(os.path.join(WT_MONOLAYER_DIR, 'top_bottom_villus_exp.csv'))

    sctr1 = plt.scatter(adata.obsm['spatial']['center_x'], adata.obsm['spatial']['center_y'],
                        c=gene_exp, s=40,  cmap='viridis')

    plt.xlim(10000,17000)
    plt.ylim(14000,21000)
    plt.gca().invert_yaxis()
    # plt.xticks()
    # plt.yticks()
    plt.title(f'{group1_name}/{group2_name} villus expression in unperturbed monolayer')
    cbar = plt.colorbar(sctr1)
    cbar.set_label(f'{group1_name}/{group2_name}\n expression')
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'bottom_top_villus_exp.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()

def expression_profile_heatmap_comparison_invivo_gene_density(gene_set):
    reconstruction = mean_gene_exp_per_zone_in_invivo_reconstruction_no_crypt()
    result_dict = load_unperturbed_monolayer_gene_densities()
    gene_density_df = pd.DataFrame({gene: values['densities'] for gene, values in result_dict.items()}).iloc[:20]
    gene_set_monolayer = [gene for gene in gene_set if gene in ORGANOID_GENE_NAMES_NOGFP] #only genes that are panneled
    invivo_raw_exp = load_TPM_LCM_intestine_atlas()
    invivo_exp, genes_LCM = get_LCM_atlas_gene_subset(invivo_raw_exp, gene_set_monolayer)
    genes_invivo_and_monolayer = [gene for gene in gene_set_monolayer if gene in invivo_exp.index]
    genes_invivo_monolayer_reconstruction = [gene for gene in genes_invivo_and_monolayer if gene in reconstruction.index]
    gene_density_df = gene_density_df[gene_set]
    reconstruction_df = reconstruction.loc[genes_invivo_monolayer_reconstruction]
    ##normalize for the comparison
    gene_density_df = gene_density_df.reset_index(drop=True).iloc[4:16,:] #was
    gene_density_normalized = ((gene_density_df - gene_density_df.min(axis=0)) / (gene_density_df.max(axis=0) - gene_density_df.min(axis=0))).T
    invivo_exp_normalized = (invivo_exp.T - invivo_exp.T.min(axis=0)) / (invivo_exp.T.max(axis=0) - invivo_exp.T.min(axis=0))
    transcript_heatmap = sns.heatmap(gene_density_normalized.iloc[:,:-1], vmin=0, vmax=1, cmap='plasma')
    colorbar = transcript_heatmap.collections[0].colorbar
    colorbar.set_label('normalized gene expression')
    plt.xlabel('iterations')
    plt.xticks(rotation=0)
    plt.ylabel('genes')
    plt.yticks()
    plt.title('transcript density profiles')
    plt.tight_layout()
    plt.show()

    sns.heatmap(invivo_exp_normalized.T, vmin=0, vmax=1, cmap='plasma')
    plt.title('invivo expression profiles')
    plt.tight_layout()
    plt.show()

    reconstruction_normalized = (reconstruction_df.T - reconstruction_df.T.min(axis=0)) / (
                reconstruction_df.T.max(axis=0) - reconstruction_df.T.min(axis=0))
    reconstruction_normalized = reconstruction_normalized = reconstruction_normalized.T.iloc[:, ::-1]
    reconstruction_normalized.rename(columns={f'V{i}_mean': f'V{i}' for i in range(1, 7)}, inplace=True)
    recon_heatmap = sns.heatmap(reconstruction_normalized, cmap='plasma')
    colorbar = recon_heatmap.collections[0].colorbar
    colorbar.set_label('normalized gene expression')
    plt.xlabel('villus top to bottom', fontsize=16)
    plt.xticks(rotation=0, fontsize=14)
    plt.ylabel('genes', fontsize=16)
    plt.yticks(fontsize=14)
    plt.title('invivo reconstruction')
    plt.tight_layout()
    plt.show()

def mean_gene_exp_per_zone_in_invivo_reconstruction_no_crypt():
    reconstruction = load_invivo_reconstruction()
    return reconstruction[['V1_mean', 'V2_mean', 'V3_mean', 'V4_mean', 'V5_mean', 'V6_mean']]