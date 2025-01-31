"""
This python file contains functions relevant for the figure which depicts that
2D intestinal organoids (monolayer) show autonomous zonated gene expression patterns, that
match those observed in the canonical zonation genes in in-vivo intestine
"""
from paper.extractedData.load_csvs import *
from paper.plotScripts.erosion_calculations import *
def plot_all_autonomous_figure_plots(saved_datasets=False):
    ### panel a is schematic ,created in Biorender.com
    ###panel b : #TODO:Yael


    ###panel c: #Todo:Yael

    ###panel d: ##Todo:Yael

    ###panel e: bottom/(bottom+top) villus gene expression pattern
    #plot_gene_groups_expression_on_wt_monolayer(TOP_GENES, BOTTOM_GENES, 'top', 'bottom', zoned=True,
    #                                            zone_x=ZONATION_EXAMPLE_ZONE_X, zone_y=ZONATION_EXAMPLE_ZONE_Y, save=False)

    ###panel f: erosion rings

    #if the erosion was already conducted, saved_datasets can be set to True, and thus plot from loaded data
    # from pickle, otherwise, calculate from scratch and save
    if saved_datasets:
         plot_erosion_rings_from_saved_components()
         result_dict = load_transcript_densities_unperturbed_monolayer()
    else:
         result_dict = calculate_eroded_transcription_densities(plot_erosion=True, save_components=False)
    ###panel g: erosion expression profiles
    plot_density_profiles(result_dict, genes_in_order_density_measure, spread_plots=True, normalize=True)
    #panel h: heatmap comparison in vivo reconstruction to monolayer, erosion measured expression from edge to interior
    expression_profile_heatmap_comparison_invivo_gene_density(genes_in_order_density_measure)

def plot_gene_groups_expression_on_wt_monolayer(gene_group1:list, gene_group2:list, group1_name:str, group2_name:str,zoned=False,zone_x=None, zone_y=None, save=False):
    """
    This function plots and saves the relative expression of one gene group to another over the unperturbed monolayer
    :param gene_group1: a list of genes in group 1
    :param gene_group2: a list of genes in group 2
    :param group1_name: title of group 1
    :param group2_name: title of group 2
    :param zoned: boolean : if to plot a subregion of the monolayer
    :param zone_x: the xrange of the subregion to plot (list with 2 values)
    :param zone_y: the y range of the subregion (list with 2 values)
    :param save: booelan: if to save the spatial signal as a csv
    """
    adata = get_unperturbed_monolayer_adata()
    #normalize to relative expression between the gene groups based on the mean expression of each group
    gene_exp = adata[:, gene_group1].X.mean(axis=1) / (
                adata[:, gene_group1].X.mean(axis=1) + adata[:, gene_group2].X.mean(axis=1))\

    #save spatial signal as a pandas dataframe
    adata.obs['relative_exp'] = gene_exp
    save_spatial_signal(adata, 'relative_exp', f'{gene_group1}_{gene_group2}_spatial_expression', sprinkled=False)

    sctr1 = plt.scatter(adata.obsm[COORDINATES][X_COORDINATES], adata.obsm[COORDINATES][Y_COORDINATES],
                        c=gene_exp, s=40,  cmap='viridis')

    if zoned:
        plt.xlim(zone_x[0], zone_x[1])
        plt.ylim(zone_y[0],zone_y[1])

    plt.gca().invert_yaxis()
    plt.title(f'{group1_name}/{group2_name} villus expression in unperturbed monolayer')
    cbar = plt.colorbar(sctr1)
    cbar.set_label(f'{group1_name}/{group2_name}\n expression')
    plt.tight_layout()
    if save:
        os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
        file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, f'{gene_group1}_{gene_group2}_spatial_expression.pdf')
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
    colorbar.set_label('Normalized Gene Expression')
    avg_iteration_width = load_iteration_average_width()
    dist_to_edge = np.round(np.arange(4, 4+gene_density_df.shape[0]) * avg_iteration_width, 2)
    plt.xticks(np.arange(gene_density_df.shape[0]), dist_to_edge)
    plt.xlabel('Distance to Monolayer Edge(μm)')
    plt.xticks(rotation=45)
    plt.ylabel('Genes')
    plt.yticks()
    plt.title('Transcript Density Profiles')
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'monolayer_transcript_density_expression_heatmap.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()

    # sns.heatmap(invivo_exp_normalized.T, vmin=0, vmax=1, cmap='plasma')
    # plt.title(r"$\it{In\ Vivo}$ Expression Profiles")
    # plt.tight_layout()
    # plt.show()

    reconstruction_normalized = (reconstruction_df.T - reconstruction_df.T.min(axis=0)) / (
                reconstruction_df.T.max(axis=0) - reconstruction_df.T.min(axis=0))
    reconstruction_normalized = reconstruction_normalized = reconstruction_normalized.T.iloc[:, ::-1]
    reconstruction_normalized.rename(columns={f'V{i}_mean': f'V{i}' for i in range(1, 7)}, inplace=True)
    recon_heatmap = sns.heatmap(reconstruction_normalized, cmap='plasma')
    colorbar = recon_heatmap.collections[0].colorbar
    colorbar.set_label('Normalized Fene expression')
    plt.xlabel('Villus Top to Bottom')
    plt.xticks(rotation=0)
    plt.ylabel('Genes')
    plt.yticks()
    plt.title(r"$\it{In\ Vivo}$ Expression Profiles Reconstructed")
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'invivo_reconstruction_expression_heatmap.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()

def mean_gene_exp_per_zone_in_invivo_reconstruction_no_crypt():
    reconstruction = load_invivo_reconstruction()
    return reconstruction[['V1_mean', 'V2_mean', 'V3_mean', 'V4_mean', 'V5_mean', 'V6_mean']]