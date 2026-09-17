"""
This python file contains functions relevant for the figure which depicts that
2D intestinal organoids (monolayer) show autonomous zonated gene expression patterns, that
match those observed in the canonical zonation genes in in-vivo intestine
"""
from paper.extractedData.load_csvs import *
from paper.plotScripts.erosion_calculations import *
def plot_all_autonomous_figure_plots(saved_datasets=False):
    """

    :param saved_datasets: if already ran with saved_datasets=True, then the necessary files
    exist in the directory and can be used for calculations and plots
    :return:
    """
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
         result_dict = calculate_eroded_transcription_densities(plot_erosion=True, save_components=True)
    ###panel g: erosion expression profiles
    plot_density_profiles(result_dict, genes_in_order_density_measure, spread_plots=True, normalize=True)
    #panel h: heatmap comparison in vivo reconstruction to monolayer, erosion measured expression from edge to interior
    expression_profile_heatmap_comparison_invivo_gene_density(genes_in_order_density_measure)
    from utils.constant import PHARMACOLOGICAL_PERTURBATIONS_PLOTS_FOLDER_PATH
    expression_profile_heatmap_invivo(EPHRIN_GENES,
                                      save_name="in_vivo_epha2_expression.pdf",
                                      output_dir=PHARMACOLOGICAL_PERTURBATIONS_PLOTS_FOLDER_PATH)

    from paper.plotScripts.crop_monolayer import plot_monolayer_raw_figures
    from paper.plotScripts.top_bottom_zonation import plot_top_bottom_villus_expression
    from paper.plotScripts.cluster_monolayer import plot_monolayer_cluster_identity

    # Panels B/D, E and C are independent: one missing input must not drop the others.
    remaining = [
        ("B,D: raw monolayer crops", plot_monolayer_raw_figures),
        ("E: top/bottom villus expression", plot_top_bottom_villus_expression),
        ("C: monolayer cluster identity", plot_monolayer_cluster_identity),
    ]
    failed = []
    for label, fn in remaining:
        print("\nGenerating Panel %s …" % label)
        try:
            fn()
        except Exception as exc:
            failed.append((label, exc))
            print("  SKIPPED Panel %s — %s: %s" % (label, type(exc).__name__, exc))
    if failed:
        import warnings
        for label, exc in failed:
            warnings.warn("Figure 1 Panel %s not generated (%s: %s)"
                          % (label, type(exc).__name__, exc), RuntimeWarning)
        print("\n%d of %d Figure 1 sub-panels could not be generated:" % (len(failed), len(remaining)))
        for label, exc in failed:
            print("  - Panel %s (%s)" % (label, type(exc).__name__))

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
    plt.close()

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
    dist_to_edge = np.round(
        (np.arange(4, 4 + gene_density_df.shape[0]) + 0.5) * RING_STEP_UM, 0).astype(int)
    plt.xticks(np.arange(gene_density_df.shape[0]), dist_to_edge)
    plt.xlabel('Distance to Monolayer Edge(μm)')
    plt.xticks(rotation=45)
    plt.ylabel('Genes')
    plt.yticks()
    plt.yticks(rotation=0)
    plt.title('Transcript Density Profiles')
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'monolayer_transcript_density_expression_heatmap.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()
    plt.close()

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
    plt.yticks(rotation=0)
    plt.title(r"$\it{In\ Vivo}$ Expression Profiles Reconstructed")
    plt.tight_layout()
    os.makedirs(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'invivo_reconstruction_expression_heatmap.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()
    plt.close()
    
def expression_profile_heatmap_invivo(
    gene_set,
    title=r"$\it{In\ Vivo}$ Expression",
    save_name="invivo_expression_heatmap.pdf",
    flip_top_to_bottom=True,
    output_dir=None,
):
    """
    Plot a heatmap of in vivo reconstruction ONLY (no monolayer checks).

    Parameters
    ----------
    gene_set : iterable of str
        Genes you'd like to plot.
    title : str
        Plot title.
    save_name : str
        PDF filename saved under AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH.
    flip_top_to_bottom : bool
        If True, reverse V1_mean..V6_mean so the x-axis reads Top→Bottom.

    Returns
    -------
    plotted_genes : list[str]
        The genes that were present in the reconstruction and plotted.
    """
    # Load in vivo reconstruction (expects V1_mean..V6_mean columns)
    reconstruction = mean_gene_exp_per_zone_in_invivo_reconstruction_no_crypt()
    recon_index = pd.Index(reconstruction.index).astype(str)

    # Keep only genes that exist in reconstruction
    requested = list(pd.Index(gene_set).astype(str))
    present = [g for g in requested if g in recon_index]
    dropped = [g for g in requested if g not in recon_index]
    if dropped:
        print(f"[INFO] {len(dropped)} genes not in reconstruction (ignored): {dropped[:10]}")
    if not present:
        print("[INFO] No requested genes found in reconstruction. Nothing to plot.")
        return []

    # Subset and keep only numeric columns (guards against stray metadata)
    df = reconstruction.loc[present]
    num_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if not num_cols:
        raise ValueError("Reconstruction table has no numeric columns to plot.")
    df = df[num_cols]

    # Optional flip: Top←→Bottom
    if flip_top_to_bottom:
        df = df.iloc[:, ::-1]

    # Nicer column labels if they look like V*_mean
    df = df.rename(columns={c: c.replace("_mean", "") for c in df.columns})

    # Normalize 0–1 per gene (row-wise)
    norm = (df.T - df.T.min(axis=0)) / (df.T.max(axis=0) - df.T.min(axis=0))
    norm = norm.T

    # Plot
    plt.figure()
    ax = sns.heatmap(norm, cmap='plasma', vmin=0, vmax=1)
    cbar = ax.collections[0].colorbar
    cbar.set_label('Normalized Gene Expression')
    plt.xlabel('Villus Top to Bottom' if flip_top_to_bottom else 'Villus Bottom to Top')
    plt.xticks(rotation=0)
    plt.ylabel('Genes')
    plt.yticks(rotation=0)
    plt.title(title)
    plt.tight_layout()

    # Save
    from utils.constant import PHARMACOLOGICAL_PERTURBATIONS_PLOTS_FOLDER_PATH
    _out_dir = output_dir if output_dir else AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH
    os.makedirs(_out_dir, exist_ok=True)
    out_path = os.path.join(_out_dir, save_name)
    plt.savefig(out_path, format='pdf')
    plt.show()
    plt.close()
    
def mean_gene_exp_per_zone_in_invivo_reconstruction_no_crypt():
    reconstruction = load_invivo_reconstruction()
    return reconstruction[['V1_mean', 'V2_mean', 'V3_mean', 'V4_mean', 'V5_mean', 'V6_mean']]