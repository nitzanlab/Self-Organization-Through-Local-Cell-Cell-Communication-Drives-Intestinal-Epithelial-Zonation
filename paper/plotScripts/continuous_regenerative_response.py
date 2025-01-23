from utils.imports import *
from utils.constant import *
from paper.extractedData.load_csvs import *
from paper.extractedData.data_preprocessing import *
def plot_all_continuous_regenerative_response_plots():
    ###panel a: schematic diagram created in https://BioRender.com
    ###panel b: raw image examples #TODO Yael

    ###panel c: L metric analysis #TODO Yael

    ###panel d: heatmap L metric in unperturbed data #TODO Yael

    ###panel e: regenerative expression in unperturbed monolayer
    #plot_regenerative_expression_across_cell_types_spatially('Msln')

    #plot_regenerative_expression_across_cell_types_spatially('Msln', x_region=UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2, y_region=UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2)
    for gene in ENTEROCYTE_GENES:
        plot_regenerative_expression_across_cell_types_spatially(gene, x_region=UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2,
                                                             y_region=UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2)

    #expression of secretory and progenitor genes
    # genes = ['Mki67'] #['Muc2','Chga','Dclk1','Lyz1','Mki67']
    # for gene in genes:
    #     plot_regenerative_expression_across_cell_types_spatially(gene, x_region=UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2, y_region=UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2)

    ###panel f: zoom in regions


    ###panel g: regenerative expression neighborhood correlation
    #plot_regenerative_gene_expression_neighborhood_similarity(ALL_CELL_TYPE_GENES)


def plot_regenerative_expression_across_cell_types_spatially(goi, x_region=UNPERTURBED_ZOOMED_IN_X,y_region=UNPERTURBED_ZOOMED_IN_Y,title=''):
    cell_by_gene = load_unperturbed_intestinal_organoid_cell_by_gene_mat()
    cell_by_gene_normed = normalize_cell_by_gene_by_cells_then_genes(cell_by_gene, ORGANOID_GENE_NAMES_NOGFP)
    cell_coords = load_unperturbed_cell_coords()
    clusters_df = pd.read_csv(os.path.join(DATA_DIR, 'cell_by_gene_cluster_annotations.csv'))
    clusters = clusters_df.sort_values(by='object_id', ascending=True)['cluster_id']
    unique_clusters = np.unique(clusters)
    num_clusters = len(unique_clusters)
    # Create a colormap for cluster edges
    edge_colormap = cm.get_cmap('tab10', num_clusters)  # Use 'tab20' for discrete colors
    edge_colors = [
        edge_colormap(cluster / num_clusters) if cluster in [4, 5] else 'none'
        for cluster in clusters
    ]
    # edge_colors = [
    #     edge_colormap(cluster / num_clusters)
    #     for cluster in clusters
    # ]

    # normalized_exp = normalize_gene_exp_for_tissue_rep(cell_by_gene_normed[goi])
    plt.scatter(cell_coords['center_x'], cell_coords['center_y'], c=cell_by_gene_normed[goi], cmap='Greens', s=20,
                edgecolors=edge_colors)
    cbar = plt.colorbar()
    cluster_dict = {0: 'regenerative', 1: 'regenerative', 2: 'regenerative', 3: 'enterocyte', 4: 'secretory',
                    5: 'progenitor', 6: 'regenerative', 7: 'unknown'}
    for cluster in [4, 5]:  # unique_clusters:
        plt.scatter([], [], edgecolor=edge_colormap(cluster / num_clusters), facecolor='none',
                    label=f'{cluster_dict[cluster]}', linewidth=1.5, s=100)

    plt.legend(title='Clusters', loc='upper left', fontsize='small', title_fontsize='medium')

    cbar.set_label(f'{goi} expression levels')
    plt.xlim(x_region[0], x_region[1])
    plt.ylim(y_region[0], y_region[1])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xticks([])
    plt.yticks([])
    plt.title(f'{goi} expression in tissue wt {title}')
    plt.gca().invert_yaxis()
    plt.show()



def plot_regenerative_gene_expression_neighborhood_similarity(genes):
    adata = get_unperturbed_monolayer_adata()
    clusters_df = pd.read_csv(os.path.join(DATA_DIR, 'cell_by_gene_cluster_annotations.csv'))
    clusters = clusters_df.sort_values(by='object_id', ascending=True)['cluster_id']
    unique_clusters = np.unique(clusters)

    cell_type_corrs = {}
    for cell_type_cluster in unique_clusters[:-1]:
        all_genes_one_type = []
        cell_indices = np.where(clusters == cell_type_cluster)[0]
        cell_neighbors_idx = adata.obsm['neighbors_idx'][cell_indices]
        for gene in genes:
            mean_neigh_exp = adata[:,gene].X[cell_neighbors_idx].mean(axis=1).flatten()
            corr = np.corrcoef(adata[:,gene].X[cell_indices].flatten(), mean_neigh_exp)[0, 1]
            all_genes_one_type.append(corr)
        cell_type_corrs[UNPERTURBED_CELL_TYPE_CLUSTERS[cell_type_cluster]] = np.array(all_genes_one_type)
    all_cell_types_corrs = pd.DataFrame(cell_type_corrs, index=genes).T
    all_cell_types_corrs.fillna(0, inplace=True)
    print(all_cell_types_corrs.columns)
    genes_enterocyte = ['Ada', 'Apoa4','Apoa1','Alpi','Sis','Aldob']
    genes_reg = ['Clu','Msln','Ahnak']

    data = [all_cell_types_corrs.loc['secretory'][genes_enterocyte], all_cell_types_corrs.loc['enterocyte'][genes_enterocyte],
            all_cell_types_corrs.loc['secretory'][genes_reg],
            all_cell_types_corrs.loc['enterocyte'][genes_reg]]
    group_titles = ['entercoyte gene correlation', 'regenerative gene correlation']
    bar_titles = ['secretory cells','enterocyte cells']
    means = [np.mean(lst) for lst in data]
    stds = [np.std(lst) for lst in data]

    # Group data for plotting
    groups = [0, 0, 1, 1]  # 0 for Group 1, 1 for Group 2
    bar_positions = np.arange(len(data))  # Bar positions
    width = 0.35  # Width of each bar

    # Create the bar plot
    fig, ax = plt.subplots(figsize=(6,4))

    bars = ax.bar(bar_positions, means, yerr=stds, capsize=5, width=width, color=['skyblue', 'lightgreen'])

    # Add bar titles
    for i, bar in enumerate(bars):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.05,
                f'{bar_titles[i % 2]}', ha='center', va='bottom', fontsize=10)

    # Add group titles
    group_ticks = [(bar_positions[groups == g].mean()) for g in np.unique(groups)]
    ax.set_xticks(group_ticks)
    ax.set_xticklabels(group_titles,fontsize=10)

    ax.set_ylabel('Mean Correlation',fontsize=10)
    ax.set_title('Mean Expression Correlations to Neighboring cells')
    plt.tight_layout()
    plt.show()

    return all_cell_types_corrs