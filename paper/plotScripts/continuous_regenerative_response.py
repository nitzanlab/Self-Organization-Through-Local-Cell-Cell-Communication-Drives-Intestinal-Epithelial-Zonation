import pandas as pd

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
    plot_gene_expression_across_cell_types_spatially(['Msln','Aldob'],colors=['Greys','Purples'], title='Msln_Aldob_expression_unperturbed')
    #plot_gene_expression_across_cell_types_spatially('Aldob', color='Purples')

    #plot_regenerative_expression_across_cell_types_spatially('Msln', x_region=UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2, y_region=UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2)
    # for gene in ENTEROCYTE_GENES:
    #     plot_regenerative_expression_across_cell_types_spatially(gene,color='Purples')

    #expression of secretory and progenitor genes
    # genes = ['Mki67'] #['Muc2','Chga','Dclk1','Lyz1','Mki67']
    # for gene in genes:
    #     plot_regenerative_expression_across_cell_types_spatially(gene, x_region=UNPERTURBED_ZOOMED_IN_X_SEC_CELL_2, y_region=UNPERTURBED_ZOOMED_IN_Y_SEC_CELL_2)

    ###panel f: zoom in regions


    ###panel g: regenerative expression neighborhood correlation
    calculate_regenerative_gene_expression_neighborhood_similarity(ALL_CELL_TYPE_GENES)


def plot_gene_expression_across_cell_types_spatially(gois,colors, x_region=UNPERTURBED_ZOOMED_IN_X, y_region=UNPERTURBED_ZOOMED_IN_Y, title=''):
    adata = get_unperturbed_monolayer_adata()
    clusters_df = pd.read_csv(os.path.join(UNPERTURBED_DIR, 'cell_by_gene_cluster_annotations.csv'))
    clusters = clusters_df.sort_values(by='object_id', ascending=True)['cluster_id']
    adata.obs['cluster'] = clusters.values

    cluster_dict = {
        0: 'regenerative', 1: 'regenerative', 2: 'regenerative', 3: 'enterocyte',
        4: 'secretory', 5: 'progenitor', 6: 'regenerative', 7: 'unknown'
    }
    adata.obs['cluster_name'] = adata.obs['cluster'].map(cluster_dict)
    signal_df = pd.DataFrame(adata[:,gois].X,columns=gois)
    signal_df['x'] = np.array(adata.obsm[COORDINATES][X_COORDINATES])
    signal_df['y'] = np.array(adata.obsm[COORDINATES][Y_COORDINATES])
    # Add cluster type
    signal_df['cluster_type'] = np.array(adata.obs['cluster_name'])
    signal_df.to_csv(os.path.join(UNPERTURBED_DIR, f'{title}.csv'))
    # Edge color customization for clusters 4 and 5
    edge_colors = [
        "#FFD700" if cluster == 4 else "#FF4500" if cluster == 5 else 'none'  # Gold for cluster 4, vibrant orange-red for cluster 5
        for cluster in clusters
    ]
    for i,goi in enumerate(gois):
        scatter = plt.scatter(
            adata.obsm[COORDINATES][X_COORDINATES], adata.obsm[COORDINATES][Y_COORDINATES],
            c=adata[:,goi].X, cmap=colors[i], s=20,
            edgecolors=edge_colors  # Apply edge color only for clusters 4 and 5
        )

        # Add color bar for expression values
        cbar = plt.colorbar(scatter)
        cbar.set_label(f'{goi} expression levels')

        # Cluster information for the legend


        # Manually add legend for clusters 4 and 5 with distinct colors
        legend_labels = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#FFD700", markersize=10,
                       label='secretory'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF4500", markersize=10,
                       label='progenitor')
        ]

        plt.legend(handles=legend_labels, title='Clusters', loc='upper left', fontsize='small', title_fontsize='medium')

        # Set axis limits, labels, and other plot properties
        plt.xlim(x_region)
        plt.ylim(y_region)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xticks([])
        plt.yticks([])
        plt.title(f'{goi} expression in tissue wt {title}')
        plt.gca().invert_yaxis()
        # Save the plot
        os.makedirs(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, exist_ok=True)
        file_name = os.path.join(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, f'{goi}_spatial_expression.pdf')
        plt.savefig(file_name, format='pdf')
        plt.show()

# def plot_gene_expression_across_cell_types_spatially(goi, x_region=UNPERTURBED_ZOOMED_IN_X, y_region=UNPERTURBED_ZOOMED_IN_Y, title='', color='Oranges'):
#     cell_by_gene = load_unperturbed_intestinal_organoid_cell_by_gene_mat()
#     cell_by_gene_normed = normalize_cell_by_gene_by_cells_then_genes(cell_by_gene, ORGANOID_GENE_NAMES_NOGFP)
#     cell_coords = load_unperturbed_cell_coords()
#     clusters_df = pd.read_csv(os.path.join(DATA_DIR, 'cell_by_gene_cluster_annotations.csv'))
#     clusters = clusters_df.sort_values(by='object_id', ascending=True)['cluster_id']
#     #pinks_colormap = LinearSegmentedColormap.from_list("pinks", ["#ffd1dc", "#ff69b4", "#ff1493", "#c71585"])
#
#     # Edge color customization for clusters 4 and 5
#     edge_colors = [
#         "#FFD700" if cluster == 4 else "#FF4500" if cluster == 5 else 'none'  # Gold for cluster 4, vibrant orange-red for cluster 5
#         for cluster in clusters
#     ]
#     scatter = plt.scatter(
#         cell_coords['center_x'], cell_coords['center_y'],
#         c=cell_by_gene_normed[goi], cmap=color, s=20,
#         edgecolors=edge_colors  # Apply edge color only for clusters 4 and 5
#     )
#
#     # Add color bar for expression values
#     cbar = plt.colorbar(scatter)
#     cbar.set_label(f'{goi} expression levels')
#
#     # Cluster information for the legend
#     cluster_dict = {
#         0: 'regenerative', 1: 'regenerative', 2: 'regenerative', 3: 'enterocyte',
#         4: 'secretory', 5: 'progenitor', 6: 'regenerative', 7: 'unknown'
#     }
#
#     # Manually add legend for clusters 4 and 5 with distinct colors
#     legend_labels = [
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#FFD700", markersize=10,
#                    label='secretory'),
#         plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#FF4500", markersize=10,
#                    label='progenitor')
#     ]
#
#     plt.legend(handles=legend_labels, title='Clusters', loc='upper left', fontsize='small', title_fontsize='medium')
#
#     # Set axis limits, labels, and other plot properties
#     plt.xlim(x_region)
#     plt.ylim(y_region)
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     plt.xticks([])
#     plt.yticks([])
#     plt.title(f'{goi} expression in tissue wt {title}')
#     plt.gca().invert_yaxis()
#     # Save the plot
#     os.makedirs(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, exist_ok=True)
#     file_name = os.path.join(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, f'{goi}_spatial_expression.pdf')
#     plt.savefig(file_name, format='pdf')
#     plt.show()



def calculate_regenerative_gene_expression_neighborhood_similarity(genes:list)->pd.DataFrame:
    """
    This function calculates per each cell type, the correlation in expression between the cells of that type
    and the mean expression of their neighboring cells per gene given in list genes
    :param genes: the genes for which to calculate similarity in expression between each cell and the mean expression
    of its neighboring cells
    :return: a pd.Dataframe of shape num cell types X num genes, where each entry is the expression correlation of the specified gene
    of the cells of that type with its neighboring cells
    """
    #load unperturbed monolayer
    adata = get_unperturbed_monolayer_adata()

    clusters_df = pd.read_csv(os.path.join(UNPERTURBED_DIR, 'cell_by_gene_cluster_annotations.csv'))
    clusters = clusters_df.sort_values(by='object_id', ascending=True)['cluster_id']
    unique_clusters = np.unique(clusters)

    cell_type_corrs = {}
    for cell_type_cluster in unique_clusters[:-1]: #for each cell type cluster, get the correlation in the mean expression of each gene in genes of the neighboring cells to the cell
        all_genes_one_type = []
        cell_indices = np.where(clusters == cell_type_cluster)[0]  #gets the indices of all of the cells in the cell_type_cluster
        cell_neighbors_idx = adata.obsm['neighbors_idx'][cell_indices] ##the indices of the cells neighboring the cell type cells
        for gene in genes:
            mean_neigh_exp = adata[:,gene].X[cell_neighbors_idx].mean(axis=1).flatten() #mean expresssion over the neighboring cells of the given gene
            corr = np.corrcoef(adata[:,gene].X[cell_indices].flatten(), mean_neigh_exp)[0, 1]
            all_genes_one_type.append(corr)
        cell_type_corrs[UNPERTURBED_CELL_TYPE_CLUSTERS[cell_type_cluster]] = np.array(all_genes_one_type)
    all_cell_types_corrs = pd.DataFrame(cell_type_corrs, index=genes).T #pd dataframe num cell types X num genes
    all_cell_types_corrs.fillna(0, inplace=True)

    plot_regenerative_neighboring_expression_correlation(all_cell_types_corrs)
    return all_cell_types_corrs

def plot_regenerative_neighboring_expression_correlation(all_cell_types_corrs):
    """
    This function plots the correlation in expression between enterocyte cells and secretory cells with their
    neighboring cells of enterocyte genes and regenerative genes
    :param all_cell_types_corrs: a pd.Dataframe of shape num cell types X num genes, where each entry is the expression correlation of the specified gene
    of the cells of that type with its neighboring cells
    """
    genes_enterocyte = ['Ada', 'Apoa4', 'Apoa1', 'Alpi', 'Sis', 'Aldob']
    genes_reg = ['Clu', 'Msln', 'Ahnak']
    data = [all_cell_types_corrs.loc['secretory'][genes_enterocyte], all_cell_types_corrs.loc['enterocyte'][genes_enterocyte],
            all_cell_types_corrs.loc['secretory'][genes_reg],
            all_cell_types_corrs.loc['enterocyte'][genes_reg]]


    group_titles = ['Enterocyte \nGenes Expression \nCorrelation', 'Regenerative\n Gene Expression\n Correlation']
    bar_titles = ['Secretory\n Cells','Enterocyte\n Cells']
    means = [np.mean(lst) for lst in data]
    stds = [np.std(lst) for lst in data]

    # Group data for plotting
    groups = [0, 0, 1, 1]  # 0 for Group 1, 1 for Group 2
    bar_positions = np.arange(len(data))/2  # Bar positions
    width = 0.35  # Width of each bar

    # Create the bar plot
    fig, ax = plt.subplots(figsize=(4,3))

    bars = ax.bar(bar_positions, means, yerr=stds, capsize=5, width=width, color=['skyblue', 'lightgreen'])

    # Add bar titles
    label_size = plt.rcParams['axes.labelsize']
    for i, bar in enumerate(bars):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + stds[i]+0.01,
                f'{bar_titles[i % 2]}', ha='center', va='bottom', fontsize=label_size)

    # Add group titles
    group_ticks = [(bar_positions[groups == g].mean()) for g in np.unique(groups)]
    ax.set_xticks(group_ticks)
    ax.set_xticklabels(group_titles, fontweight='bold')
    ax.set_ylim(0,0.9)
    ax.set_ylabel('Mean Correlation')
    ax.set_title('Mean Expression Correlations to Neighboring cells')
    plt.tight_layout()
    os.makedirs(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, exist_ok=True)
    file_name = os.path.join(CONTINUOUS_REGENERATIVE_RESPONSE_PLOTS_FOLDER_PATH, f'regenerative_expression_across_cell_types.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()
