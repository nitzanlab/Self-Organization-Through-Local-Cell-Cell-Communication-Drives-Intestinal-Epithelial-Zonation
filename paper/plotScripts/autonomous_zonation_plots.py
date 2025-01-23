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
    plot_gene_groups_expression_on_wt_monolayer(['Ada','Apoa4','Apoa1'], ['Sis','Alpi'], 'top', 'bottom')

    ###panel f: erosion rings
    # plot_erosion_rings()
    # plot_zoomed_in_erosion_rings()

    ###panel g: erosion expression profiles

    ###panel h: heatmap comparison in vivo reconstruction to monolayer, erosion measured expression from edge to interior
    pass


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
    file_name = os.path.join(AUTONOMOUS_ZONATION_PLOTS_FOLDER_PATH, 'bottom_top_villus_exp.pdf')
    plt.savefig(file_name, format='pdf')
    plt.show()

