def normalize_cell_by_gene_by_cells_then_genes(cell_by_gene:pd.DataFrame, genes):
    cell_by_gene_subset = ad.AnnData(cell_by_gene[genes])
    sc.pp.normalize_total(cell_by_gene_subset, target_sum=1e4)
    normalized_df = pd.DataFrame(cell_by_gene_subset.X, columns=cell_by_gene[genes].columns)
    return normalized_df