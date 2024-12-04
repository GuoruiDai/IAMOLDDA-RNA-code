import os
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet
import torch
import random
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def load_and_filter_data(h5_file_path, demuxlet_file_path):
    # Load aggregated h5 data and process cell barcode suffix
    adata = sc.read_10x_h5(h5_file_path)
    adata.obs.index = adata.obs.index.str.replace('-2$', '-1', regex=True) # change second set of barcode with suffix -1
    adata = adata[~adata.obs.index.duplicated(keep='last')]
    assert not adata.obs.index.duplicated().any() # assert no duplicate barcodes

    # swap var_name with gene_ids column such that, ensembl id is used as var_name, and mixed gene id is under gene_ids
    _adata = adata.copy()
    _adata.var_names_make_unique()  # remove duplicates in the mixed gene id column first
    _converted_gene_id_array = _adata.var_names
    _raw_ensembl_id_array = _adata.var['gene_ids'].values
    _adata.var_names = _raw_ensembl_id_array
    _adata.var['gene_ids'] = _converted_gene_id_array
    adata = _adata
    adata.var_names_make_unique()  # remove duplicates in ensembl id
    print(f"Raw cell counts: {adata.n_obs} | Raw gene counts: {adata.n_vars}")

    # Filter out non-singlets from demxulet output and assign patient_id to adata
    demux = pd.read_csv(demuxlet_file_path, sep='\t') # Load demuxlet results
    demux_singlets = demux[demux['BEST'].str.startswith('SNG')] # Only keep the singlets from demux matrix
    barcode_patient_dict_raw = dict(zip(demux_singlets['BARCODE'], demux_singlets['BEST'])) # Dictionary of cell barcode to patient id
    _pid_prefixes = ['SNG-1269-1-', 'SNG-1269-2-'] # get rid of the unnecessary prefixes
    barcode_patient_dict = {
        key: value[len(_pid_prefixes[0]):] if value.startswith(_pid_prefixes[0]) else
            value[len(_pid_prefixes[1]):] if value.startswith(_pid_prefixes[1]) else
            value
        for key, value in barcode_patient_dict_raw.items()
    }
    adata = adata[adata.obs_names.isin(barcode_patient_dict.keys())] # Filter adata with mapped barcode-patient pair
    _adata = adata.copy()
    _adata.obs['patient_id'] = _adata.obs_names.map(barcode_patient_dict) # Add patient_id obs to adata
    adata = _adata

    # Filter out cells with low gene count and genes with low cell count
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_cells(adata, min_counts=20)
    sc.pp.filter_genes(adata, min_cells=20)
    sc.pp.filter_genes(adata, min_counts=10)

    # Filter out doublets with Scrublet
    scrub = scrublet.Scrublet(adata.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, 
        min_cells=3, 
        min_gene_variability_pctl=85, 
        n_prin_comps=30
    )
    adata = adata[~predicted_doublets]

    # Filter out cells with excessive mitochondria, ribosome and hemoglobin RNA
    mt_threshold = 5
    ribo_threshold = 5
    hb_threshold = 5
    _adata_qc_copy = adata.copy()
    _adata_qc_copy.var['mt'] = _adata_qc_copy.var_names.str.startswith('MT-')
    _adata_qc_copy.var['ribo'] = _adata_qc_copy.var_names.str.startswith(("RPS", "RPL"))
    _adata_qc_copy.var['hb'] = _adata_qc_copy.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(_adata_qc_copy, qc_vars=['mt', 'ribo', 'hb'], percent_top=None, log1p=False, inplace=True)
    adata = adata[(_adata_qc_copy.obs['pct_counts_mt'] < mt_threshold) & 
                (_adata_qc_copy.obs['pct_counts_ribo'] < ribo_threshold) & 
                (_adata_qc_copy.obs['pct_counts_hb'] < hb_threshold)]
    print(f"Filtered cell counts: {adata.n_obs} | Filtered gene counts: {adata.n_vars}")

    # Assign patients to phenotype group
    phenotype_dict = {
        '1123': 'normal',
        '1061': 'normal',
        '1163': 'normal',
        '1111': 'normal',
        '2284': 'abnormal',
        '1053': 'abnormal',
        '1248': 'abnormal',
        '1235': 'abnormal',
    }
    _adata = adata.copy()
    _adata.obs['phenotype'] = _adata.obs['patient_id'].apply(
        lambda x: next(phenotype for key, phenotype in phenotype_dict.items() if key in x)
    )
    adata = _adata

    # Data normalization
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # sc.pp.scale(adata)
    return adata




def analyze_cluster_differential_expression(adata, top_k_genes: int = 20):
    adata = adata.copy()
    sc.tl.rank_genes_groups(adata, 
                          groupby='leiden', 
                          n_genes=adata.n_vars, 
                          method='wilcoxon')
    
    # Extract top genes for each cluster
    cluster_top_genes_dict = {}
    for cluster in adata.obs['leiden'].unique():
        gene_indices = adata.uns['rank_genes_groups']['names'][cluster][:top_k_genes]
        gene_symbol = [adata.var['gene_ids'][idx] for idx in gene_indices]
        log2fc = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][:top_k_genes]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster][:top_k_genes]
        
        cluster_genes = list(zip(gene_symbol, log2fc, pvals_adj))
        cluster_genes.sort(key=lambda x: x[1], reverse=True)
        cluster_top_genes_dict[cluster] = cluster_genes
    return cluster_top_genes_dict



def leiden_and_umap(input_adata, pca_n_comps, leiden_res):
    adata = input_adata.copy()
    sc.pp.pca(adata, n_comps=pca_n_comps)
    sc.pp.neighbors(adata, use_rep='X_pca')
    sc.tl.leiden(adata, resolution=leiden_res, flavor='igraph', n_iterations=2)
    sc.tl.umap(adata, min_dist=0.65)
    sc.pl.umap(adata, color='leiden', legend_loc='on data', legend_fontoutline=2)
    return adata



def plot_multiple_gene_expression_umap(adata, gene_symbols: list, n_cols: int = 4):
    # Calculate number of rows needed
    n_rows = int(np.ceil(len(gene_symbols) / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10 * n_cols, 8 * n_rows), dpi=40)
    cmap = colors.LinearSegmentedColormap.from_list("custom", ['#d7d8d9', '#729dd6', '#0856bd'], N=256)
    
    # Plot each gene
    for idx, gene_symbol in enumerate(gene_symbols):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]
        
        # Get expression values
        expression_values = adata[:, adata.var['gene_ids']==gene_symbol].X.toarray().squeeze()
        
        # Create scatter plot
        ax.set_facecolor('white')
        scatter = ax.scatter(
            adata.obsm['X_umap'][:, 0],
            adata.obsm['X_umap'][:, 1],
            c=expression_values,
            cmap=cmap,
            s=1.4
        )
        
        # Customize appearance
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(1.0)
        
        plt.colorbar(scatter, ax=ax)
        ax.set_title(f'{gene_symbol}', fontsize=30)

    # Hide empty subplots
    for idx in range(len(gene_symbols), n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    plt.show()




def run_DEG_by_phenotype(input_adata, p_adj_cutoff=0.05):
    adata = input_adata.copy()
    
    # Exclude ensembl IDs
    non_ensembl_mask = ~adata.var['gene_ids'].str.match(r'^ENSG\d{11}$', na=False)
    adata_filtered = adata[:, non_ensembl_mask]
    
    # Initialize .uns dictionary to avoid view warning
    adata_filtered.uns = dict(adata_filtered.uns)
    
    # Run differential expression analysis on the filtered dataset
    sc.tl.rank_genes_groups(adata_filtered,
                         groupby='phenotype',
                         n_genes=adata_filtered.n_vars,
                         method='wilcoxon')
    
    # Run DEG for the 'abnormal' group
    deg_df = sc.get.rank_genes_groups_df(adata_filtered, group='abnormal')
    deg_df['gene_ids'] = adata_filtered.var.loc[deg_df['names'], 'gene_ids'].values
    filtered_df = deg_df[deg_df['pvals_adj'] <= p_adj_cutoff]
    adata_filtered.uns['deg_results'] = filtered_df
    return adata_filtered



def plot_DEG_results(adata, n_top_genes=10, dot_plot_title='', blacklist_genes=[]):
    # Remove blacklisted genes
    deg_df = adata.uns['deg_results']
    deg_df = deg_df[~deg_df['gene_ids'].isin(blacklist_genes)]
    adata_for_plotting = adata[:, ~adata.var['gene_ids'].isin(blacklist_genes)]

    # Plot heatmap
    sc.pl.rank_genes_groups_heatmap(
        adata_for_plotting,
        n_genes=n_top_genes,
        show_gene_labels=True,
        gene_symbols='gene_ids',
        dendrogram=False,
        swap_axes=True,
        show=True,
    )

    # Plot dotplot
    top_upregulated = deg_df.nlargest(n_top_genes, 'logfoldchanges')
    top_downregulated = deg_df.nsmallest(n_top_genes, 'logfoldchanges')
    plot_data = pd.concat([top_upregulated, top_downregulated])
    plot_data = plot_data.sort_values('logfoldchanges', ascending=False)
    
    # Set the figure size and DPI
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
    
    # Set white background
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    
    # Create custom colormap
    colors = ['#FF0000', '#800080', '#0000FF']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom", colors, N=256)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=0.1)
    
    # Create scatter plot
    sns.scatterplot(
        x='logfoldchanges',
        y='gene_ids',
        hue='pvals_adj',
        palette=cmap,
        hue_norm=norm,
        s=120,
        data=plot_data,
        legend=False,
        ax=ax
    )
    
    # Customize grid
    ax.yaxis.grid(True, linestyle='-', color='lightgray', alpha=0.5)
    ax.xaxis.grid(False)  # Remove vertical grid lines
    
    # Add center vertical line
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)
    
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.5)
    
    # Set title and labels
    ax.set_title(dot_plot_title)
    ax.set_xlabel('Log Fold Change')
    ax.set_ylabel('')
    
    # Create colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    cax = fig.add_axes([1, 0.3, 0.02, 0.4])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Adjusted p-value')
    plt.show()

