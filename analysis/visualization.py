import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import networkx as nx
import os
import base64
from io import BytesIO


class Visualizer:
    """Generates all visualizations for the project."""

    def __init__(self, output_dir='static/images'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        sns.set_style('whitegrid')
        plt.rcParams['figure.dpi'] = 150

    def _save_plot(self, filename, dpi=150):
        """Save plot to file and return path."""
        filepath = os.path.join(self.output_dir, filename)
        plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
        plt.close()
        return filepath

    def _fig_to_base64(self):
        """Convert current figure to base64 string."""
        buf = BytesIO()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode('utf-8')
        plt.close()
        return img_str

    def plot_volcano(self, results_df, title='Volcano Plot',
                     lfc_threshold=1.0, p_threshold=0.05, save=True):
        """Generate volcano plot of differential expression."""
        fig, ax = plt.subplots(figsize=(10, 8))

        up = results_df[results_df['regulation'] == 'Upregulated']
        down = results_df[results_df['regulation'] == 'Downregulated']
        ns = results_df[results_df['regulation'] == 'Not Significant']

        ax.scatter(ns['log2_fold_change'], -np.log10(ns['p_adjusted'] + 1e-300),
                   c='gray', alpha=0.4, s=8, label=f'Not Significant ({len(ns)})')
        ax.scatter(up['log2_fold_change'], -np.log10(up['p_adjusted'] + 1e-300),
                   c='red', alpha=0.6, s=15, label=f'Upregulated ({len(up)})')
        ax.scatter(down['log2_fold_change'], -np.log10(down['p_adjusted'] + 1e-300),
                   c='blue', alpha=0.6, s=15, label=f'Downregulated ({len(down)})')

        top_genes = pd.concat([up.nsmallest(5, 'p_adjusted'),
                               down.nsmallest(5, 'p_adjusted')])
        for _, row in top_genes.iterrows():
            ax.annotate(row.name,
                        (row['log2_fold_change'],
                         -np.log10(row['p_adjusted'] + 1e-300)),
                        fontsize=7, alpha=0.8,
                        arrowprops=dict(arrowstyle='->', alpha=0.5))

        ax.axvline(-lfc_threshold, color='black', linestyle='--', alpha=0.3)
        ax.axvline(lfc_threshold, color='black', linestyle='--', alpha=0.3)
        ax.axhline(-np.log10(p_threshold), color='black', linestyle='--', alpha=0.3)

        ax.set_xlabel('Log2 Fold Change', fontsize=12, fontweight='bold')
        ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(True, alpha=0.2)

        if save:
            return self._save_plot('volcano_plot.png')
        return self._fig_to_base64()

    def plot_heatmap(self, df, n_genes=50, title='Gene Expression Heatmap',
                     save=True):
        """Generate heatmap of top differentially expressed genes."""
        if n_genes < len(df):
            plot_df = df.head(n_genes)
        else:
            plot_df = df

        fig, ax = plt.subplots(figsize=(12, max(8, n_genes * 0.2)))
        sns.heatmap(plot_df.select_dtypes(include=[np.number]),
                    cmap='RdBu_r', center=0, robust=True,
                    xticklabels=True, yticklabels=True,
                    cbar_kws={'label': 'Expression Level'},
                    ax=ax)
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        plt.yticks(fontsize=6)
        plt.xticks(rotation=45, fontsize=6)

        if save:
            return self._save_plot('heatmap.png')
        return self._fig_to_base64()

    def plot_ppi_network(self, G, hub_genes=None, title='PPI Network',
                         save=True):
        """Generate PPI network visualization."""
        fig, ax = plt.subplots(figsize=(14, 12))

        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)

        degrees = dict(G.degree())
        node_sizes = [max(degrees[n] * 30, 50) for n in G.nodes()]
        node_colors = ['red' if n in hub_genes.index.tolist()
                       else 'lightblue' for n in G.nodes()]

        nx.draw_networkx_edges(G, pos, alpha=0.15, edge_color='gray', ax=ax)
        nx.draw_networkx_nodes(G, pos, node_size=node_sizes,
                               node_color=node_colors, alpha=0.8, ax=ax)

        if hub_genes is not None:
            hub_labels = {g: g for g in hub_genes.index[:10]
                          if g in G.nodes()}
            nx.draw_networkx_labels(G, pos, labels=hub_labels,
                                    font_size=8, font_weight='bold', ax=ax)

        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axis('off')

        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='red', markersize=10,
                       label='Hub Genes'),
            plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor='lightblue', markersize=10,
                       label='Other Genes'),
        ]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10)

        if save:
            return self._save_plot('ppi_network.png')
        return self._fig_to_base64()

    def plot_centrality_barplot(self, centrality_df, n=10,
                                title='Top Hub Genes - Degree Centrality',
                                save=True):
        """Generate barplot of top hub genes."""
        top_n = centrality_df.head(n)
        fig, ax = plt.subplots(figsize=(10, 6))
        bars = ax.barh(range(len(top_n)), top_n['Degree'].values,
                       color='steelblue', edgecolor='navy')
        ax.set_yticks(range(len(top_n)))
        ax.set_yticklabels(top_n.index, fontsize=10)
        ax.set_xlabel('Degree Centrality', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.invert_yaxis()

        for bar, val in zip(bars, top_n['Degree'].values):
            ax.text(val + 0.01, bar.get_y() + bar.get_height() / 2,
                    f'{val:.3f}', va='center', fontsize=9)

        if save:
            return self._save_plot('hub_genes_barplot.png')
        return self._fig_to_base64()

    def plot_go_enrichment(self, go_df, ontology='BP', n=10,
                           title=None, save=True):
        """Generate GO enrichment dot plot."""
        if go_df.empty:
            return None
        top_n = go_df.head(n)
        if title is None:
            title = f'GO Enrichment - {ontology}'

        fig, ax = plt.subplots(figsize=(10, max(6, len(top_n) * 0.4)))
        colors = -np.log10(top_n['p_adjusted'].values + 1e-300)
        sizes = top_n['overlap'].values * 20

        scatter = ax.scatter(top_n['fold_enrichment'].values,
                             range(len(top_n)),
                             c=colors, s=sizes, cmap='viridis_r',
                             edgecolor='black', alpha=0.7)
        ax.set_yticks(range(len(top_n)))
        ax.set_yticklabels(top_n['term'].values, fontsize=9)
        ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.axvline(1, color='red', linestyle='--', alpha=0.3)

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('-Log10 Adjusted P-value', fontsize=10)

        if save:
            return self._save_plot(f'go_enrichment_{ontology}.png')
        return self._fig_to_base64()

    def plot_kegg_enrichment(self, kegg_df, n=15, save=True):
        """Generate KEGG pathway enrichment plot."""
        if kegg_df.empty:
            return None
        top_n = kegg_df.head(n)

        fig, ax = plt.subplots(figsize=(12, max(6, len(top_n) * 0.4)))
        colors = -np.log10(top_n['p_adjusted'].values + 1e-300)
        sizes = top_n['overlap'].values * 25

        scatter = ax.scatter(top_n['fold_enrichment'].values,
                             range(len(top_n)),
                             c=colors, s=sizes, cmap='plasma_r',
                             edgecolor='black', alpha=0.7)
        ax.set_yticks(range(len(top_n)))
        ax.set_yticklabels(top_n['term'].values, fontsize=9)
        ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
        ax.set_title('KEGG Pathway Enrichment', fontsize=14, fontweight='bold')
        ax.axvline(1, color='red', linestyle='--', alpha=0.3)

        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('-Log10 Adjusted P-value', fontsize=10)

        if save:
            return self._save_plot('kegg_enrichment.png')
        return self._fig_to_base64()

    def plot_expression_distribution(self, df, title='Expression Distribution',
                                     save=True):
        """Generate expression distribution plot."""
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        ax1 = axes[0]
        for col in df.columns[:10]:
            ax1.hist(df[col], bins=50, alpha=0.3, label=col[:15])
        ax1.set_xlabel('Expression Level', fontsize=11)
        ax1.set_ylabel('Frequency', fontsize=11)
        ax1.set_title('Sample Expression Distributions', fontsize=12)
        ax1.legend(fontsize=6, loc='upper right')

        ax2 = axes[1]
        ax2.boxplot([df[col].values for col in df.columns[:20]],
                    labels=[col[:10] for col in df.columns[:20]])
        ax2.set_xlabel('Samples', fontsize=11)
        ax2.set_ylabel('Expression Level', fontsize=11)
        ax2.set_title('Expression Boxplots', fontsize=12)
        plt.xticks(rotation=45, fontsize=7)

        fig.suptitle(title, fontsize=14, fontweight='bold')

        if save:
            return self._save_plot('expression_distribution.png')
        return self._fig_to_base64()

    def plot_ma_plot(self, results_df, title='MA Plot', save=True):
        """Generate MA plot."""
        fig, ax = plt.subplots(figsize=(10, 8))

        mean_expr = (results_df['case_mean'] + results_df['control_mean']) / 2
        lfc = results_df['log2_fold_change']
        colors = results_df['regulation'].map({
            'Upregulated': 'red', 'Downregulated': 'blue',
            'Not Significant': 'gray'
        })

        ax.scatter(mean_expr, lfc, c=colors, alpha=0.4, s=8)
        ax.axhline(0, color='black', linestyle='-', alpha=0.3)
        ax.set_xlabel('Mean Expression', fontsize=12, fontweight='bold')
        ax.set_ylabel('Log2 Fold Change', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.2)

        if save:
            return self._save_plot('ma_plot.png')
        return self._fig_to_base64()

    def plot_pca(self, df, title='PCA Plot', save=True):
        """Generate PCA plot."""
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        components = pca.fit_transform(df.T.values)
        pca_df = pd.DataFrame(components, index=df.columns,
                              columns=['PC1', 'PC2'])
        pca_df['variance'] = pca.explained_variance_ratio_

        fig, ax = plt.subplots(figsize=(10, 8))
        ax.scatter(pca_df['PC1'], pca_df['PC2'], c='steelblue',
                   s=50, alpha=0.7, edgecolor='navy')
        for sample in pca_df.index:
            ax.annotate(sample[:15], (pca_df.loc[sample, 'PC1'],
                                       pca_df.loc[sample, 'PC2']),
                        fontsize=7, alpha=0.7)

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)',
                      fontsize=12, fontweight='bold')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)',
                      fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.2)

        if save:
            return self._save_plot('pca_plot.png')
        return self._fig_to_base64()

    def generate_all_plots(self, de_results, ppi_network, centrality_df,
                           go_results, kegg_results, processed_df):
        """Generate all plots for the project."""
        plots = {}

        plots['volcano'] = self.plot_volcano(de_results)
        plots['heatmap'] = self.plot_heatmap(
            de_results.sort_values('p_adjusted'),
            n_genes=50
        )
        plots['expression_dist'] = self.plot_expression_distribution(processed_df)

        hub_genes = centrality_df.head(10) if centrality_df is not None else None
        plots['ppi'] = self.plot_ppi_network(ppi_network, hub_genes)
        plots['hub_barplot'] = self.plot_centrality_barplot(centrality_df)

        for ont, df in go_results.items():
            if not df.empty:
                plots[f'go_{ont}'] = self.plot_go_enrichment(df, ontology=ont)

        if kegg_results is not None and not kegg_results.empty:
            plots['kegg'] = self.plot_kegg_enrichment(kegg_results)

        plots['ma'] = self.plot_ma_plot(de_results)
        plots['pca'] = self.plot_pca(processed_df)

        return plots
