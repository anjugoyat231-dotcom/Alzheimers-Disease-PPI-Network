import networkx as nx
import pandas as pd
import numpy as np
import json


class PPIAnalyzer:
    """Protein-Protein Interaction network analysis."""

    def __init__(self):
        self.graph = nx.Graph()
        self.centrality_results = None
        self.hub_genes = None
        self.summary = {}

    def load_string_interactions(self, filepath):
        """Load STRING database interactions from CSV."""
        df = pd.read_csv(filepath)
        required_cols = ['protein1', 'protein2', 'combined_score']
        if all(c in df.columns for c in required_cols):
            interactions = df
        elif 'node1' in df.columns and 'node2' in df.columns:
            df = df.rename(columns={
                'node1': 'protein1', 'node2': 'protein2',
                'score': 'combined_score'
            })
            interactions = df
        else:
            interactions = df
            if 'combined_score' not in interactions.columns:
                score_cols = [c for c in interactions.columns
                              if 'score' in c.lower()]
                if score_cols:
                    interactions['combined_score'] = interactions[score_cols[0]]
                else:
                    interactions['combined_score'] = 1.0
        self.interactions = interactions
        return interactions

    def create_network_from_genes(self, gene_list, interactions_df=None,
                                  score_threshold=0.4):
        """Create PPI network from gene list and interactions."""
        G = nx.Graph()
        for gene in gene_list:
            G.add_node(gene)

        if interactions_df is not None:
            for _, row in interactions_df.iterrows():
                p1 = row.get('protein1', '')
                p2 = row.get('protein2', '')
                score = row.get('combined_score', 1.0)
                p1_clean = p1.split('.')[1] if '.' in str(p1) else str(p1)
                p2_clean = p2.split('.')[1] if '.' in str(p2) else str(p2)

                if p1_clean in gene_list and p2_clean in gene_list:
                    if score >= score_threshold:
                        G.add_edge(p1_clean, p2_clean, weight=score)

        self.graph = G
        self.summary['nodes'] = int(G.number_of_nodes())
        self.summary['edges'] = int(G.number_of_edges())
        self.summary['density'] = round(nx.density(G), 4)
        return G

    def compute_degree_centrality(self):
        """Compute degree centrality for all nodes."""
        dc = nx.degree_centrality(self.graph)
        return dc

    def compute_closeness_centrality(self):
        """Compute closeness centrality."""
        cc = nx.closeness_centrality(self.graph)
        return cc

    def compute_betweenness_centrality(self, normalized=True):
        """Compute betweenness centrality."""
        bc = nx.betweenness_centrality(self.graph, normalized=normalized)
        return bc

    def compute_eigenvector_centrality(self, max_iter=1000):
        """Compute eigenvector centrality."""
        try:
            ec = nx.eigenvector_centrality(self.graph, max_iter=max_iter)
        except nx.PowerIterationFailedConvergence:
            ec = nx.eigenvector_centrality_numpy(self.graph)
        return ec

    def compute_mcc(self):
        """Compute Maximal Clique Centrality (approximation using clustering)."""
        mcc_scores = {}
        for node in self.graph.nodes():
            neighbors = list(self.graph.neighbors(node))
            if len(neighbors) < 2:
                mcc_scores[node] = 0
                continue
            subgraph = self.graph.subgraph(neighbors + [node])
            try:
                cliques = list(nx.find_cliques(subgraph))
                max_size = max([len(c) for c in cliques]) if cliques else 1
                mcc_scores[node] = max_size - 1
            except nx.NetworkXError:
                mcc_scores[node] = 0
        return mcc_scores

    def compute_all_centralities(self):
        """Compute all centrality measures."""
        degree = self.compute_degree_centrality()
        closeness = self.compute_closeness_centrality()
        betweenness = self.compute_betweenness_centrality()
        eigenvector = self.compute_eigenvector_centrality()
        mcc = self.compute_mcc()

        results_df = pd.DataFrame({
            'Degree': degree,
            'Closeness': closeness,
            'Betweenness': betweenness,
            'Eigenvector': eigenvector,
            'MCC': mcc
        })
        results_df = results_df.fillna(0)
        results_df['Degree_Rank'] = results_df['Degree'].rank(ascending=False)
        results_df['Closeness_Rank'] = results_df['Closeness'].rank(ascending=False)
        results_df['Betweenness_Rank'] = results_df['Betweenness'].rank(ascending=False)
        results_df['Eigenvector_Rank'] = results_df['Eigenvector'].rank(ascending=False)
        results_df['MCC_Rank'] = results_df['MCC'].rank(ascending=False)
        results_df['Average_Rank'] = results_df[
            ['Degree_Rank', 'Closeness_Rank', 'Betweenness_Rank',
             'Eigenvector_Rank', 'MCC_Rank']
        ].mean(axis=1)
        results_df = results_df.sort_values('Average_Rank')
        self.centrality_results = results_df
        return results_df

    def get_hub_genes(self, n=10):
        """Get top N hub genes based on average rank."""
        if self.centrality_results is None:
            self.compute_all_centralities()
        hub_genes = self.centrality_results.head(n)
        self.hub_genes = hub_genes
        self.summary['hub_genes_count'] = int(len(hub_genes))
        self.summary['top_hub_genes'] = list(hub_genes.index[:10])
        return hub_genes

    def get_network_stats(self):
        """Get comprehensive network statistics."""
        if self.graph.number_of_nodes() == 0:
            return self.summary
        self.summary['nodes'] = int(self.graph.number_of_nodes())
        self.summary['edges'] = int(self.graph.number_of_edges())
        self.summary['density'] = round(nx.density(self.graph), 4)
        self.summary['avg_degree'] = round(
            sum(dict(self.graph.degree()).values()) / self.graph.number_of_nodes(), 2
        )
        components = list(nx.connected_components(self.graph))
        self.summary['components'] = len(components)
        self.summary['largest_component_size'] = len(max(components, key=len))
        if self.graph.number_of_nodes() > 0:
            self.summary['avg_clustering'] = round(
                nx.average_clustering(self.graph), 4
            )
        return self.summary

    def export_cytoscape(self, filepath):
        """Export network to Cytoscape-compatible CSV."""
        edges = list(self.graph.edges(data=True))
        edge_df = pd.DataFrame([
            {'source': u, 'target': v, **d}
            for u, v, d in edges
        ])
        edge_df.to_csv(filepath, index=False)
        return filepath

    def export_graphml(self, filepath):
        """Export network as GraphML."""
        nx.write_graphml(self.graph, filepath)
        return filepath
