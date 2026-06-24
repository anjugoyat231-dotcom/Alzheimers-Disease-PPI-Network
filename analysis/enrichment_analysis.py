import pandas as pd
import numpy as np
from scipy import stats
import json


class EnrichmentAnalyzer:
    """Gene Ontology and KEGG pathway enrichment analysis."""

    def __init__(self):
        self.go_results = {}
        self.kegg_results = None
        self.summary = {}

    # Built-in reference gene sets for Alzheimer's Disease
    GO_BP_TERMS = {
        'Neurotransmitter secretion': ['APP', 'PSEN1', 'PSEN2', 'BACE1', 'SNCA', 'SYT1', 'SNAP25', 'VAMP2'],
        'Synaptic transmission': ['APP', 'PSEN1', 'PSEN2', 'GRIN1', 'GRIN2B', 'DLG4', 'SYNGAP1', 'SHANK3'],
        'Amyloid-beta clearance': ['APP', 'BACE1', 'PSEN1', 'PSEN2', 'ADAM10', 'BIN1', 'CLU', 'APOE'],
        'Tau protein phosphorylation': ['MAPT', 'GSK3B', 'CDK5', 'DYRK1A', 'CK1', 'PP2A', 'MARK1', 'TTBK2'],
        'Oxidative stress response': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'TXNRD1', 'NFE2L2', 'HMOX1'],
        'Neuroinflammation': ['IL6', 'TNF', 'IL1B', 'TLR4', 'MYD88', 'NFKB1', 'RELA', 'CCL2'],
        'Apoptosis': ['BCL2', 'BAX', 'BAD', 'CASP3', 'CASP9', 'CYCS', 'APAF1', 'BID'],
        'Autophagy': ['ATG5', 'ATG7', 'BECN1', 'SQSTM1', 'LAMP2', 'ULK1', 'PINK1', 'PARK2'],
        'Mitochondrial dysfunction': ['MTOR', 'PGC1A', 'TFAM', 'NDUFV1', 'SDHA', 'CYTB', 'COX1', 'ATP5A1'],
        'Neurogenesis': ['DCX', 'NES', 'SOX2', 'PAX6', 'TUBB3', 'MAP2', 'NGF', 'BDNF'],
        'Synaptic plasticity': ['CREB1', 'BDNF', 'CAMK2A', 'HOMER1', 'ARC', 'EGR1', 'FOS', 'JUN'],
        'Protein misfolding': ['HSPA1A', 'HSP90', 'DNAJB1', 'BAG3', 'HSPB1', 'STUB1', 'UBQLN1', 'VCP'],
        'Cholesterol metabolism': ['APOE', 'CLU', 'ABCA1', 'APOC1', 'LDLR', 'LRP1', 'CYP46A1', 'SOAT1'],
        'Glucose metabolism': ['SLC2A1', 'HK2', 'PKM', 'LDHA', 'PDHA1', 'ACO2', 'IDH3A', 'OGDH'],
        'Calcium signaling': ['CALM1', 'CAMK2A', 'CAMK2B', 'PPP3CA', 'ITPR1', 'RYR2', 'CALB1', 'NCALD']
    }

    GO_MF_TERMS = {
        'Acetylcholine binding': ['CHRNA1', 'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6', 'CHRNA7', 'CHRNA9'],
        'Tau protein binding': ['MAPT', 'FKBP4', 'PIN1', 'PP2A', 'GSK3B', 'CDK5', 'DYRK1A', 'MARK1'],
        'Amyloid-beta binding': ['APP', 'BACE1', 'PSEN1', 'PSEN2', 'ADAM10', 'BIN1', 'CLU', 'APOE'],
        'Kinase activity': ['GSK3B', 'CDK5', 'DYRK1A', 'MAPK1', 'MAPK3', 'AKT1', 'CK1', 'PKA'],
        'Protease activity': ['BACE1', 'ADAM10', 'PSEN1', 'CASP3', 'MMP2', 'MMP9', 'IDE', 'NEP'],
        'Antioxidant activity': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'TXNRD1', 'GSTP1', 'MGST1'],
        'DNA binding': ['TP53', 'RELA', 'NFKB1', 'CREB1', 'STAT3', 'FOXO1', 'FOXO3', 'PPARG'],
        'Receptor activity': ['GRIN1', 'GRIN2B', 'TLR4', 'TNFRSF1A', 'EGFR', 'IGF1R', 'NTRK1', 'BDNF'],
        'Ion channel activity': ['GRIN1', 'GRIN2B', 'GLRA1', 'CHRNA7', 'SCN1A', 'CACNA1C', 'KCNJ3', 'HCN1'],
        'Chaperone activity': ['HSPA1A', 'HSP90', 'DNAJB1', 'HSPB1', 'BAG3', 'HSPD1', 'HSPE1', 'CLU'],
        'Acetyltransferase activity': ['EP300', 'CREBBP', 'KAT2A', 'KAT2B', 'HDAC1', 'HDAC2', 'HDAC6', 'SIRT1'],
        'Phosphatase activity': ['PP2A', 'PPP3CA', 'PTPN1', 'DUSP1', 'PTPRB', 'PTPRC', 'PP1', 'SHP2']
    }

    GO_CC_TERMS = {
        'Synapse': ['DLG4', 'SYNGAP1', 'SHANK3', 'HOMER1', 'GRIN1', 'GRIN2B', 'SNAP25', 'VAMP2'],
        'Mitochondrion': ['MTOR', 'PGC1A', 'TFAM', 'NDUFV1', 'SDHA', 'CYTB', 'COX1', 'ATP5A1'],
        'Amyloid plaque': ['APP', 'BACE1', 'PSEN1', 'PSEN2', 'ADAM10', 'BIN1', 'CLU', 'APOE'],
        'Neurofibrillary tangle': ['MAPT', 'GSK3B', 'CDK5', 'PIN1', 'DYRK1A', 'PP2A', 'FKBP4', 'MARK1'],
        'Nucleus': ['TP53', 'RELA', 'NFKB1', 'CREB1', 'STAT3', 'FOXO1', 'FOXO3', 'PPARG'],
        'Cell membrane': ['GRIN1', 'GRIN2B', 'TLR4', 'TNFRSF1A', 'EGFR', 'IGF1R', 'NTRK1', 'CHRNA7'],
        'Cytoplasm': ['AKT1', 'GSK3B', 'BCL2', 'BAX', 'CASP3', 'MAPK1', 'MAPK3', 'SOD1'],
        'Endoplasmic reticulum': ['PSEN1', 'PSEN2', 'BACE1', 'APP', 'CALM1', 'ITPR1', 'RYR2', 'HSPA1A'],
        'Lysosome': ['LAMP2', 'SQSTM1', 'BECN1', 'ATG5', 'ATG7', 'CTSD', 'CTSB', 'CTTNB'],
        'Extracellular space': ['APOE', 'CLU', 'TNF', 'IL6', 'IL1B', 'CCL2', 'APP', 'BDNF'],
        'Proteasome': ['PSMB1', 'PSMB5', 'PSMD1', 'PSMD2', 'UBQLN1', 'VCP', 'UCHL1', 'PARK2'],
        'Golgi apparatus': ['APP', 'BACE1', 'PSEN1', 'CLU', 'APOE', 'VAMP2', 'COG1', 'COG2']
    }

    KEGG_PATHWAYS = {
        'Alzheimer disease': ['APP', 'PSEN1', 'PSEN2', 'BACE1', 'MAPT', 'APOE', 'GSK3B', 'CDK5'],
        'Parkinson disease': ['SNCA', 'PARK2', 'PINK1', 'DJ1', 'LRRK2', 'UCHL1', 'HTRA2', 'VCP'],
        'Neurotrophin signaling': ['NGF', 'BDNF', 'NTRK1', 'NTRK2', 'AKT1', 'MAPK1', 'MAPK3', 'CREB1'],
        'Apoptosis': ['BCL2', 'BAX', 'BAD', 'CASP3', 'CASP9', 'CYCS', 'APAF1', 'TP53'],
        'MAPK signaling': ['MAPK1', 'MAPK3', 'MAPK14', 'MAPK8', 'AKT1', 'EGFR', 'TNF', 'IL1B'],
        'PI3K-Akt signaling': ['AKT1', 'MTOR', 'GSK3B', 'BCL2', 'BAD', 'CASP9', 'EGFR', 'IGF1R'],
        'NF-kappa B signaling': ['NFKB1', 'RELA', 'TLR4', 'MYD88', 'TNF', 'IL1B', 'BCL2', 'BAX'],
        'mTOR signaling': ['MTOR', 'AKT1', 'RICTOR', 'RPTOR', 'TSC1', 'TSC2', 'RHEB', 'S6K1'],
        'cAMP signaling': ['CREB1', 'BDNF', 'CALM1', 'CAMK2A', 'PKA', 'GRIN1', 'GRIN2B', 'MAPK1'],
        'Calcium signaling': ['CALM1', 'CAMK2A', 'CAMK2B', 'PPP3CA', 'ITPR1', 'RYR2', 'GRIN1', 'GRIN2B'],
        'Oxidative phosphorylation': ['NDUFV1', 'SDHA', 'CYTB', 'COX1', 'ATP5A1', 'UCP2', 'UCP3', 'SOD1'],
        'Autophagy': ['ATG5', 'ATG7', 'BECN1', 'SQSTM1', 'LAMP2', 'ULK1', 'PINK1', 'PARK2'],
        'Synaptic vesicle cycle': ['SNAP25', 'VAMP2', 'SYT1', 'SYN1', 'RAB3A', 'RIMS1', 'UNC13', 'CLTA'],
        'Neuroactive ligand-receptor interaction': ['GRIN1', 'GRIN2B', 'CHRNA7', 'TLR4', 'TNFRSF1A', 'NTRK1', 'NTRK2', 'GLRA1'],
        'Cholinergic synapse': ['CHRNA1', 'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6', 'CHRNA7', 'ACHE'],
        'Dopaminergic synapse': ['DRD1', 'DRD2', 'COMT', 'MAOB', 'TH', 'DDC', 'SLC6A3', 'SNAP25'],
        'Glutamatergic synapse': ['GRIN1', 'GRIN2B', 'GRIA1', 'GRIA2', 'GRM1', 'GRM5', 'SLC1A2', 'SLC1A3'],
        'Long-term potentiation': ['CREB1', 'CAMK2A', 'CAMK2B', 'MAPK1', 'MAPK3', 'PKA', 'GRIN1', 'GRIN2B'],
        'Long-term depression': ['PPP3CA', 'GRIN1', 'GRIN2B', 'MAPK1', 'MAPK3', 'GSK3B', 'PKA', 'PP2A'],
        'Huntington disease': ['HTT', 'BDNF', 'CASP3', 'CASP9', 'SOD1', 'TFAM', 'PGC1A', 'HSPA1A']
    }

    def _hypergeom_test(self, overlap, gene_set_size, total_genes, pathway_size):
        """Perform hypergeometric test for enrichment."""
        if overlap == 0:
            return 1.0
        p_val = stats.hypergeom.sf(
            overlap - 1, total_genes, pathway_size, gene_set_size
        )
        return min(p_val, 1.0)

    def _calculate_enrichment(self, gene_list, term_dict, total_genes):
        """Calculate enrichment for a set of GO terms."""
        gene_set = set(gene_list)
        results = []
        for term, genes in term_dict.items():
            pathway_genes = set(genes)
            overlap = len(gene_set & pathway_genes)
            if overlap == 0:
                continue
            p_value = self._hypergeom_test(
                overlap, len(gene_set), total_genes, len(pathway_genes)
            )
            fold_enrichment = (overlap / len(gene_set)) / (len(pathway_genes) / total_genes)
            results.append({
                'term': term,
                'overlap': overlap,
                'gene_set_size': len(gene_set),
                'pathway_size': len(pathway_genes),
                'overlap_genes': list(gene_set & pathway_genes),
                'p_value': round(p_value, 6),
                'fold_enrichment': round(fold_enrichment, 3),
                '-log10_pval': round(-np.log10(p_value + 1e-300), 3)
            })
        return sorted(results, key=lambda x: x['p_value'])

    def run_go_enrichment(self, gene_list, ontology='BP'):
        """Run GO enrichment analysis for specified ontology."""
        total_genes = len(gene_list)
        if total_genes == 0:
            return []

        if ontology == 'BP':
            results = self._calculate_enrichment(gene_list, self.GO_BP_TERMS, total_genes)
        elif ontology == 'MF':
            results = self._calculate_enrichment(gene_list, self.GO_MF_TERMS, total_genes)
        elif ontology == 'CC':
            results = self._calculate_enrichment(gene_list, self.GO_CC_TERMS, total_genes)
        else:
            results = []

        df = pd.DataFrame(results)
        if not df.empty and 'p_value' in df.columns:
            from statsmodels.stats.multitest import multipletests
            _, df['p_adjusted'], _, _ = multipletests(df['p_value'], method='fdr_bh')
            df['p_adjusted'] = df['p_adjusted'].round(6)
            df['-log10_padj'] = round(-np.log10(df['p_adjusted'] + 1e-300), 3)

        self.go_results[ontology] = df
        return df

    def run_all_go(self, gene_list):
        """Run GO enrichment for all ontologies."""
        self.run_go_enrichment(gene_list, 'BP')
        self.run_go_enrichment(gene_list, 'MF')
        self.run_go_enrichment(gene_list, 'CC')
        return self.go_results

    def run_kegg_enrichment(self, gene_list):
        """Run KEGG pathway enrichment analysis."""
        total_genes = len(gene_list)
        if total_genes == 0:
            return pd.DataFrame()

        results = self._calculate_enrichment(gene_list, self.KEGG_PATHWAYS, total_genes)
        df = pd.DataFrame(results)
        if not df.empty and 'p_value' in df.columns:
            from statsmodels.stats.multitest import multipletests
            _, df['p_adjusted'], _, _ = multipletests(df['p_value'], method='fdr_bh')
            df['p_adjusted'] = df['p_adjusted'].round(6)
            df['-log10_padj'] = round(-np.log10(df['p_adjusted'] + 1e-300), 3)

        self.kegg_results = df
        return df

    def get_summary(self):
        """Get enrichment analysis summary."""
        summary = {}
        for ont, df in self.go_results.items():
            if not df.empty:
                summary[f'go_{ont.lower()}_terms'] = int(len(df))
                summary[f'go_{ont.lower()}_significant'] = int(
                    (df['p_adjusted'] < 0.05).sum()
                )
        if self.kegg_results is not None and not self.kegg_results.empty:
            summary['kegg_pathways'] = int(len(self.kegg_results))
            summary['kegg_significant'] = int(
                (self.kegg_results['p_adjusted'] < 0.05).sum()
            )
        self.summary = summary
        return summary
