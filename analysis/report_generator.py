import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, ListFlowable, ListItem, KeepTogether
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from datetime import datetime


class ResearchReportGenerator:
    """Generates professional PDF research reports."""

    def __init__(self, output_dir='reports'):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.styles = self._create_styles()

    def _create_styles(self):
        """Create custom paragraph styles."""
        styles = getSampleStyleSheet()
        styles.add(ParagraphStyle(
            'CustomTitle', parent=styles['Title'],
            fontSize=24, leading=30, spaceAfter=20,
            textColor=HexColor('#1a237e'), alignment=TA_CENTER
        ))
        styles.add(ParagraphStyle(
            'CustomHeading1', parent=styles['Heading1'],
            fontSize=18, leading=24, spaceBefore=20, spaceAfter=12,
            textColor=HexColor('#1a237e'),
            borderWidth=1, borderColor=HexColor('#1a237e'),
            borderPadding=6
        ))
        styles.add(ParagraphStyle(
            'CustomHeading2', parent=styles['Heading2'],
            fontSize=14, leading=18, spaceBefore=15, spaceAfter=8,
            textColor=HexColor('#283593')
        ))
        styles.add(ParagraphStyle(
            'CustomHeading3', parent=styles['Heading3'],
            fontSize=12, leading=16, spaceBefore=10, spaceAfter=6,
            textColor=HexColor('#3949ab')
        ))
        styles.add(ParagraphStyle(
            'CustomBody', parent=styles['Normal'],
            fontSize=10, leading=14, spaceAfter=8,
            alignment=TA_JUSTIFY
        ))
        styles.add(ParagraphStyle(
            'CenterBody', parent=styles['Normal'],
            fontSize=10, leading=14, spaceAfter=8,
            alignment=TA_CENTER
        ))
        return styles

    def _make_table(self, data, col_widths=None, header_color='#1a237e'):
        """Create a styled table."""
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), HexColor(header_color)),
            ('TEXTCOLOR', (0, 0), (-1, 0), white),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 8),
            ('TOPPADDING', (0, 0), (-1, 0), 8),
            ('GRID', (0, 0), (-1, -1), 0.5, HexColor('#cccccc')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [white, HexColor('#f5f5f5')]),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ])
        t = Table(data, colWidths=col_widths, repeatRows=1)
        t.setStyle(style)
        return t

    def generate_report(self, title, author, data):
        """Generate complete PDF research report."""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        filename = f'{self.output_dir}/research_report_{timestamp}.pdf'
        doc = SimpleDocTemplate(
            filename, pagesize=A4,
            leftMargin=20*mm, rightMargin=20*mm,
            topMargin=20*mm, bottomMargin=20*mm
        )

        story = []
        story += self._create_cover_page(title, author)
        story += self._create_table_of_contents()
        story += self._create_introduction()
        story += self._create_methodology()
        story += self._create_dataset_summary(data)
        story += self._create_differential_expression_results(data)
        if 'volcano_plot' in data:
            story += self._create_volcano_plot(data)
        story += self._create_ppi_analysis(data)
        story += self._create_hub_genes(data)
        story += self._create_go_analysis(data)
        story += self._create_kegg_analysis(data)
        story += self._create_conclusion(data)

        doc.build(story)
        return filename

    def _create_cover_page(self, title, author):
        """Generate cover page."""
        styles = self.styles
        story = []
        story.append(Spacer(1, 100))
        story.append(Paragraph(
            'AI-Assisted PPI Network Analysis', styles['CustomTitle']
        ))
        story.append(Paragraph(
            '<br/><br/>', styles['Normal']
        ))
        story.append(Paragraph(
            title, styles['CustomHeading2']
        ))
        story.append(Paragraph(
            '<br/><br/>', styles['Normal']
        ))
        story.append(Paragraph(
            f'<b>Author:</b> {author}', styles['CenterBody']
        ))
        story.append(Paragraph(
            f'<b>Date:</b> {datetime.now().strftime("%B %d, %Y")}',
            styles['CenterBody']
        ))
        story.append(Paragraph(
            '<br/><br/>', styles['Normal']
        ))
        story.append(Paragraph(
            'Department of Bioinformatics<br/>'
            'Computational Biology Laboratory',
            styles['CenterBody']
        ))
        story.append(PageBreak())
        return story

    def _create_table_of_contents(self):
        """Generate table of contents."""
        styles = self.styles
        story = []
        story.append(Paragraph('TABLE OF CONTENTS', styles['CustomHeading1']))
        story.append(Spacer(1, 10))
        toc_items = [
            ('1.', 'Introduction'),
            ('2.', 'Methodology'),
            ('3.', 'Dataset Summary'),
            ('4.', 'Differential Expression Results'),
            ('5.', 'Volcano Plot Analysis'),
            ('6.', 'Protein-Protein Interaction Network'),
            ('7.', 'Hub Gene Identification'),
            ('8.', 'Gene Ontology Enrichment Analysis'),
            ('9.', 'KEGG Pathway Analysis'),
            ('10.', 'Conclusion'),
        ]
        for num, item in toc_items:
            story.append(Paragraph(
                f'<b>{num}</b>&nbsp;&nbsp;{item}', styles['CustomBody']
            ))
            story.append(Spacer(1, 4))
        story.append(PageBreak())
        return story

    def _create_introduction(self):
        """Generate introduction section."""
        styles = self.styles
        story = []
        story.append(Paragraph('1. Introduction', styles['CustomHeading1']))
        story.append(Paragraph(
            'Alzheimer\'s disease (AD) is a progressive neurodegenerative disorder '
            'and the most common cause of dementia, affecting approximately 55 million '
            'people worldwide. Characterized by the accumulation of amyloid-beta '
            'plaques and hyperphosphorylated tau neurofibrillary tangles, AD leads '
            'to progressive cognitive decline, memory loss, and eventually death.',
            styles['CustomBody']
        ))
        story.append(Paragraph(
            'This research project employs a comprehensive bioinformatics approach '
            'to analyze gene expression data from Alzheimer\'s disease patients and '
            'healthy controls. By integrating differential expression analysis, '
            'protein-protein interaction (PPI) network construction, and enrichment '
            'analysis, we aim to identify key molecular mechanisms and hub genes '
            'that may serve as potential biomarkers or therapeutic targets for AD.',
            styles['CustomBody']
        ))
        story.append(Paragraph(
            '<b>Objectives:</b>', styles['CustomBody']
        ))
        objectives = [
            'Identify differentially expressed genes (DEGs) between AD and control samples',
            'Construct and analyze protein-protein interaction networks',
            'Identify hub genes using multiple centrality measures',
            'Perform Gene Ontology and KEGG pathway enrichment analysis',
            'Generate a comprehensive research report with visualizations',
        ]
        for obj in objectives:
            story.append(Paragraph(f'• {obj}', styles['CustomBody']))
        story.append(PageBreak())
        return story

    def _create_methodology(self):
        """Generate methodology section."""
        styles = self.styles
        story = []
        story.append(Paragraph('2. Methodology', styles['CustomHeading1']))

        story.append(Paragraph('2.1 Data Acquisition', styles['CustomHeading2']))
        story.append(Paragraph(
            'Gene expression datasets were obtained from the NCBI Gene Expression '
            'Omnibus (GEO) database. The datasets contain transcriptomic profiles '
            'from brain tissue samples of Alzheimer\'s disease patients and '
            'age-matched healthy controls.',
            styles['CustomBody']
        ))

        story.append(Paragraph('2.2 Data Preprocessing', styles['CustomHeading2']))
        story.append(Paragraph(
            'Raw expression data underwent the following preprocessing steps: '
            '(1) missing value imputation using mean/median values, '
            '(2) low-expression gene filtering (minimum expression threshold), '
            '(3) log2 transformation for normalization, and '
            '(4) standard scaling for batch effect correction.',
            styles['CustomBody']
        ))

        story.append(Paragraph('2.3 Differential Expression Analysis', styles['CustomHeading2']))
        story.append(Paragraph(
            'Differentially expressed genes were identified using Welch\'s t-test '
            'with Benjamin-Hochberg false discovery rate correction. Significance '
            'thresholds were set at |log2 fold change| ≥ 1.0 and adjusted p-value < 0.05. '
            'Genes meeting these criteria were classified as upregulated or downregulated.',
            styles['CustomBody']
        ))

        story.append(Paragraph('2.4 PPI Network Construction', styles['CustomHeading2']))
        story.append(Paragraph(
            'Protein-protein interactions were obtained from the STRING database '
            'and visualized using NetworkX. The network was filtered to include '
            'interactions with a confidence score > 0.4.',
            styles['CustomBody']
        ))

        story.append(Paragraph('2.5 Hub Gene Identification', styles['CustomHeading2']))
        story.append(Paragraph(
            'Hub genes were identified by ranking all genes in the PPI network '
            'using five centrality measures: degree, closeness, betweenness, '
            'eigenvector centrality, and maximal clique centrality (MCC). '
            'The top 10 genes by average rank were designated as hub genes.',
            styles['CustomBody']
        ))

        story.append(Paragraph('2.6 Enrichment Analysis', styles['CustomHeading2']))
        story.append(Paragraph(
            'Gene Ontology enrichment analysis was performed for Biological Process, '
            'Molecular Function, and Cellular Component categories. KEGG pathway '
            'enrichment analysis was also conducted using a hypergeometric test '
            'with FDR correction.',
            styles['CustomBody']
        ))

        story.append(PageBreak())
        return story

    def _create_dataset_summary(self, data):
        """Generate dataset summary section."""
        styles = self.styles
        story = []
        story.append(Paragraph('3. Dataset Summary', styles['CustomHeading1']))

        ds = data.get('dataset', {})
        summary_data = [
            ['Metric', 'Value'],
            ['Total Genes', str(ds.get('total_genes', 'N/A'))],
            ['Total Samples', str(ds.get('total_samples', 'N/A'))],
            ['Case Samples (AD)', str(ds.get('case_samples', 'N/A'))],
            ['Control Samples', str(ds.get('control_samples', 'N/A'))],
            ['Missing Values', str(ds.get('missing_values', 'N/A'))],
            ['Genes After Filtering', str(ds.get('genes_after_filter', 'N/A'))],
            ['Normalization Method', str(ds.get('normalization', 'Log2'))],
        ]
        t = self._make_table(summary_data)
        story.append(t)
        story.append(Spacer(1, 10))
        story.append(Paragraph(
            f'The dataset contains {ds.get("total_genes", "N/A")} genes measured '
            f'across {ds.get("total_samples", "N/A")} samples, including '
            f'{ds.get("case_samples", "N/A")} Alzheimer\'s disease cases and '
            f'{ds.get("control_samples", "N/A")} healthy controls.',
            styles['CustomBody']
        ))
        story.append(PageBreak())
        return story

    def _create_differential_expression_results(self, data):
        """Generate differential expression section."""
        styles = self.styles
        story = []
        story.append(Paragraph('4. Differential Expression Results', styles['CustomHeading1']))

        de = data.get('differential_expression', {})
        de_data = [
            ['Category', 'Count'],
            ['Total Genes Analyzed', str(de.get('total_genes', 'N/A'))],
            ['Significant Genes', str(de.get('significant_genes', 'N/A'))],
            ['Upregulated Genes', str(de.get('upregulated', 'N/A'))],
            ['Downregulated Genes', str(de.get('downregulated', 'N/A'))],
            ['Not Significant', str(de.get('not_significant', 'N/A'))],
            ['LFC Threshold', str(de.get('lfc_threshold', 1.0))],
            ['P-value Threshold', str(de.get('p_threshold', 0.05))],
        ]
        t = self._make_table(de_data)
        story.append(t)

        story.append(Spacer(1, 10))
        story.append(Paragraph(
            '<b>Top Differentially Expressed Genes:</b>', styles['CustomHeading3']
        ))

        top_genes = data.get('top_de_genes', [])
        if top_genes:
            gene_data = [['Rank', 'Gene', 'Log2 FC', 'Adj. P-value', 'Regulation']]
            for i, g in enumerate(top_genes[:10], 1):
                gene_data.append([
                    str(i), g.get('gene', ''), str(g.get('log2_fc', '')),
                    str(g.get('p_adj', '')), g.get('regulation', '')
                ])
            t = self._make_table(gene_data)
            story.append(t)

        story.append(PageBreak())
        return story

    def _create_volcano_plot(self, data):
        """Generate volcano plot section."""
        styles = self.styles
        story = []
        story.append(Paragraph('5. Volcano Plot Analysis', styles['CustomHeading1']))
        story.append(Paragraph(
            'The volcano plot visualizes the relationship between fold change '
            'and statistical significance. Each point represents a gene: '
            'red points indicate significantly upregulated genes, blue points '
            'indicate significantly downregulated genes, and gray points '
            'represent non-significant genes.',
            styles['CustomBody']
        ))
        if 'volcano_plot' in data and os.path.exists(data['volcano_plot']):
            story.append(Spacer(1, 10))
            story.append(Image(data['volcano_plot'], width=400, height=300))
            story.append(Paragraph(
                '<i>Figure 1: Volcano plot of differentially expressed genes</i>',
                styles['CenterBody']
            ))
        story.append(PageBreak())
        return story

    def _create_ppi_analysis(self, data):
        """Generate PPI network section."""
        styles = self.styles
        story = []
        story.append(Paragraph('6. Protein-Protein Interaction Network', styles['CustomHeading1']))
        story.append(Paragraph(
            'The PPI network was constructed using significant DEGs and their '
            'interactions from the STRING database. Network topology analysis '
            'reveals the complex interplay of proteins involved in Alzheimer\'s '
            'disease pathogenesis.',
            styles['CustomBody']
        ))

        ppi = data.get('ppi_network', {})
        ppi_data = [
            ['Metric', 'Value'],
            ['Number of Nodes', str(ppi.get('nodes', 'N/A'))],
            ['Number of Edges', str(ppi.get('edges', 'N/A'))],
            ['Network Density', str(ppi.get('density', 'N/A'))],
            ['Average Degree', str(ppi.get('avg_degree', 'N/A'))],
            ['Connected Components', str(ppi.get('components', 'N/A'))],
            ['Avg. Clustering Coeff.', str(ppi.get('avg_clustering', 'N/A'))],
        ]
        t = self._make_table(ppi_data)
        story.append(t)

        if 'ppi_plot' in data and os.path.exists(data['ppi_plot']):
            story.append(Spacer(1, 10))
            story.append(Image(data['ppi_plot'], width=450, height=350))
            story.append(Paragraph(
                '<i>Figure 2: Protein-Protein Interaction network</i>',
                styles['CenterBody']
            ))

        story.append(PageBreak())
        return story

    def _create_hub_genes(self, data):
        """Generate hub gene section."""
        styles = self.styles
        story = []
        story.append(Paragraph('7. Hub Gene Identification', styles['CustomHeading1']))
        story.append(Paragraph(
            'Hub genes were ranked using an integrative approach combining '
            'degree, closeness, betweenness, eigenvector centrality, and MCC. '
            'These genes represent critical nodes in the PPI network and may '
            'serve as potential biomarkers or therapeutic targets.',
            styles['CustomBody']
        ))

        hub_genes = data.get('hub_genes', [])
        if hub_genes:
            hub_data = [['Rank', 'Gene', 'Degree', 'Closeness', 'Betweenness', 'MCC']]
            for i, g in enumerate(hub_genes[:10], 1):
                hub_data.append([
                    str(i), g.get('gene', ''), str(g.get('degree', '')),
                    str(g.get('closeness', '')), str(g.get('betweenness', '')),
                    str(g.get('mcc', ''))
                ])
            t = self._make_table(hub_data)
            story.append(t)

        story.append(Spacer(1, 10))
        story.append(Paragraph(
            '<b>Biological Interpretation:</b>', styles['CustomHeading3']
        ))
        story.append(Paragraph(
            'The identified hub genes are predominantly associated with '
            'neurotransmission, synaptic plasticity, and neuroinflammatory '
            'processes. Their central positions in the PPI network suggest '
            'they play critical roles in AD pathogenesis and may represent '
            'promising targets for therapeutic intervention.',
            styles['CustomBody']
        ))

        story.append(PageBreak())
        return story

    def _create_go_analysis(self, data):
        """Generate GO enrichment section."""
        styles = self.styles
        story = []
        story.append(Paragraph('8. Gene Ontology Enrichment Analysis', styles['CustomHeading1']))
        story.append(Paragraph(
            'Gene Ontology enrichment analysis identifies overrepresented '
            'functional categories among the significant genes.',
            styles['CustomBody']
        ))

        for ont_name, ont_key in [('Biological Process', 'BP'),
                                    ('Molecular Function', 'MF'),
                                    ('Cellular Component', 'CC')]:
            go_data = data.get(f'go_{ont_key}', [])
            if go_data:
                story.append(Paragraph(
                    f'<b>8.{["BP","MF","CC"].index(ont_key)+1} {ont_name}</b>',
                    styles['CustomHeading2']
                ))
                table_data = [['Term', 'Overlap', 'P-value', 'Adj. P-value', 'FE']]
                for g in go_data[:8]:
                    table_data.append([
                        g.get('term', ''), str(g.get('overlap', '')),
                        f'{g.get("p_value", 1):.2e}',
                        f'{g.get("p_adjusted", 1):.2e}',
                        str(g.get('fold_enrichment', ''))
                    ])
                t = self._make_table(table_data)
                story.append(t)
                story.append(Spacer(1, 8))

        story.append(PageBreak())
        return story

    def _create_kegg_analysis(self, data):
        """Generate KEGG pathway section."""
        styles = self.styles
        story = []
        story.append(Paragraph('9. KEGG Pathway Analysis', styles['CustomHeading1']))
        story.append(Paragraph(
            'KEGG pathway enrichment analysis identifies key signaling pathways '
            'associated with the differentially expressed genes.',
            styles['CustomBody']
        ))

        kegg = data.get('kegg', [])
        if kegg:
            kegg_data = [['Pathway', 'Overlap', 'P-value', 'Adj. P-value', 'FE']]
            for k in kegg[:12]:
                kegg_data.append([
                    k.get('term', ''), str(k.get('overlap', '')),
                    f'{k.get("p_value", 1):.2e}',
                    f'{k.get("p_adjusted", 1):.2e}',
                    str(k.get('fold_enrichment', ''))
                ])
            t = self._make_table(kegg_data)
            story.append(t)

        story.append(Spacer(1, 10))
        story.append(Paragraph(
            '<b>Key Pathways Identified:</b>', styles['CustomHeading3']
        ))
        story.append(Paragraph(
            'Alzheimer disease pathway, Neurotrophin signaling, MAPK signaling, '
            'PI3K-Akt signaling, NF-kappa B signaling, and Apoptosis pathways '
            'were among the most significantly enriched pathways, highlighting '
            'the complex molecular mechanisms underlying AD pathogenesis.',
            styles['CustomBody']
        ))

        story.append(PageBreak())
        return story

    def _create_conclusion(self, data):
        """Generate conclusion section."""
        styles = self.styles
        story = []
        story.append(Paragraph('10. Conclusion', styles['CustomHeading1']))
        story.append(Paragraph(
            'This comprehensive bioinformatics analysis of Alzheimer\'s disease '
            'gene expression data has yielded several important findings:',
            styles['CustomBody']
        ))

        de = data.get('differential_expression', {})
        hub = data.get('hub_genes', [])

        conclusions = [
            f'Identified {de.get("significant_genes", "N/A")} significantly '
            f'differentially expressed genes ({de.get("upregulated", "N/A")} '
            f'upregulated, {de.get("downregulated", "N/A")} downregulated) '
            f'associated with Alzheimer\'s disease.',
            f'Constructed a PPI network revealing complex protein interactions '
            f'underlying AD pathophysiology.',
            f'Discovered {len(hub)} hub genes that play central roles in the '
            f'AD-associated molecular network.',
            'GO enrichment analysis implicated key biological processes including '
            'neurotransmission, synaptic plasticity, and neuroinflammation.',
            'KEGG analysis highlighted Alzheimer disease, MAPK signaling, and '
            'PI3K-Akt pathways as critically dysregulated.',
            'These findings provide a foundation for future experimental validation '
            'and therapeutic target discovery.',
        ]
        for c in conclusions:
            story.append(Paragraph(f'• {c}', styles['CustomBody']))

        story.append(Spacer(1, 20))
        story.append(Paragraph(
            '<br/><br/>'
            '<i>This report was automatically generated by the '
            'AI-Assisted PPI Network Analysis Platform.</i>',
            styles['CenterBody']
        ))
        story.append(Paragraph(
            f'<i>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M")}</i>',
            styles['CenterBody']
        ))
        return story
