import os
import json
import pandas as pd
import numpy as np
from flask import (
    Flask, render_template, request, jsonify,
    send_file, redirect, url_for, session
)
from flask_cors import CORS
from werkzeug.utils import secure_filename
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

from analysis.preprocessing import DataPreprocessor
from analysis.differential_expression import DifferentialExpressionAnalyzer
from analysis.ppi_analysis import PPIAnalyzer
from analysis.enrichment_analysis import EnrichmentAnalyzer
from analysis.visualization import Visualizer
from analysis.report_generator import ResearchReportGenerator

app = Flask(__name__)
app.secret_key = os.urandom(24).hex()
CORS(app)

UPLOAD_FOLDER = 'uploads'
PROCESSED_FOLDER = 'data/processed'
RAW_FOLDER = 'data/raw'
REPORTS_FOLDER = 'reports'
ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt', 'gz'}

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(PROCESSED_FOLDER, exist_ok=True)
os.makedirs(RAW_FOLDER, exist_ok=True)
os.makedirs(REPORTS_FOLDER, exist_ok=True)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024


def allowed_file(filename):
    return '.' in filename and \
        filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def get_session_data():
    return {
        'dataset_loaded': session.get('dataset_loaded', False),
        'preprocessed': session.get('preprocessed', False),
        'de_analyzed': session.get('de_analyzed', False),
        'ppi_analyzed': session.get('ppi_analyzed', False),
        'enrichment_done': session.get('enrichment_done', False),
        'current_file': session.get('current_file', ''),
    }


@app.route('/')
def index():
    return render_template('index.html',
                         session_data=get_session_data())


@app.route('/dashboard')
def dashboard():
    sdata = get_session_data()
    stats = {}
    if os.path.exists('reports/analysis_state.json'):
        with open('reports/analysis_state.json') as f:
            stats = json.load(f)
    return render_template('dashboard.html',
                         session_data=sdata,
                         stats=stats)


@app.route('/upload', methods=['GET', 'POST'])
def upload():
    sdata = get_session_data()
    if request.method == 'POST':
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            session['current_file'] = filepath
            session['dataset_loaded'] = True
            session['preprocessed'] = False
            session['de_analyzed'] = False
            session['ppi_analyzed'] = False
            session['enrichment_done'] = False

            df = pd.read_csv(filepath)
            preview = df.head(10).to_html(classes='table table-striped',
                                          index=False)
            return jsonify({
                'success': True,
                'filename': filename,
                'rows': len(df),
                'cols': len(df.columns),
                'columns': list(df.columns),
                'preview': preview
            })
        return jsonify({'error': 'File type not allowed'}), 400
    return render_template('upload.html', session_data=sdata)


@app.route('/preprocess', methods=['GET', 'POST'])
def preprocess():
    sdata = get_session_data()
    if not session.get('dataset_loaded'):
        return redirect(url_for('upload'))

    if request.method == 'POST':
        data = request.get_json()
        preprocessor = DataPreprocessor()
        filepath = session.get('current_file', '')

        if not os.path.exists(filepath):
            return jsonify({'error': 'Dataset file not found'}), 400

        try:
            df = preprocessor.run_pipeline(
                filepath,
                missing_strategy=data.get('missing_strategy', 'mean'),
                min_expression=float(data.get('min_expression', 1.0)),
                min_samples=float(data.get('min_samples', 0.5)),
                normalization=data.get('normalization', 'log2'),
                scaling=data.get('scaling', 'standard')
            )

            processed_path = os.path.join(PROCESSED_FOLDER,
                                          'processed_data.csv')
            preprocessor.save_processed_data(processed_path)
            session['processed_file'] = processed_path
            session['preprocessed'] = True

            summary = preprocessor.get_summary()
            return jsonify({'success': True, 'summary': summary})
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('preprocess.html', session_data=sdata)


@app.route('/differential_expression', methods=['GET', 'POST'])
def differential_expression():
    sdata = get_session_data()
    if not session.get('preprocessed'):
        return redirect(url_for('preprocess'))

    if request.method == 'POST':
        data = request.get_json()
        analyzer = DifferentialExpressionAnalyzer()
        processed_file = session.get('processed_file', '')

        if not os.path.exists(processed_file):
            return jsonify({'error': 'Processed data not found'}), 400

        try:
            df = pd.read_csv(processed_file, index_col=0)
            results = analyzer.run_analysis(
                df,
                condition_col=data.get('condition_col', 'condition'),
                case_group=data.get('case_group', 'AD'),
                control_group=data.get('control_group', 'Control'),
                lfc_threshold=float(data.get('lfc_threshold', 1.0)),
                p_threshold=float(data.get('p_threshold', 0.05))
            )

            results_path = os.path.join(PROCESSED_FOLDER,
                                        'de_results.csv')
            analyzer.export_results(results_path)
            sig_path = os.path.join(PROCESSED_FOLDER,
                                    'significant_genes.csv')
            analyzer.export_results(sig_path, significant_only=True)
            session['de_results'] = results_path
            session['sig_genes'] = sig_path
            session['de_analyzed'] = True

            visualizer = Visualizer()
            visualizer.plot_volcano(results)

            summary = analyzer.get_summary()
            top_genes = analyzer.get_top_genes(20)
            top_list = []
            if top_genes is not None:
                for gene, row in top_genes.iterrows():
                    top_list.append({
                        'gene': gene,
                        'log2_fc': round(row.get('log2_fold_change', 0), 3),
                        'p_adj': f'{row.get("p_adjusted", 1):.2e}',
                        'regulation': row.get('regulation', '')
                    })

            return jsonify({
                'success': True,
                'summary': summary,
                'top_genes': top_list
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('differential_expression.html',
                         session_data=sdata)


@app.route('/ppi_network', methods=['GET', 'POST'])
def ppi_network():
    sdata = get_session_data()
    if not session.get('de_analyzed'):
        return redirect(url_for('differential_expression'))

    if request.method == 'POST':
        data = request.get_json() or {}
        analyzer = PPIAnalyzer()
        sig_file = session.get('sig_genes', '')

        if not os.path.exists(sig_file):
            return jsonify({'error': 'Significant genes file not found'}), 400

        try:
            sig_df = pd.read_csv(sig_file, index_col=0)
            gene_list = sig_df.index.tolist()

            interactions_file = data.get('interactions_file', '')
            interactions_df = None
            if interactions_file and os.path.exists(interactions_file):
                interactions_df = analyzer.load_string_interactions(
                    interactions_file
                )

            G = analyzer.create_network_from_genes(
                gene_list, interactions_df,
                score_threshold=float(data.get('score_threshold', 0.4))
            )

            centrality = analyzer.compute_all_centralities()
            hub_genes = analyzer.get_hub_genes(n=10)
            network_stats = analyzer.get_network_stats()

            analyzer.export_cytoscape(
                os.path.join(PROCESSED_FOLDER, 'cytoscape_network.csv')
            )
            analyzer.export_graphml(
                os.path.join(PROCESSED_FOLDER, 'network.graphml')
            )
            centrality.to_csv(
                os.path.join(PROCESSED_FOLDER, 'centrality_results.csv')
            )

            session['ppi_analyzed'] = True

            visualizer = Visualizer()
            visualizer.plot_ppi_network(G, hub_genes)
            visualizer.plot_centrality_barplot(centrality)

            hub_list = []
            for gene, row in hub_genes.iterrows():
                hub_list.append({
                    'gene': gene,
                    'degree': round(row.get('Degree', 0), 4),
                    'closeness': round(row.get('Closeness', 0), 4),
                    'betweenness': round(row.get('Betweenness', 0), 4),
                    'eigenvector': round(row.get('Eigenvector', 0), 4),
                    'mcc': round(row.get('MCC', 0), 4)
                })

            return jsonify({
                'success': True,
                'network_stats': network_stats,
                'hub_genes': hub_list,
                'total_genes_analyzed': len(gene_list)
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('ppi_network.html', session_data=sdata)


@app.route('/enrichment', methods=['GET', 'POST'])
def enrichment():
    sdata = get_session_data()
    if not session.get('ppi_analyzed'):
        return redirect(url_for('ppi_network'))

    if request.method == 'POST':
        analyzer = EnrichmentAnalyzer()
        sig_file = session.get('sig_genes', '')

        if not os.path.exists(sig_file):
            return jsonify({'error': 'Significant genes file not found'}), 400

        try:
            sig_df = pd.read_csv(sig_file, index_col=0)
            gene_list = sig_df.index.tolist()

            go_results = analyzer.run_all_go(gene_list)
            kegg_results = analyzer.run_kegg_enrichment(gene_list)
            summary = analyzer.get_summary()

            session['enrichment_done'] = True

            visualizer = Visualizer()
            for ont, df in go_results.items():
                if not df.empty:
                    visualizer.plot_go_enrichment(df, ontology=ont)
            if kegg_results is not None and not kegg_results.empty:
                visualizer.plot_kegg_enrichment(kegg_results)

            def serialize_go(df):
                if df.empty:
                    return []
                return df.head(15).to_dict('records')

            return jsonify({
                'success': True,
                'summary': summary,
                'go_bp': serialize_go(go_results.get('BP', pd.DataFrame())),
                'go_mf': serialize_go(go_results.get('MF', pd.DataFrame())),
                'go_cc': serialize_go(go_results.get('CC', pd.DataFrame())),
                'kegg': (kegg_results.head(15).to_dict('records')
                         if kegg_results is not None and not kegg_results.empty
                         else [])
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('enrichment.html', session_data=sdata)


@app.route('/report', methods=['GET', 'POST'])
def report():
    sdata = get_session_data()
    if not session.get('enrichment_done'):
        return redirect(url_for('enrichment'))

    if request.method == 'POST':
        data = request.get_json() or {}
        report_gen = ResearchReportGenerator(output_dir=REPORTS_FOLDER)

        report_data = _collect_report_data()
        report_data['volcano_plot'] = 'static/images/volcano_plot.png'
        report_data['ppi_plot'] = 'static/images/ppi_network.png'

        try:
            filename = report_gen.generate_report(
                title=data.get('title', 'Alzheimer\'s Disease PPI Network Analysis'),
                author=data.get('author', 'Bioinformatics Research Team'),
                data=report_data
            )
            return jsonify({
                'success': True,
                'report_file': filename
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 500

    return render_template('report.html', session_data=sdata)


def _collect_report_data():
    """Collect all analysis results for report generation."""
    data = {}

    if os.path.exists('reports/analysis_state.json'):
        with open('reports/analysis_state.json') as f:
            data = json.load(f)

    if session.get('de_results') and os.path.exists(session['de_results']):
        de_df = pd.read_csv(session['de_results'], index_col=0)
        if 'log2_fold_change' in de_df.columns:
            data['differential_expression'] = {
                'total_genes': len(de_df),
                'significant_genes': int(
                    de_df.get('significant', de_df['regulation'] != 'Not Significant').sum()
                ),
                'upregulated': int(
                    (de_df['regulation'] == 'Upregulated').sum()
                ),
                'downregulated': int(
                    (de_df['regulation'] == 'Downregulated').sum()
                ),
                'not_significant': int(
                    (de_df['regulation'] == 'Not Significant').sum()
                )
            }
            top_de = []
            for gene, row in de_df.sort_values('p_adjusted').head(10).iterrows():
                top_de.append({
                    'gene': gene,
                    'log2_fc': round(row.get('log2_fold_change', 0), 3),
                    'p_adj': f'{row.get("p_adjusted", 1):.2e}',
                    'regulation': row.get('regulation', '')
                })
            data['top_de_genes'] = top_de

    if session.get('processed_file') and os.path.exists(session['processed_file']):
        proc_df = pd.read_csv(session['processed_file'], index_col=0)
        data['dataset'] = {
            'total_genes': proc_df.shape[0],
            'total_samples': proc_df.shape[1],
            'normalization': 'Log2 + Standard Scaling'
        }

    if session.get('sig_genes') and os.path.exists(session['sig_genes']):
        sig_df = pd.read_csv(session['sig_genes'], index_col=0)
        data['significant_gene_count'] = len(sig_df)
        data['significant_genes_list'] = sig_df.index.tolist()

    if os.path.exists('data/processed/centrality_results.csv'):
        cent_df = pd.read_csv('data/processed/centrality_results.csv',
                              index_col=0)
        data['hub_genes'] = []
        for gene, row in cent_df.head(10).iterrows():
            data['hub_genes'].append({
                'gene': gene,
                'degree': round(row.get('Degree', 0), 4),
                'closeness': round(row.get('Closeness', 0), 4),
                'betweenness': round(row.get('Betweenness', 0), 4),
                'mcc': round(row.get('MCC', 0), 4)
            })

    go_files = [
        ('go_BP', 'GO_BP'),
        ('go_MF', 'GO_MF'),
        ('go_CC', 'GO_CC'),
    ]
    for key, ont in go_files:
        go_path = f'data/processed/enrichment_{ont}.csv'
        if os.path.exists(go_path):
            go_df = pd.read_csv(go_path)
            data[f'go_{key}'] = go_df.to_dict('records')

    kegg_path = 'data/processed/enrichment_KEGG.csv'
    if os.path.exists(kegg_path):
        kegg_df = pd.read_csv(kegg_path)
        data['kegg'] = kegg_df.to_dict('records')

    return data


@app.route('/api/save_analysis_state', methods=['POST'])
def save_analysis_state():
    """Save current analysis state to JSON file."""
    data = request.get_json()
    with open('reports/analysis_state.json', 'w') as f:
        json.dump(data, f, indent=2)
    return jsonify({'success': True})


@app.route('/api/get_analysis_state')
def get_analysis_state():
    """Get saved analysis state."""
    if os.path.exists('reports/analysis_state.json'):
        with open('reports/analysis_state.json') as f:
            return jsonify(json.load(f))
    return jsonify({})


@app.route('/api/reset_analysis', methods=['POST'])
def reset_analysis():
    """Reset all analysis state."""
    for key in ['dataset_loaded', 'preprocessed', 'de_analyzed',
                'ppi_analyzed', 'enrichment_done', 'current_file',
                'processed_file', 'de_results', 'sig_genes']:
        session.pop(key, None)
    if os.path.exists('reports/analysis_state.json'):
        os.remove('reports/analysis_state.json')
    return jsonify({'success': True})


@app.route('/api/export/<data_type>')
def export_data(data_type):
    """Export analysis results."""
    export_map = {
        'de_results': session.get('de_results', ''),
        'significant_genes': session.get('sig_genes', ''),
        'centrality': 'data/processed/centrality_results.csv',
        'cytoscape': 'data/processed/cytoscape_network.csv',
        'network': 'data/processed/network.graphml',
    }
    filepath = export_map.get(data_type, '')
    if filepath and os.path.exists(filepath):
        return send_file(filepath, as_attachment=True)
    return jsonify({'error': 'File not found'}), 404


@app.route('/api/gene_search')
def gene_search():
    """Search for specific genes in analysis results."""
    query = request.args.get('q', '').strip().upper()
    results = []
    sig_file = session.get('sig_genes', '')
    if query and sig_file and os.path.exists(sig_file):
        sig_df = pd.read_csv(sig_file, index_col=0)
        matches = [g for g in sig_df.index if query in g.upper()]
        for gene in matches[:20]:
            row = sig_df.loc[gene]
            results.append({
                'gene': gene,
                'log2_fc': round(row.get('log2_fold_change', 0), 3),
                'p_adj': f'{row.get("p_adjusted", 1):.2e}',
                'regulation': row.get('regulation', '')
            })
    return jsonify(results)


@app.route('/api/network_data')
def network_data():
    """Get PPI network data as JSON for frontend visualization."""
    import networkx as nx
    gml_path = 'data/processed/network.graphml'
    if os.path.exists(gml_path):
        G = nx.read_graphml(gml_path)
        nodes = [{'id': n, 'label': n} for n in G.nodes()]
        edges = [{'source': u, 'target': v} for u, v in G.edges()]
        return jsonify({'nodes': nodes, 'edges': edges})
    return jsonify({'nodes': [], 'edges': []})


@app.route('/api/chart_data/<chart_type>')
def chart_data(chart_type):
    """Get chart data for frontend visualizations."""
    data = {}

    if chart_type == 'de_summary' and session.get('de_results'):
        de_df = pd.read_csv(session['de_results'], index_col=0)
        data = {
            'upregulated': int((de_df['regulation'] == 'Upregulated').sum()),
            'downregulated': int((de_df['regulation'] == 'Downregulated').sum()),
            'not_significant': int(
                (de_df['regulation'] == 'Not Significant').sum()
            ),
        }

    elif chart_type == 'hub_genes':
        cent_path = 'data/processed/centrality_results.csv'
        if os.path.exists(cent_path):
            cent_df = pd.read_csv(cent_path, index_col=0).head(10)
            data = {
                'genes': list(cent_df.index),
                'degree': cent_df['Degree'].round(4).tolist(),
                'closeness': cent_df['Closeness'].round(4).tolist(),
                'betweenness': cent_df['Betweenness'].round(4).tolist(),
            }

    elif chart_type == 'go_enrichment':
        for ont in ['BP', 'MF', 'CC']:
            go_path = f'data/processed/enrichment_GO_{ont}.csv'
            if os.path.exists(go_path):
                go_df = pd.read_csv(go_path).head(10)
                data[ont] = {
                    'terms': go_df['term'].tolist(),
                    'neg_log10_p': go_df['-log10_padj'].tolist()
                    if '-log10_padj' in go_df.columns
                    else go_df['-log10_pval'].tolist(),
                    'fold_enrichment': go_df['fold_enrichment'].tolist(),
                }

    return jsonify(data)


@app.route('/api/download_report')
def download_report():
    """Download the latest generated report."""
    report_dir = REPORTS_FOLDER
    if os.path.exists(report_dir):
        reports = sorted(
            [f for f in os.listdir(report_dir) if f.endswith('.pdf')],
            reverse=True
        )
        if reports:
            return send_file(
                os.path.join(report_dir, reports[0]),
                as_attachment=True,
                download_name='alzheimer_ppi_report.pdf'
            )
    return jsonify({'error': 'No report found'}), 404


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
