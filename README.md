# AI-Assisted PPI Network Analysis for Alzheimer's Disease

A comprehensive bioinformatics platform for Alzheimer's Disease gene expression analysis, protein-protein interaction (PPI) network construction, hub gene identification, and enrichment analysis.

## Project Overview

This end-to-end research platform processes Alzheimer's Disease gene expression data to:

1. Perform differential gene expression analysis
2. Construct and analyze PPI networks
3. Identify hub genes using multiple centrality measures
4. Perform GO and KEGG enrichment analysis
5. Generate interactive visualizations
6. Produce publication-ready PDF research reports

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Backend | Python, Flask |
| Data Science | Pandas, NumPy, SciPy |
| Machine Learning | scikit-learn, statsmodels |
| Network Analysis | NetworkX |
| Visualization | Matplotlib, Seaborn, Chart.js |
| Frontend | HTML5, CSS3, JavaScript, Bootstrap 5 |
| Report Generation | ReportLab |
| Database | File-based CSV/JSON |

## Project Structure

```
alzheimer_ppi_project/
├── app.py                      # Flask application entry point
├── requirements.txt            # Python dependencies
├── README.md                   # Documentation
├── .env                        # Environment variables
│
├── analysis/                   # Analysis modules
│   ├── __init__.py
│   ├── preprocessing.py        # Data preprocessing pipeline
│   ├── differential_expression.py  # DE analysis
│   ├── ppi_analysis.py         # PPI network analysis
│   ├── enrichment_analysis.py  # GO/KEGG enrichment
│   ├── visualization.py        # Plot generation
│   └── report_generator.py     # PDF report generation
│
├── data/                       # Data storage
│   ├── raw/                    # Raw uploaded datasets
│   └── processed/              # Processed analysis results
│
├── templates/                  # HTML templates
│   ├── base.html               # Base layout
│   ├── index.html              # Landing page
│   ├── upload.html             # Data upload page
│   ├── preprocess.html         # Preprocessing page
│   ├── differential_expression.html  # DE analysis page
│   ├── ppi_network.html        # PPI network page
│   ├── enrichment.html         # Enrichment analysis page
│   ├── dashboard.html          # Research dashboard
│   └── report.html             # Report generation page
│
├── static/                     # Static assets
│   ├── css/
│   │   └── style.css           # Custom styles
│   ├── js/
│   │   └── main.js             # Frontend JavaScript
│   └── images/                 # Generated plots
│
├── uploads/                    # User-uploaded files
├── reports/                    # Generated PDF reports
└── docs/                       # Documentation
```

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Git (optional, for version control)

### Local Installation

```bash
# 1. Clone or download the repository
git clone https://github.com/yourusername/alzheimer_ppi_project.git
cd alzheimer_ppi_project

# 2. Create virtual environment (recommended)
python -m venv venv

# 3. Activate virtual environment
# Windows:
venv\Scripts\activate
# macOS/Linux:
source venv/bin/activate

# 4. Install dependencies
pip install -r requirements.txt

# 5. Run the application
python app.py

# 6. Open browser and navigate to:
http://localhost:5000
```

### Dataset Format

Upload CSV files with the following format:
- **Rows**: Genes (first column = gene symbols/IDs)
- **Columns**: Samples
- Sample names should contain group identifiers (e.g., `AD_1`, `Control_1`)

**Recommended GEO Datasets:**
- GSE5281 - AD vs Control, Entorhinal Cortex (87 samples)
- GSE36980 - AD vs Control, Hippocampus (26 samples)
- GSE33000 - AD vs Control, Prefrontal Cortex (310 samples)

## Usage Guide

### 1. Data Upload
- Navigate to the Upload page
- Select a GEO dataset CSV file
- System validates and previews the data

### 2. Data Preprocessing
- Configure missing value imputation strategy
- Set expression filtering thresholds
- Choose normalization (log2/quantile)
- Select scaling method

### 3. Differential Expression Analysis
- Specify case and control group identifiers
- Set fold change and p-value thresholds
- View volcano plot and significant gene table
- Export results to CSV

### 4. PPI Network Analysis
- Set STRING confidence score threshold
- Upload custom interaction file (optional)
- Build and analyze the PPI network
- Identify hub genes using 5 centrality measures

### 5. Enrichment Analysis
- Run GO enrichment (BP, MF, CC)
- Run KEGG pathway enrichment
- View enriched terms and pathways

### 6. Dashboard
- View comprehensive KPI cards
- Interactive charts and visualizations
- Gene search functionality
- Export analysis results

### 7. Report Generation
- Configure report settings
- Generate professional PDF report
- All sections auto-populated

## Deployment

### Render (Python Backend)

1. Push code to GitHub
2. Create a new Web Service on Render
3. Connect GitHub repository
4. Set:
   - Build Command: `pip install -r requirements.txt`
   - Start Command: `gunicorn app:app`
5. Deploy

### Netlify (Frontend)

1. Build frontend assets
2. Create `netlify.toml` in root
3. Connect GitHub repository
4. Set build command and publish directory
5. Deploy

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Home page |
| `/dashboard` | GET | Research dashboard |
| `/upload` | GET/POST | Upload dataset |
| `/preprocess` | GET/POST | Run preprocessing |
| `/differential_expression` | GET/POST | Run DE analysis |
| `/ppi_network` | GET/POST | Run PPI analysis |
| `/enrichment` | GET/POST | Run enrichment |
| `/report` | GET/POST | Generate report |
| `/api/save_analysis_state` | POST | Save state |
| `/api/get_analysis_state` | GET | Load state |
| `/api/reset_analysis` | POST | Reset analysis |
| `/api/export/<type>` | GET | Export results |
| `/api/gene_search` | GET | Search genes |
| `/api/chart_data/<type>` | GET | Chart data |
| `/api/download_report` | GET | Download PDF |

## Features

### Data Preprocessing
- Missing value handling (mean, median, zero, drop)
- Low-expression gene filtering
- Log2 transformation and quantile normalization
- Standard and min-max scaling

### Differential Expression
- Fold change and log2 fold change computation
- Welch's t-test for significance
- Benjamini-Hochberg FDR correction
- Volcano plot and heatmap visualization

### PPI Network Analysis
- STRING database integration
- Multiple centrality measures:
  - Degree Centrality
  - Closeness Centrality
  - Betweenness Centrality
  - Eigenvector Centrality
  - Maximal Clique Centrality (MCC)
- Network topology statistics
- Cytoscape-compatible export

### Enrichment Analysis
- GO Biological Process enrichment
- GO Molecular Function enrichment
- GO Cellular Component enrichment
- KEGG pathway enrichment
- Hypergeometric test with FDR correction

### Research Dashboard
- Interactive KPI cards
- Real-time chart updates
- Gene search functionality
- Multi-tab visualization viewer

### Report Generation
- Professional PDF formatting
- Auto-populated sections
- Embedded figures and tables
- Publication-ready output

## Output Files

Generated during analysis:
- `data/processed/processed_data.csv` - Normalized expression data
- `data/processed/de_results.csv` - Differential expression results
- `data/processed/significant_genes.csv` - Significant DEGs
- `data/processed/centrality_results.csv` - Centrality measures
- `data/processed/cytoscape_network.csv` - Cytoscape-compatible export
- `data/processed/network.graphml` - NetworkX GraphML file
- `static/images/*.png` - Generated visualizations
- `reports/research_report_*.pdf` - Generated PDF reports
- `reports/analysis_state.json` - Pipeline state data

## License

This project is intended for educational and research purposes.

## Citation

If using this platform for research, please cite:
```
AI-Assisted PPI Network Analysis and Hub Gene Identification 
for Alzheimer's Disease. Bioinformatics Research Platform, 2024.
```

## Contact

For questions, suggestions, or contributions:
- GitHub: [repository URL]
- Email: [contact email]
