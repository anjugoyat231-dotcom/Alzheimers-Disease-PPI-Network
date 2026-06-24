# Installation Guide

## System Requirements

- **OS**: Windows 10+, macOS 12+, Ubuntu 20.04+
- **Python**: 3.8 - 3.11
- **RAM**: 4GB minimum, 8GB recommended
- **Storage**: 500MB for project + dependencies

## Step-by-Step Installation

### 1. Install Python

Download Python 3.8+ from [python.org](https://python.org)
- Windows: Check "Add Python to PATH" during installation
- macOS/Linux: Use package manager (brew, apt)

Verify:
```bash
python --version
pip --version
```

### 2. Create Project Directory

```bash
mkdir alzheimer_ppi_project
cd alzheimer_ppi_project
```

### 3. Set Up Virtual Environment

**Windows:**
```bash
python -m venv venv
venv\Scripts\activate
```

**macOS/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
```

### 4. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

### 5. Run Application

```bash
python app.py
```

### 6. Access Web Interface

Open http://localhost:5000 in your browser.

## Troubleshooting

### Port Already in Use
```bash
# Change port in app.py or use:
python -m flask run --port=5001
```

### Missing Dependencies
```bash
pip install -r requirements.txt --force-reinstall
```

### Matplotlib Backend Issues
Set environment variable:
```bash
# Windows PowerShell
$env:MPLBACKEND="Agg"
# macOS/Linux
export MPLBACKEND="Agg"
```
