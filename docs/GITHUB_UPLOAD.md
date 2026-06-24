# GitHub Upload Guide

## Step 1: Create GitHub Repository

1. Go to https://github.com
2. Click "+" → "New repository"
3. Name: `alzheimer-ppi-analysis`
4. Description: "AI-Assisted PPI Network Analysis for Alzheimer's Disease"
5. Set to Public or Private
6. Do NOT initialize with README (we have ours)
7. Click "Create repository"

## Step 2: Push Local Code

```bash
# Navigate to project
cd alzheimer_ppi_project

# Initialize git
git init

# Add all files
git add .

# Create .gitignore 
echo "venv/
__pycache__/
*.pyc
.DS_Store
.env
uploads/*
!uploads/.gitkeep
reports/*.pdf
static/images/*.png" > .gitignore

git add .gitignore

# Commit
git commit -m "Initial commit: AD PPI Network Analysis platform"

# Add remote
git remote add origin https://github.com/YOUR_USERNAME/alzheimer-ppi-analysis.git

# Push
git push -u origin main
```

## Step 3: Verify

Visit `https://github.com/YOUR_USERNAME/alzheimer-ppi-analysis` to verify all files are uploaded.
