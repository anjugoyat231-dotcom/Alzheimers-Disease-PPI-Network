# Render Deployment Guide

## Step 1: Prepare Repository

Ensure your project is on GitHub.

Create a `gunicorn_config.py` in the project root:
```python
import multiprocessing

workers = multiprocessing.cpu_count() * 2 + 1
bind = "0.0.0.0:10000"
timeout = 120
```

## Step 2: Deploy on Render

1. Go to https://dashboard.render.com
2. Click "New +" → "Web Service"
3. Connect your GitHub account
4. Select `alzheimer-ppi-analysis`
5. Configure:

   | Setting | Value |
   |---------|-------|
   | Name | `alzheimer-ppi-analysis` |
   | Environment | `Python 3` |
   | Build Command | `pip install -r requirements.txt` |
   | Start Command | `gunicorn app:app` |
   | Plan | Free |

6. Under "Advanced" add environment variables:
   - `PYTHON_VERSION`: `3.11.0`
   - `FLASK_ENV`: `production`

7. Click "Create Web Service"

## Step 3: Access

Your app will be available at:
`https://alzheimer-ppi-analysis.onrender.com`

## Troubleshooting

- Build fails: Check `requirements.txt` for compatibility
- Timeout: Increase `gunicorn` timeout in start command
- Disk space: Use free plan (512MB is sufficient)
