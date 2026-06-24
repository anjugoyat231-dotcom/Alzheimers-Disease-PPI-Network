# Netlify Deployment Guide (Frontend)

## Step 1: Prepare Frontend

Create `netlify.toml` in project root:
```toml
[build]
  command = "echo 'No build step required'"
  publish = "."

[[redirects]]
  from = "/api/*"
  to = "https://your-render-app.onrender.com/api/:splat"
  status = 200
```

## Step 2: Deploy on Netlify

1. Go to https://app.netlify.com
2. Click "Add new site" → "Import an existing project"
3. Connect to GitHub
4. Select `alzheimer-ppi-analysis`
5. Configure:

   | Setting | Value |
   |---------|-------|
   | Build command | (leave empty) |
   | Publish directory | `.` |
   | Functions directory | (leave empty) |

6. Click "Deploy site"

## Step 3: Configure Backend Proxy

Set environment variables in Netlify:
- `API_URL`: Your Render backend URL

## Step 4: Custom Domain (Optional)

1. Go to Site settings → Domain management
2. Add custom domain
3. Update DNS records
