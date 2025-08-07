# scRNA-seq Shiny App

A lightweight Shiny application for visualizing single-cell RNA-seq data stored in `.h5ad` format, using Python's `anndata` via `reticulate`.

## 🧰 Features

* Load and explore `.h5ad` files with UMAP/plotly visualizations
* Powered by `reticulate`, `anndata`, and `scanpy` (Python), and `shiny`, `plotly`, and `patchwork` (R)
* Dockerized for easy deployment

---

## 📁 Project Structure

```
.
├── app.R                      # Main Shiny app
├── processed_data/           # Processed .h5ad files
│   └── pbmc10k_subset.h5ad   # Subsampled example
├── yml/
│   └── shiny_app_env.yml     # Conda environment for app
├── Dockerfile                # For reproducible container builds
└── README.md
```

---

## 🚀 Quickstart (Docker)

### 1. Build the image

```bash
sudo docker build -t scrna-shiny-app .
```

### 2. Run the container

```bash
sudo docker run --rm -p 3838:3838 scrna-shiny-app
```

Then open your browser at [http://localhost:3838](http://localhost:3838)

---

## 🥪 Development Mode

If you'd like to run directly:

```bash
Rscript app.R
```

Make sure the conda env is activated, or manually call:

```r
reticulate::use_condaenv("shiny_app_env", required = TRUE)
```

---

## 📦 Dependencies

### R (installed via conda)

* `shiny`, `reticulate`, `plotly`, `dplyr`, `patchwork`, `cowplot`

### Python (inside conda env)

* `anndata`, `scanpy`, `rpy2`, etc.

See [`shiny_app_env.yml`](yml/shiny_app_env.yml) for details.

---

## 📬 Contact

For questions, issues or suggestions, feel free to open an [Issue](https://github.com/yourusername/scrna-shiny-app/issues).

---

## 📜 License

MIT License