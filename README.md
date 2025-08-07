# scRNA-seq Shiny App

A Shiny application for visualizing single-cell RNA-seq data stored in `.h5ad` format, using Python's `anndata` via `reticulate`. Now powered by WebGL for fast rendering and exploration.

---

## 🧰 Features

* 📊 Interactive UMAP visualization using [**`regl-scatterplot`**](https://github.com/flekschas/regl-scatterplot) (WebGL)
* 🎯 Dynamic subdivision of UMAPs by **gene expression** (select multiple genes and compare side-by-side)
* ✨ MAGIC Imputation support with **animated transitions** between raw and imputed coordinates
* 🧠 Efficient memory usage for large datasets
* 🔗 Integrates `reticulate`, `anndata`, `scanpy`, and `MAGIC` (Python), with `shiny` and `regl-scatterplot` (R/JS)
* 🐳 Dockerized for reproducible and portable deployment

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
R -e "shiny::runApp('.', port = 3838)"
```

Make sure the conda env is activated, or manually call:

```r
reticulate::use_condaenv("shiny_app_env", required = TRUE)
```

---

## 📦 Dependencies

### R (via conda)

* `shiny`, `reticulate`, `dplyr`, `jsonlite`, `shinyjs`

### Python (via conda)

* `anndata`, `scanpy`, `magic-impute`, `rpy2`, `numpy`, `pandas`

See [`shiny_app_env.yml`](yml/shiny_app_env.yml) for full environment details.

---

## 🧪 Advanced Features

* **Gene Expression Viewer** — select multiple genes and view UMAPs side-by-side
* **MAGIC Coordinate Toggle** — seamlessly animate between raw and imputed UMAP layouts
* **Memory Optimization** — lazy-loading and efficient plotting of large datasets (in developemnt)
* **Custom WebGL rendering** — using `regl-scatterplot` with support for interactivity (hover, lasso/box select)

---

## 📬 Contact

For questions, issues or suggestions, feel free to open an [Issue](https://github.com/george123ya/scRNAseq_Visualizer/issues).

---

## 📜 License

MIT License