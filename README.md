# scRNA-seq Shiny App

A Shiny application for visualizing single-cell RNA-seq data stored in Zarr archives. This application utilizes a Python Flask proxy for efficient data streaming and writing, coupled with a Shiny frontend powered by WebGL for fast rendering.

---

## 🧰 Features

* 📊 Interactive UMAP visualization using [**`regl-scatterplot`**](https://github.com/flekschas/regl-scatterplot) (WebGL)
* ☁️ **Cloud-Native**: Streams data directly from AWS S3 or Google Cloud Storage via Zarr
* 🎯 Dynamic subdivision of UMAPs by **gene expression** (select multiple genes and compare side-by-side)
* ✨ MAGIC Imputation support with **animated transitions** between raw and imputed coordinates
* 💾 **Persistent State**: Save filtering masks and analysis results back to the cloud (requires credentials)
* 🐳 Dockerized for reproducible and portable deployment

---

## ⚙️ Configuration

The application looks for a `cloud_config.json` file in the root directory to locate datasets.

### Data Sources (`cloud_config.json`)
You can specify individual Zarr datasets via `urls` or entire buckets/folders to scan via `folders`.

```json
{
  "folders": [
    "[https://storage.googleapis.com/scrna-seqbrowser/](https://storage.googleapis.com/scrna-seqbrowser/)"
  ],
  "urls": [
    "[https://scrnaseq-browser.s3.us-east-2.amazonaws.com/corrected_doublet_decontX_mt_Duo_allergy.zarr](https://scrnaseq-browser.s3.us-east-2.amazonaws.com/corrected_doublet_decontX_mt_Duo_allergy.zarr)",
    "[https://scrnaseq-browser.s3.us-east-2.amazonaws.com/strict_epdsc_annotated_data_csc_full.zarr](https://scrnaseq-browser.s3.us-east-2.amazonaws.com/strict_epdsc_annotated_data_csc_full.zarr)"
  ]
}
````

-----

## 🔐 Authentication & Write Access

The app uses a Python proxy (`proxy/app.py`) to handle data operations.

  * **Read-Only Mode:** If no credentials are provided, the app works in read-only mode (viewing data works, but saving results/masks is disabled).
  * **Write Mode:** To enable writing data back to the Zarr store, you must provide `proxy/credentials.json`.

**⚠️ Security Note:** Ensure `proxy/credentials.json` is in your `.gitignore`. Never commit real keys to GitHub.

### Credentials Format (`proxy/credentials.json`)

```json
{
    "default": {
        "aws_access_key_id": "YOUR_ACCESS_KEY",
        "aws_secret_access_key": "YOUR_SECRET_KEY"
    },
    "gcp_cred": {
        "gcp_service_account_key": "./path-to-gcp-key.json"
    }
}
```

  * **default**: Used for AWS S3 buckets (and S3-compatible storage).
  * **gcp\_cred**: Used for Google Cloud Storage specific authentication.

-----

## 🚀 Quickstart (Docker)

Running via Docker is recommended as it handles both the R Shiny server and the Python Proxy automatically.

### 1\. Build the image

```bash
docker build -t scrna-shiny-app .
```

### 2\. Run the container

```bash
docker run --rm -p 3838:3838 scrna-shiny-app
```

Then open your browser at [http://localhost:3838](https://www.google.com/search?q=http://localhost:3838).

-----

## 🥪 Local Development

To run the app locally without Docker, you must run the backend (Proxy) and frontend (Shiny) separately.

### 1\. Start the Data Proxy (Python)

This handles reading/writing to the Zarr archives.

```bash
# Activate your environment first
python proxy/app.py
```

*This will typically start the API on `localhost:8080` (or similar).*

### 2\. Start the Shiny App (R)

In a separate terminal window:

```bash
R -e "shiny::runApp('.', host='0.0.0.0', port=3838)"
```

*Note: Ensure your conda environment is activated or configured in `reticulate`.*

-----

## 📦 Dependencies

### R (via conda)

  * `shiny`, `reticulate`, `dplyr`, `jsonlite`, `shinyjs`

### Python (via conda)

  * `flask`, `fsspec`, `zarr`, `pandas`, `anndata`, `scanpy`, `magic-impute`

See [`yml/shiny_app_env.yml`](https://www.google.com/search?q=yml/shiny_app_env.yml) for full environment details.

-----

## 📬 Contact

For questions, issues or suggestions, feel free to open an [Issue](https://github.com/george123ya/scRNAseq_Visualizer/issues).

-----

## 📜 License

MIT License