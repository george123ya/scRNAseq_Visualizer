FROM continuumio/miniconda3

# Set environment paths
ENV CONDA_ENV=shiny_app_env
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/envs/$CONDA_ENV/bin:$PATH
ENV CONDA_DEFAULT_ENV=$CONDA_ENV

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gdebi-core \
    wget \
    curl \
    bzip2 \
    build-essential \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install mamba to speed up env creation
RUN conda install -n base -c conda-forge mamba -y

RUN pip install gdown

# Copy and create conda environment
COPY yml/shiny_app_env.yml /tmp/
RUN mamba env create -f /tmp/shiny_app_env.yml

RUN conda info --envs
RUN ls -l /opt/conda/envs/

# Copy app and necessary data
COPY . /home/shiny-app/
# COPY processed_data/pbmc10k_subset.h5ad /home/shiny-app/processed_data/
# COPY processed_data/relaxed_epdsc_annotated_data.h5 /home/shiny-app/processed_data/

RUN mkdir -p /home/shiny-app/processed_data

WORKDIR /home/shiny-app

# Expose Shiny Server port
EXPOSE 3838

# Run the R Shiny app
# CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate shiny_app_env && Rscript app.R"]
CMD ["R", "-e", "shiny::runApp('.', host='0.0.0.0', port=3838)"]