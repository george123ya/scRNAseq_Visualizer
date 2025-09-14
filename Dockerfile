# Use lightweight micromamba base image with Debian bullseye
FROM mambaorg/micromamba:1.5.8-bullseye

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    bzip2 \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    pandoc \
 && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER


# Copy environment YAML file
COPY yml/shiny_app_env.yml /tmp/env.yml

# Create conda env with micromamba and clean caches
RUN micromamba create -y -n shiny_app_env -f /tmp/env.yml && \
    micromamba clean --all --yes

# Install peakRAM and remotes inside the environment
RUN micromamba run -n shiny_app_env R -e "install.packages(c('peakRAM', 'remotes'), repos='https://cloud.r-project.org')"

# Install your GitHub package using remotes
RUN micromamba run -n shiny_app_env R -e "remotes::install_github('george123ya/reglScatterplotR', dependencies=TRUE)"

# Copy Shiny app code
COPY . /home/shiny-app
WORKDIR /home/shiny-app

# Expose Shiny app port
EXPOSE 3838

# Run the Shiny app
CMD ["sh", "-c", "micromamba run -n shiny_app_env R -e \"shiny::runApp('.', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', 3838)))\""]
