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
    libhiredis-dev \
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

# Set flags for hiredis
# ENV PKG_CFLAGS="-I/usr/include/hiredis"
# ENV PKG_LIBS="-L/usr/lib/x86_64-linux-gnu -lhiredis"

USER root

RUN ln -s /usr/lib/x86_64-linux-gnu/libhiredis.so \
       /opt/conda/envs/shiny_app_env/lib/libhiredis.so

USER $MAMBA_USER

# Install R packages that need hiredis
RUN micromamba run -n shiny_app_env R -e "install.packages(c('redux'), repos='https://cloud.r-project.org')"

# Install remotes::install_github("bwlewis/future.redis")
RUN micromamba run -n shiny_app_env R -e "remotes::install_github('bwlewis/future.redis', ref='HEAD', dependencies=TRUE, upgrade='always', force=TRUE)" && \
    micromamba run -n shiny_app_env R -e "cat('Installed future.redis commit: ', packageDescription('future.redis')\$RemoteSha, '\n')"

# Install your GitHub package using remotes
RUN micromamba run -n shiny_app_env R -e "remotes::install_github('george123ya/reglScatterplotR', ref='HEAD', dependencies=TRUE, upgrade='always', force=TRUE)" && \
    micromamba run -n shiny_app_env R -e "cat('Installed reglScatterplot commit: ', packageDescription('reglScatterplot')\$RemoteSha, '\n')"

# Copy Shiny app code
COPY . /home/shiny-app
WORKDIR /home/shiny-app

# Expose Shiny app port
EXPOSE 3838

USER root
RUN apt-get update && apt-get install -y redis-server && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

# Run the Shiny app
CMD sh -c "\
    redis-server --daemonize yes && \
    micromamba run -n shiny_app_env python /home/shiny-app/proxy/app.py & \
    micromamba run -n shiny_app_env R -e \"shiny::runApp('.', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', 3838)))\" \
"