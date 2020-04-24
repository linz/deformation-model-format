FROM ubuntu:18.04
# Docker to get current gdal for building grids (needs writing to GeoTIFF layers)
RUN apt-get update && \
    apt-get install -y python3 software-properties-common && \
    add-apt-repository -y ppa:ubuntugis/ubuntugis-unstable && \
    apt-get install -y gdal-bin python3-gdal && \
    rm -rf /var/lib/apt/lists/*
COPY src/* /usr/bin/
