#!/bin/sh
docker run --rm -v`pwd`:/data deformation_model_builder:latest /data/build_grids.sh /data /usr/bin
sudo chown ccrook:ccrook *.tif
mv *.tif ..

