#!/bin/sh
mkdir src
cp ../../tools/* src
docker build --tag=deformation_model_builder:latest .
rm -rf src
