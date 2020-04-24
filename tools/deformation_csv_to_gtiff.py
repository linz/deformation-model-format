#!/usr/bin/python3
###############################################################################
# $Id$
#
#  Project:
#  Purpose:  Convert deformation model .CSV files to GTIFF format
#  Author:   Chris Crook <ccrook@linz.govt.nz>
#
###############################################################################
#  Copyright (c) 2019, Land Information New Zealand (www.linz.govt.nz)
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################
# This software adapted from ntv2_to_gtiff.py by Even Rouault

synopsis = """
Creates a GeoTIFF formatted nest grid file based on a set of CSV files 
formatted as defined for the LINZ deformation format (url).

The files to include are defined in an input file which includes the 
interpolation crs defined as EPSG:code, the number of rows in the grid,
the grid dimensions (displacement:error), and for each grid the number of
rows and columns.

Assumes grids are ordered with child grids after parent grids.

This is formatted as a JSON file with structure:

  {
  "crs": "EPSG:4959",
  "type": "3d:none",
  "copyright": "Land Information New Zealand (2013): Released under Create Commons Attribution 4.0 International",
  "description": "Darfield earthquake",
  "version": "YYYYMMDD",
  "grids": [
    {
      "filename": "patch_c1_20100904/grid_c1_L1.csv",
      "extent": [168.1,176.2,-46.75,-40.375],
      "size": [55,52]
    },
    ...
    ]
    }
"""

from distutils.version import LooseVersion
from osgeo import gdal
from osgeo import osr
from cloud_optimize_gtiff import generate_optimized_file
import argparse
import copy
import datetime
import csv
import json
import math
import os
import re
import struct
import numpy as np

gdal_min_version = "3.0.0"

if LooseVersion(gdal.__version__) < LooseVersion(gdal_min_version):
    raise RuntimeError("This program requires GDAL version {0}".format(gdal_min_version))


def get_args(args=None):
    parser = argparse.ArgumentParser(description="Convert JSON csv file list to GeoTIFF nested grid.")
    parser.add_argument("source", help="Source JSON file list")
    parser.add_argument("dest", help="Destination GeoTIFF file")
    parser.add_argument("-m", "--model-dir", help="Base directory of deformation model")
    parser.add_argument("-c", "--compact-metadata", action="store_true", help="Reduce size of metadata in GeoTIFF directory")
    parser.add_argument("-n", "--no-optimize", action="store_true", help="Don't create cloud optimized files")
    parser.add_argument(
        "-s", "--split-child-grids", action="store_true", help="Split child grids that overlap multiple parents"
    )
    parser.add_argument(
        "-i", "--individual-grids", action="store_true", help="Write separate grid files for each child (for testing)"
    )
    parser.add_argument(
        "--uint16-encoding",
        dest="uint16_encoding",
        action="store_true",
        help="Use uint16 storage with linear scaling/offseting",
    )
    return parser.parse_args()


class Grid:

    SPLIT_TOL = 1.0e-6

    def __init__(self, spec):
        self.filename = spec["filename"]
        self.sourcefile = self.filename
        self.name = self.filename.replace("/", "_")
        if self.name.endswith(".csv"):
            self.name = self.name[:-4]
        extent = spec["extent"]
        size = spec["size"]
        self.xmin = float(extent[0])
        self.xmax = float(extent[2])
        self.ymin = float(extent[1])
        self.ymax = float(extent[3])
        self.nx = int(size[0])
        self.ny = int(size[1])
        self.sourcenx = self.nx
        self.sourceny = self.ny
        self.offsetx = 0
        self.offsety = 0
        self.dx = (self.xmax - self.xmin) / (self.nx - 1)
        self.dy = (self.ymax - self.ymin) / (self.ny - 1)
        self.parent = None
        self.children = []

    def addChild(self, other):
        other.parent = self
        self.children.append(other)

    def isSubgridOf(self, other):
        return self.xmin >= other.xmin and self.xmax <= other.xmax and self.ymin >= other.ymin and self.ymax <= other.ymax

    def overlaps(self, other):
        return self.xmin < other.xmax and self.xmax > other.xmin and self.ymin < other.ymax and self.ymax > other.ymin

    def splitOn(self, other):
        # Split self to cover other, assumes that self overlaps other,
        # but not self.isSubgridOf(other)
        splitx = None
        splity = None
        if other.xmin > self.xmin and other.xmin < self.xmax:
            splitx = other.xmin
        elif other.xmax > self.xmin and other.xmax < self.xmax:
            splitx = other.xmax
        if other.ymin > self.ymin and other.ymin < self.ymax:
            splity = other.ymin
        elif other.ymax > self.ymin and other.ymax < self.ymax:
            splity = other.ymax
        assert splitx is not None or splity is not None
        grids = [self]
        if splitx:
            sx = (splitx - self.xmin) / self.dx
            nsx = int(round(sx))
            test = abs(sx - nsx)
            if test > Grid.SPLIT_TOL:
                raise RuntimeError(
                    "Cannot split grid {0} over {1} - cell x not aligned ({2})".format(self.name, other.name, test)
                )
            grid1 = copy.deepcopy(self)
            grid1.xmax = splitx
            grid1.nx = nsx + 1
            grid2 = copy.deepcopy(self)
            grid2.xmin = splitx
            grid2.nx -= nsx
            grid2.offsetx += nsx
            grids = [grid1, grid2]
        if splity:
            source = grids
            grids = []
            for g in source:
                sy = (splity - g.ymin) / g.dy
                nsy = int(round(sy))
                test = abs(sy - nsy)
                if test > Grid.SPLIT_TOL:
                    raise RuntimeError(
                        "Cannot split grid {0} over {1} - cell y not aligned ({2})".format(self.name, other.name, test)
                    )
                grid1 = copy.deepcopy(g)
                grid1.ymax = splity
                grid1.ny = nsy + 1
                grid2 = copy.deepcopy(g)
                grid2.ymin = splity
                grid2.ny -= nsy
                grid2.offsety += nsy
                grids.append(grid1)
                grids.append(grid2)
        for ng, g in enumerate(grids):
            g.name = g.name + str(ng)
        return grids


def create_unoptimized_file(sourcefilename, basefilename, args, basedir="", tmpext=""):
    gridspec = json.loads(open(sourcefilename).read())
    subgrids = [Grid(g) for g in gridspec["grids"]]
    if basedir == "":
        basedir = os.path.dirname(sourcefilename)

    # Determine source CSV fields and target bands for GeoTIFF file

    displacement_type, error_type = gridspec["type"].upper().split(":")
    csv_fields = []
    bands = []
    if displacement_type in ("HORIZONTAL", "3D"):
        csv_fields.append("de")
        csv_fields.append("dn")
        bands.append("east_offset")
        bands.append("north_offset")
    if displacement_type in ("VERTICAL", "3D"):
        csv_fields.append("du")
        bands.append("vertical_offset")
    if error_type in ("HORIZONTAL", "3D"):
        csv_fields.append("eh")
        bands.append("horizontal_uncertainty")
    if error_type in ("VERTICAL", "3D"):
        csv_fields.append("ev")
        bands.append("vertical_uncertainty")
    nbands = len(bands)
    usecols = list(range(2, len(bands) + 2))
    interpolation_method = gridspec.get("interpolation_method", "bilinear")
    if interpolation_method not in ("bilinear", "geocentric_bilinear"):
        raise RuntimeError("Invalid interpolation method in {0}".format(sourcefilename))

    # Prepare the subgrids for inclusion into the data set
    # Set the name for each grid, check the source file exists ...

    modeldir = args.model_dir
    splitgrids = []
    filenames = {}
    names = set()
    while len(subgrids) > 0:
        subgrid = subgrids.pop(0)
        filename = subgrid.filename
        if modeldir is None:
            filename = os.path.join(basedir, filename)
        else:
            filename = os.path.join(modeldir, "model", filename)
        if not os.path.exists(filename):
            raise RuntimeError("Source file {0} is missing".format(filename))
        try:
            csvr = csv.reader(open(filename))
            source_fields = next(csvr)
            csvr = None
        except:
            raise RuntimeError("Error opening CSV file {0}".format(filename))
        if source_fields[:2] != ["lon", "lat"] or source_fields[2:] != csv_fields:
            raise RuntimeError("Invalid fields in CSV {0}: should be lon, lat, {1}".format(filename, ", ".join(csv_fields)))
        subgrid.sourcefile = filename

        parent = None
        split = False
        for testgrid in reversed(splitgrids):
            if subgrid.isSubgridOf(testgrid):
                parent = testgrid
                break
            elif subgrid.overlaps(testgrid):
                if args.split_child_grids:
                    splits = subgrid.splitOn(testgrid)
                    subgrids[0:0] = splits
                    split = True
                    break
                else:
                    raise RuntimeError("Invalid overlap of grids {0} and {1}".format(subgrid.name, testgrid.name))
        if split:
            continue
        if parent is None or args.individual_grids:
            filenames[subgrid] = basefilename
        if parent is not None:
            parent.addChild(subgrid)
        name = subgrid.name
        while name in names:
            name = re.sub(r"(\d+)$", lambda m: str(int(m.group(1) or "0") + 1), name)
        names.add(name)
        subgrid.name = name
        splitgrids.append(subgrid)

    subgrids = splitgrids

    if len(filenames) > 1:
        template = re.sub(r"((?:\.\w+)?)$", r"{0}\1", basefilename)
        for nfile, parent in enumerate(filenames):
            filenames[parent] = template.format(nfile + 1)
    opened = []

    # Compile the grids

    srcfile = None
    srcdata = None

    for idx_ifd, subgrid in enumerate(subgrids):

        # Read the data from the CSV file

        csvfile = subgrid.sourcefile
        if csvfile != srcfile:
            srcfile = csvfile
            srcdata = np.genfromtxt(csvfile, names=True, delimiter=",", usecols=usecols)
            # Shape is (nrows=nlat values,ncols = n lon values)
            srcdata = srcdata.reshape((subgrid.sourceny, subgrid.sourcenx)).T
        data = srcdata[subgrid.offsetx : subgrid.offsetx + subgrid.nx, subgrid.offsety : subgrid.offsety + subgrid.ny]

        tmp_ds = gdal.GetDriverByName("GTiff").Create(
            "/vsimem/tmp", subgrid.nx, subgrid.ny, nbands, gdal.GDT_Float32 if not args.uint16_encoding else gdal.GDT_UInt16
        )

        if idx_ifd == 0 or args.individual_grids:
            description = gridspec["description"]
            copyright = gridspec["copyright"]
            version = gridspec["version"]
            griddate = "{0}:{1}:{2} 00:00:00".format(version[:4], version[4:6], version[6:8])
            tmp_ds.SetMetadataItem("TIFFTAG_IMAGEDESCRIPTION", description)
            tmp_ds.SetMetadataItem("TIFFTAG_COPYRIGHT", copyright)
            tmp_ds.SetMetadataItem("TIFFTAG_DATETIME", griddate)
            tmp_ds.SetMetadataItem("recommended_interpolation_method", interpolation_method)

        if idx_ifd == 0 or not args.compact_metadata or args.individual_grids:
            tmp_ds.SetMetadataItem("TYPE", "DEFORMATION_MODEL")
            tmp_ds.SetMetadataItem("DISPLACEMENT_TYPE", displacement_type)
            tmp_ds.SetMetadataItem("UNCERTAINTY_TYPE", error_type)
            for i, band in enumerate(bands):
                tmp_ds.GetRasterBand(i + 1).SetDescription(band)
                tmp_ds.GetRasterBand(i + 1).SetUnitType("metre")

        # Calculate the GeoTransform
        transform = [subgrid.xmin - subgrid.dx / 2.0, subgrid.dx, 0, subgrid.ymin - subgrid.dy / 2.0, 0, subgrid.dy]

        src_crs = osr.SpatialReference()
        src_crs.SetFromUserInput(gridspec["crs"])
        tmp_ds.SetSpatialRef(src_crs)
        tmp_ds.SetGeoTransform(transform)
        tmp_ds.SetMetadataItem("AREA_OR_POINT", "Point")

        tmp_ds.SetMetadataItem("grid_name", subgrid.name)
        if subgrid.parent is not None:
            tmp_ds.SetMetadataItem("parent_grid_name", subgrid.parent.name)
        if len(subgrid.children) > 0:
            tmp_ds.SetMetadataItem("number_of_nested_grids", str(len(subgrid.children)))

        for i, field in enumerate(csv_fields):
            bdata = data[field].copy()
            assert bdata.shape == (subgrid.nx, subgrid.ny)
            # Might be useful to reconsider how to handle scaling if we
            # decide to employ it. eg to try 1.0e-6, 1.0e-5 ... to better
            # round to actual values.  Worst case (?) 10m earth movement
            # would be rounded to 1mm.  Also check rounding of source data,
            # eg if already rounded to 0.0001.
            if args.uint16_encoding:
                bmin = bdata.min()
                bmax = bdata.max()
                scale = (bmax - bmin) / 65535
                bdata = (bdata - bmin) / scale
                tmp_ds.GetRasterBand(i + 1).SetOffset(min)
                tmp_ds.GetRasterBand(i + 1).SetScale(scale)
            tmp_ds.GetRasterBand(i + 1).WriteArray(bdata.T)

        options = [
            "PHOTOMETRIC=MINISBLACK",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3" if not args.uint16_encoding else "PREDICTOR=2",
            "INTERLEAVE=BAND",
            "GEOTIFF_VERSION=1.1",
        ]
        if tmp_ds.RasterXSize > 256 and tmp_ds.RasterYSize > 256:
            options.append("TILED=YES")
        else:
            options.append("BLOCKYSIZE=" + str(tmp_ds.RasterYSize))

        rootgrid = subgrid
        while rootgrid not in filenames and rootgrid.parent is not None:
            rootgrid = rootgrid.parent
        filename = filenames.get(rootgrid)
        tmpfilename = os.path.join(basedir, filename + tmpext)
        if filename not in opened:
            opened.append(filename)
            gdal.Unlink(tmpfilename)
        else:
            options.append("APPEND_SUBDATASET=YES")

        assert gdal.GetDriverByName("GTiff").CreateCopy(tmpfilename, tmp_ds, options=options)

    return opened


def build_deformation_gtiff(source, target, args, basedir=""):
    tmpext = "" if args.no_optimize else ".tmp"
    gridfiles = create_unoptimized_file(source, target, args, tmpext=tmpext, basedir=basedir)
    for gridfile in gridfiles:
        gridpath = os.path.join(basedir, gridfile)
        tmpfile = gridpath + tmpext
        generate_optimized_file(tmpfile, gridpath)
        gdal.Unlink(tmpfile)
    return gridfiles


if __name__ == "__main__":
    args = get_args()
    gridfiles = build_deformation_gtiff(args.source, args.dest, args)
    for gridfile in gridfiles:
        print("Built", gridfile)
