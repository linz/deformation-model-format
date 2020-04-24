#!/usr/bin/python3

from osgeo import gdal
import numpy as np
import math
from abc import ABC, abstractmethod


DisplacementFields = {
    "NONE": [],
    "HORIZONTAL": ["east_offset", "north_offset"],
    "VERTICAL": ["vertical_offset"],
    "3D": ["east_offset", "north_offset", "vertical_offset"],
}

UncertaintyFields = {
    "NONE": [],
    "HORIZONTAL": ["horizontal_uncertainty"],
    "VERTICAL": ["vertical_uncertainty"],
    "3D": ["horizontal_uncertainty", "vertical_uncertainty"],
}

DegreeFields = ["east_offset", "north_offset"]
BilinearMethod = "bilinear"
BilinearGeocentricMethod = "geocentric_bilinear"
BilinearGeocentricTypes = ["HORIZONTAL", "3D"]


class RangeError(ValueError):
    pass


class GridError(ValueError):
    pass


class DeformationGrid(ABC):
    class Extent:
        def __init__(self, minlon, minlat, maxlon, maxlat, geographic=True):
            self.minlon = minlon
            self.minlat = minlat
            self.maxlon = maxlon
            self.maxlat = maxlat
            self.geographic = geographic

        def wrapLongitude(self, lon):
            if self.geographic:
                while lon < self.minlon:
                    lon += 360.0
                while lon > self.maxlon:
                    lon -= 360.0
            return lon

        def containsPoint(self, lon, lat):
            lon = self.wrapLongitude(lon)
            return lon >= self.minlon and lon <= self.maxlon and lat >= self.minlat and lat <= self.maxlat

        def containsExtent(self, other):
            return (
                self.minlon <= other.minlon
                and self.maxlon >= other.maxlon
                and self.minlat <= other.minlat
                and self.maxlat >= other.maxlat
            )

        def overlaps(self, other):
            return (
                self.minlon < other.maxlon
                and self.maxlon > other.minlon
                and self.minlat < other.maxlat
                and self.maxlat > other.minlat
            )

    class Size:
        def __init__(self, ncol, nrow):
            self.ncol = ncol
            self.nrow = nrow

    class SubGrid(ABC):
        def __init__(self, id, size, extent, method, displacement_type, uncertainty_type, degrees_units, geographic_grid=True):
            self.id = id
            self.size = size
            assert size.ncol > 1 and size.nrow > 1, "Number of rows and columns must be greater than 1"
            self.extent = extent
            self.lonsize = (extent.maxlon - extent.minlon) / (size.ncol - 1)
            self.latsize = (extent.maxlat - extent.minlat) / (size.nrow - 1)
            self.children = []
            self.parent = []
            self.method = None
            if method == BilinearMethod:
                self.calcValue = self.calcValueBilinear
            elif method == BilinearGeocentricMethod:
                self.calcValue = self.calcValueGeocentricBilinear
                if degrees_units:
                    raise RuntimeError(
                        "Cannot use {0} interpolation method with degrees units".format(BilinearGeocentricMethod)
                    )
            else:
                raise ValueError(
                    'Interpolation method must be "{0}" or "{1}"'.format(BilinearMethod, BilinearGeocentricMethod)
                )
            if degrees_units and not geographic_grid:
                raise ValueError("Invalid grid - degrees units incompatible with projection CRS")
            self.displacement_type = displacement_type
            self.uncertainty_type = uncertainty_type
            self.geographic_grid = geographic_grid
            self.degrees_units = degrees_units
            self.data = None
            self.enxyz = None

        def addChildSubgrid(self, grid):
            if not self.extent.containsExtent(grid.extent):
                raise GridError("Grid {0} is not contained in grid {1}".format(grid.id, self.id))
            for c in self.children:
                if c.extent.containsExtent(grid.extent):
                    c.addChildSubgrid(grid)
                    return
                if c.extent.overlaps(grid.extent):
                    raise GridError("Grid {0} overlaps grid {1}".format(grid.id, self.id))

        def interpolationFactors(self, lon, lat):
            lon = self.extent.wrapLongitude(lon)
            fcol = (lon - self.extent.minlon) / self.lonsize
            frow = (lat - self.extent.minlat) / self.latsize
            col = max(0, min(self.size.ncol - 2, int(fcol)))
            row = max(0, min(self.size.nrow - 2, int(frow)))
            fcol -= col
            frow -= row
            indices = ((col, row), (col, row + 1), (col + 1, row), (col + 1, row + 1))
            factors = np.array([(1.0 - fcol) * (1.0 - frow), (1.0 - fcol) * frow, fcol * (1.0 - frow), fcol * frow])
            return indices, factors

        def getLonLat(self, col, row):
            return (self.extent.minlon + col * self.lonsize, self.extent.minlat + row * self.latsize)

        def calcValueBilinear(self, lon, lat):
            indices, factors = self.interpolationFactors(lon, lat)
            nodevalues = np.array([self.getValue(col, row) for col, row in indices])
            value = np.sum(nodevalues * factors.reshape(4, 1), axis=0)
            return value

        def calcValueGeocentricBilinear(self, lon, lat):
            enxyz = self.enxyz
            if enxyz is None:
                latvals = np.linspace(self.extent.minlat, self.extent.maxlat, self.size.ncol)
                latvals = np.radians(latvals)
                coslat = np.cos(latvals)
                sinlat = np.sin(latvals)
                londiff = np.radians(self.lonsize / 2.0)
                coslon = np.cos(londiff)
                sinlon = -np.sin(londiff)
                # Vectors in terms of X axis through mid longitude of cell col,row
                # evec1, nvec[row] are vectors at col,row
                # evec2, nvec[row] are vectors at col+1,row
                evec1 = np.array([-sinlon, coslon, 0.0])
                evec2 = np.array([sinlon, coslon, 0.0])
                nvec = np.transpose(np.array([-coslon * sinlat, -sinlon * sinlat, coslat]))
                self.enxyz = (evec1, evec2, nvec)
            else:
                evec1, evec2, nvec = enxyz
            indices, factors = self.interpolationFactors(lon, lat)
            nodevalues = np.array([self.getValue(col, row) for col, row in indices])
            # NOTE: This order must match row,col valuesin in interpolationFactors
            col, row = indices[0]
            envecs = np.array([[evec1, nvec[row]], [evec1, nvec[row + 1]], [evec2, nvec[row]], [evec2, nvec[row + 1]]])
            envecs[2, 1, 1] = -nvec[row, 1]
            envecs[3, 1, 1] = -nvec[row + 1, 1]
            # envecs is now east and north vectors at each node, array (4,2,3)
            # Calculate east and north vectors at calculation point relative to
            # meridian through centre of cell.
            # bottom left lon lat
            clon, clat = self.getLonLat(col, row)
            rlon = math.radians(lon - clon - self.lonsize / 2.0)
            rlat = math.radians(lat)
            coslon = math.cos(rlon)
            sinlon = math.sin(rlon)
            coslat = math.cos(rlat)
            sinlat = math.sin(rlat)
            # en vectors at calculation point, array 2,3
            envec = np.array([[-sinlon, coslon, 0.0], [-coslon * sinlat, -sinlon * sinlat, coslat]])
            # Then multiply envec at each node by corresponding displacement value and
            # sum on 0/1 axes to compute xyz displacement at calculation point
            dispvec = np.sum(nodevalues[:, :2].reshape(4, 2, 1) * envecs * factors.reshape((4, 1, 1)), axis=(0, 1))
            # dot product to calculate east/north vector at calc point
            result = dispvec.reshape((1, 3)).dot(envec.T).flatten()
            if len(nodevalues) > 2:
                hvalue = np.sum(nodevalues[:, 2:] * factors.reshape(4, 1), axis=0)
                result = np.hstack([result, hvalue])
            #!!!NOTE Need to handle degrees_units here - not allowed
            return result

        @abstractmethod
        def getValue(self, col, row):
            pass

    def __init__(self, gridfile, geographic=True):
        self.gridfile = gridfile
        self.loaded = False
        self.basegrid = None
        self.geographic = geographic

    @abstractmethod
    def _loadGrid(self):
        pass

    def load(self):
        if not self.loaded:
            self._loadGrid()
            self.loaded = True

    def addSubgrid(self, subgrid):
        if subgrid.displacement_type not in BilinearGeocentricTypes and subgrid.method == BilinearGeocentricMethod:
            subgrid.method = BilinearMethod
        if subgrid.method == BilinearGeocentricMethod and subgrid.degrees_units:
            raise RuntimeError("Cannot use {0} interpolation method with degrees units".format(BilinearGeocentricMethod))

        if self.basegrid is None:
            self.basegrid = subgrid
        else:
            if subgrid.method != self.basegrid.method:
                raise RuntimeError(
                    "Interpolation method {1} not consistent with {2}".format(subgrid.method, self.basegrid.method)
                )
            if subgrid.displacement_type != self.basegrid.displacement_type:
                raise RuntimeError(
                    "Displacement type {1} not consistent with {2}".format(
                        subgrid.displacement_type, self.basegrid.displacement_type
                    )
                )
            if subgrid.uncertainty_type != self.basegrid.uncertainty_type:
                raise RuntimeError(
                    "Uncertainty type {1} not consistent with {2}".format(
                        subgrid.uncertainty_type, self.basegrid.uncertainty_type
                    )
                )
            self.basegrid.addChildSubgrid(subgrid)

    def selectSubgrid(self, lon, lat):
        if not self.loaded:
            self.load()
        if not self.basegrid.extent.containsPoint(lon, lat):
            return None
        grid = self.basegrid
        while True:
            subgrid = None
            for child in grid.children:
                if child.extent.containsPoint(lon, lat):
                    subgrid = child
                    break
            if subgrid is None:
                break
            grid = subgrid
        return grid

    def calcValue(self, lon, lat):
        subgrid = self.selectSubgrid(lon, lat)
        if subgrid is None:
            raise RangeError()
        return subgrid.calcValue(lon, lat)

    def interpolation_method(self):
        self.load()
        return self.basegrid.method

    def displacement_type(self):
        self.load()
        return self.basegrid.displacement_type

    def uncertainty_type(self):
        self.load()
        return self.basegrid.uncertainty_type


class DeformationGridGeoTIFF(DeformationGrid):
    class SubGrid(DeformationGrid.SubGrid):
        def __init__(self, id, dataset, scaling, method, disptype, unctype, degrees_units, geographic=True):
            size = DeformationGrid.Size(dataset.RasterXSize, dataset.RasterYSize)
            (lon0, dlondcol, dlondrow, lat0, dlatdcol, dlatdrow) = dataset.GetGeoTransform()
            self.transpose = False
            if dlondcol == 0.0 and dlatdrow == 0:
                self.transpose = True
                dlondcol, dlondrow = dlondrow, dlondcol
                dlatdcol, dlatdrow = dlatdrow, dlatdcol
            assert dlondrow == 0.0 and dlatdcol == 0.0, "Can only handle grids aligned E-W and N-S"
            # Offset to centre of raster cell
            lon0 += dlondcol * 0.5
            lat0 += dlatdrow * 0.5
            lon1 = lon0 + dlondcol * (size.ncol - 1)
            lat1 = lat0 + dlatdrow * (size.nrow - 1)
            extent = DeformationGrid.Extent(min(lon0, lon1), min(lat0, lat1), max(lon0, lon1), max(lat0, lat1), geographic)
            self.reverseCols = lon1 < lon0
            self.reverseRows = lat1 < lat0
            self.scaling = scaling
            self.data = None
            self.gridname = dataset.GetMetadataItem("grid_name")
            self.parentgrid = dataset.GetMetadataItem("arent_grid_name")
            self.nchildren = int(dataset.GetMetadataItem("number_of_nested_grids") or "0")
            DeformationGrid.SubGrid.__init__(self, id, size, extent, method, disptype, unctype, degrees_units, geographic)

        def _loadData(self):
            ds = gdal.Open(self.id)
            arrays = []
            for i in range(ds.RasterCount):
                band = ds.GetRasterBand(i + 1)
                arrays.append(band.ReadAsArray())
            data = np.array(arrays)
            if self.transpose:
                data = np.transpose(data, (1, 2, 0))
            else:
                data = np.transpose(data, (2, 1, 0))
            if self.reverseCols:
                data = np.flip(data, 0)
            if self.reverseRows:
                data = np.flip(data, 1)
            self.data = data

        def getValue(self, col, row):
            if self.data is None:
                self._loadData()
            value = self.data[col, row]
            if self.scaling is not None:
                value = value * self.scaling[:, 0] + self.scaling[:, 1]
            return value

    def raiseError(self, error):
        raise GridError(error)

    def _getMetadata(self, ds, item, default=None):
        value = ds.GetMetadataItem(item)
        if value == None:
            value = default
        return value

    def _loadGrid(self):
        gridfile = self.gridfile
        gtiff = gdal.Open(gridfile)
        assert gtiff.GetDriver().ShortName == "GTiff", "Grid file {0} is not a GeoTIFF file".format(gridfile)
        subdatasets = gtiff.GetSubDatasets()
        if len(subdatasets) == 0:
            subdatasets = [(self.gridfile, None)]

        # Load the subgrid metadata
        gridtype = ""
        disptype = ""
        unctype = ""
        degrees_units = None
        interpolation_method = BilinearMethod
        first = True
        for subds in subdatasets:
            ds = gdal.Open(subds[0])
            subgridtype = self._getMetadata(ds, "TYPE", gridtype)
            if subgridtype != "DEFORMATION_MODEL":
                self.raiseError("{0} TYPE is not DEFORMATION_MODEL".format(subds[0]))
                continue
            subdisptype = self._getMetadata(ds, "DISPLACEMENT_TYPE", disptype)
            subunctype = self._getMetadata(ds, "UNCERTAINTY_TYPE", unctype)
            subinterp = self._getMetadata(ds, "recommended_interpolation_method", interpolation_method)
            if subdisptype not in DisplacementFields:
                self.raiseError("{0} DISPLACEMENT_TYPE {1} not valid".format(subds[0], disptype))
            if subunctype not in UncertaintyFields:
                self.raiseError("{0} UNCERTAINTY_TYPE {1} not valid".format(subds[0], unctype))
            if first:
                first = False
                disptype = subdisptype
                unctype = subunctype
                interpolation_method = subinterp
                fields = DisplacementFields.get(disptype, []) + UncertaintyFields.get(unctype, [])

            if subdisptype not in BilinearGeocentricTypes and subinterp == BilinearGeocentricMethod:
                subinterp = BilinearMethod

            nband = ds.RasterCount
            if nband != len(fields):
                raise self.raiseError(
                    "{0} has incorrect number of fields ({1} instead of {2})".format(subds[0], nband, len(fields))
                )

            scaling = []
            usescale = False
            for i in range(nband):
                band = ds.GetRasterBand(i + 1)
                banddesc = band.GetDescription()
                bandunit = band.GetUnitType()
                scale = band.GetScale()
                if scale is None:
                    scale = 1.0
                offset = band.GetOffset()
                if offset is None:
                    offset = 0.0
                usescale = usescale or scale != 1.0 or offset != 0.0
                scaling.append([scale, offset])
                if banddesc != fields[i]:
                    self.raiseError("{0} band {1} description {2} should be {3}".format(subds[0], i + 1, banddesc, fields[i]))
                unit = "metre"
                if banddesc in DegreeFields:
                    if degrees_units is None:
                        degrees_units = bandunit == "degree"
                    if degrees_units:
                        unit = "degree"
                if bandunit != unit:
                    self.raiseError("{0} band {1} unit {2} should be {3}".format(subds[0], i + 1, bandunit, unit))
            if not usescale:
                scaling = None
            try:

                subgrid = self.SubGrid(
                    subds[0], ds, scaling, subinterp, subdisptype, subunctype, degrees_units, self.geographic
                )
                self.addSubgrid(subgrid)
            except GridError as ex:
                self.raiseError("Error loading grid {0}: {1}".format(subds[0], str(ex)))
            except Exception as ex:
                self.raiseError("Error loading grid {0}: {1}".format(subds[0], str(ex)))

    @staticmethod
    def main(argv):
        import argparse

        parser = argparse.ArgumentParser(description="Calculate deformation in a deformation model GeoTIFF file")
        parser.add_argument("geotiff_file", help="Name of the deformation model GeoTIFF file")
        parser.add_argument("--point", nargs=2, help="lon,lat to evaluate")
        parser.add_argument("--projection", action="store_true", help="Grid is based on projection CRS")
        args = parser.parse_args(argv)
        grid = DeformationGridGeoTIFF(args.geotiff_file)
        grid.load()
        if args.point:
            (lon, lat) = (float(x) for x in args.point)
            print(grid.calcValue(lon, lat))
        # print(grid.basegrid.getValue(0, 0))


if __name__ == "__main__":
    import sys

    DeformationGridGeoTIFF.main(sys.argv[1:])
