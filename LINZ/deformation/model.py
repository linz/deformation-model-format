#!/usr/bin/python3
#
# DeformationModelJson: Script to load and calculate a deformation model formatted
# with a JSON master file and GeoTIFF grid files.

import hashlib
import json
import math
import numpy as np
import os
import os.path
import re
import sys
from collections import namedtuple
from datetime import datetime
from grid import DeformationGridGeoTIFF, RangeError
from dict_object import DictObject, Field, FormatError

DisplacementFields = {
    "none": [],
    "horizontal": ["east_offset", "north_offset"],
    "vertical": ["vertical_offset"],
    "3d": ["east_offset", "north_offset", "vertical_offset"],
}

UncertaintyFields = {
    "none": [],
    "horizontal": ["horizontal_uncertainty"],
    "vertical": ["vertical_uncertainty"],
    "3d": ["horizontal_uncertainty", "vertical_uncertainty"],
}

DegreeFields = ["east_offset", "north_offset"]
BilinearMethod = "bilinear"
BilinearGeocentricMethod = "geocentric_bilinear"
BilinearGeocentricTypes = ["horizontal", "3d"]

E_OFFSET_ORDINATE = 0
N_OFFSET_ORDINATE = 1
U_OFFSET_ORDINATE = 2
H_UNCERTAINTY_ORDINATE = 3
V_UNCERTAINTY_ORDINATE = 4

VarianceIndices = np.array([H_UNCERTAINTY_ORDINATE, V_UNCERTAINTY_ORDINATE])

DisplacementIndices = {
    "none": [],
    "horizontal": [E_OFFSET_ORDINATE, N_OFFSET_ORDINATE],
    "vertical": [U_OFFSET_ORDINATE],
    "3d": [E_OFFSET_ORDINATE, N_OFFSET_ORDINATE, U_OFFSET_ORDINATE],
}

UncertaintyIndices = {
    "none": [],
    "horizontal": [H_UNCERTAINTY_ORDINATE],
    "vertical": [V_UNCERTAINTY_ORDINATE],
    "3d": [H_UNCERTAINTY_ORDINATE, V_UNCERTAINTY_ORDINATE],
}


class Context:
    def __init__(self, sourcefile, check=False):
        self.sourcefile = sourcefile
        self.basedir = os.path.dirname(sourcefile)
        self.check = check
        self.errors = []

    def raiseError(self, error):
        if self.check:
            self.errors.append(error)
        else:
            raise ValueError(error)


class TimeFunction:

    types = dict()

    @staticmethod
    def factory(value, context=None):
        funcdef = DictObject([Field("type", TimeFunction.types), Field("parameters", dict)], value, context=context)
        result = funcdef.type(funcdef.parameters, context=context)
        return result

    @staticmethod
    def parseTime(value):
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")

    @staticmethod
    def decimalYear(value):
        if type(value) == str:
            value = TimeFunction.parseTime(value)
        if type(value) == datetime:
            r0 = datetime(value.year, 1, 1)
            r1 = datetime(value.year + 1, 1, 1)
            value = value.year + (value - r0).total_seconds() / (r1 - r0).total_seconds()
        if type(value) != float:
            raise ValueError("Invalid value {0} for date/time".format(value))
        return value

    def valueAt(self, epoch):
        raise RuntimeError("TimeFunction.valueAt must be overridden")


class ConstantFunction(TimeFunction):
    def __init__(self, value, context=None):
        pass

    def valueAt(self, epoch):
        return 1.0


TimeFunction.types["constant"] = ConstantFunction


class VelocityFunction(DictObject, TimeFunction):

    fields = [Field("reference_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def valueAt(self, epoch):
        return epoch - self.reference_epoch


TimeFunction.types["velocity"] = VelocityFunction


class StepFunction(DictObject):

    fields = [Field("step_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def valueAt(self, epoch):
        return 1.0 if epoch >= self.step_epoch else 0.0


TimeFunction.types["step"] = StepFunction


class ReverseStepFunction(DictObject):

    fields = [Field("step_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def valueAt(self, epoch):
        return -1.0 if epoch < self.step_epoch else 0.0


TimeFunction.types["reverse_step"] = ReverseStepFunction


class PiecewiseFunctionPoint(DictObject):

    fields = [Field("epoch", TimeFunction.decimalYear), Field("scale_factor", float)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class PiecewiseFunction(DictObject):

    ZERO = 0
    CONSTANT = 1
    LINEAR = 2

    Extrapolation = {"zero": ZERO, "constant": CONSTANT, "linear": LINEAR}

    fields = [
        Field("model", [PiecewiseFunctionPoint]),
        Field("before_first", lambda v: PiecewiseFunction.Extrapolation[v]),
        Field("after_last", lambda v: PiecewiseFunction.Extrapolation[v]),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        if len(self.model) < 1:
            raise ValueError("Piecewise function must contain at least one function point")
        for p0, p1 in zip(self.model[:-1], self.model[1:]):
            if p1.epoch < p0.epoch:
                raise ValueError("Piecewise function epoch {0} cannot be less than {1}".format(p1.epoch, p0.epoch))
        model = self.model
        if self.before_first == self.LINEAR and (len(model) < 2 or model[0].epoch == model[1].epoch):
            raise ValueError("Cannot use linear extrapolation on first two points of piecewise function")
        if self.after_last == self.LINEAR and (len(model) < 2 or model[-2].epoch == model[-1].epoch):
            raise ValueError("Cannot use linear extrapolation on last two points of piecewise function")

    def linearExtrapolation(self, r0, r1, epoch):
        if r0.epoch == r1.epoch:
            return r1.scale_factor
        return ((epoch - r0.epoch) * r1.scale_factor + (r1.epoch - epoch) * r0.scale_factor) / (r1.epoch - r0.epoch)

    def valueAt(self, epoch):
        if epoch <= self.model[0].epoch:
            if self.before_first == self.ZERO:
                return 0.0
            elif self.before_first == self.CONSTANT:
                return self.model[0].scale_factor
            else:
                return self.linearExtrapolation(self.refpoint[0], self.refpoint[1], epoch)
        elif epoch > self.model[1].epoch:
            if self.after_last == self.ZERO:
                return 0.0
            elif self.after_last == self.CONSTANT:
                return self.model[-1].scale_factor
            else:
                return self.linearExtrapolation(self.refpoint[-2], self.refpoint[-1], epoch)
        for r0, r1 in zip(self.model[:-1], self.model[1:]):
            if epoch <= r1.epoch:
                return self.linearExtrapolation(r0, r1, epoch)
        raise RuntimeError("Failed to find epoch {0} in piecewise time function".format(epoch))


TimeFunction.types["piecewise"] = PiecewiseFunction


class ExponentialFunction(DictObject, TimeFunction):

    fields = [
        Field("reference_epoch", TimeFunction.decimalYear),
        Field("end_epoch", TimeFunction.decimalYear, optional=True),
        Field("relaxation_constant", float, min=0.0),
        Field("before_scale_factor", float),
        Field("initial_scale_factor", float),
        Field("final_scale_factor", float),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def valueAt(self, epoch):
        if self.end_epoch is not None and epoch > self.end_epoch:
            epoch = self.end_epoch
        if epoch < self.start_epoch:
            return self.before_scale_factor
        return self.initial_scale_factor + (
            (self.final_scale_factor - self.initial_scale_factor)
            * (1.0 - math.exp((self.ref_epoch - epoch) / self.relaxation_constant))
        )


TimeFunction.types["exponential"] = ExponentialFunction


class TimeExtent(DictObject):

    fields = [Field("first", TimeFunction.decimalYear), Field("last", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def contains(self, epoch):
        return epoch >= self.first and epoch < self.last


class Extent:
    def __init__(self, value, context=None):
        if "type" not in value or value["type"] != "bbox":
            raise ValueError("Invalid extent type {0} - must be bbox".format(value.get(type, "undefined")))
        if type(value.get("parameters")) != dict or "bbox" not in value["parameters"]:
            raise ValueError("Invalid extent - requires bbox parameter")
        try:
            extents = value["parameters"]["bbox"]
            xmin, ymin, xmax, ymax = (float(x) for x in extents)
            if xmax < xmin or ymax < ymin:
                raise ValueError("Error in bounding box - max value less than min value")
            self.xmin = xmin
            self.xmax = xmax
            self.ymin = ymin
            self.ymax = ymax
        except:
            raise ValueError("Invalid bounding box - must be [xmin,ymin,xmax,ymax]")

    def contains(self, x, y):
        if x < self.xmin or x >= self.xmax:
            return False
        if y < self.ymin or y >= self.ymax:
            return False
        return True


class SourceFile:
    def __init__(self, filename, md5=None, context=None):
        self.filename = filename
        self.basedir = context.basedir if context is not None else ""
        self.filepath = os.path.join(self.basedir, self.filename)
        if not os.path.isfile(self.filepath):
            raise ValueError("File {0} is missing".format(self.filepath))
        if context is not None and context.check and md5 is not None:
            if self.calcMd5() != md5:
                raise ValueError("File {0} md5 checksum is not correct")

    def calcMd5(self):
        md5 = hashlib.md5()
        with open(self.filepath, "rb") as gf:
            while True:
                buffer = gf.read(2048)
                if len(buffer) == 0:
                    break
                md5.update(buffer)
        return md5.hexdigest()


class SpatialModel(DictObject):

    fields = [
        Field("type", str, regex="^GeoTIFF$"),
        Field("interpolation_method", str, regex="^(bilinear|geocentric_bilinear)$"),
        Field("filename", str),
        Field("md5_checksum", str),
    ]

    @staticmethod
    def zero():
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        sourcefile = SourceFile(self.filename, self.md5_checksum, context)
        self.grid = DeformationGridGeoTIFF(sourcefile.filepath)
        self.loaded = False

    def load(self):
        if not self.loaded:
            self.loaded = True
            self.grid.load()
            dispindices = DisplacementIndices.get(self.grid.displacement_type().lower(), [])
            uncindices = UncertaintyIndices.get(self.grid.uncertainty_type().lower(), [])
            self.gridindices = dispindices + uncindices

    def check(self, displacement_type, uncertainty_type, context):
        self.load()
        grid_disptype = self.grid.displacement_type().lower()
        grid_unctype = self.grid.uncertainty_type().lower()

        if grid_disptype != displacement_type:
            raise ValueError("Grid displacement type {0} doesn't match expected {1}".format(grid_disptype, displacement_type))
        if grid_unctype != uncertainty_type:
            raise ValueError("Grid uncertainty type {0} doesn't match expected {1}".format(grid_unctype, uncertainty_type))

    def valueAt(self, lon, lat):
        self.load()
        gridvalue = self.grid.calcValue(lon, lat)
        result = self.zero()
        result[self.gridindices] = gridvalue
        return result


class Component(DictObject):

    displacementTypes = {t: t for t in DisplacementFields}
    uncertaintyTypes = {t: t for t in UncertaintyFields}

    fields = [
        Field("description", str),
        Field("extent", Extent),
        Field("displacement_type", displacementTypes),
        Field("uncertainty_type", uncertaintyTypes),
        Field("horizontal_uncertainty", float, min=0.0),
        Field("vertical_uncertainty", float, min=0.0),
        Field("spatial_model", SpatialModel),
        Field("time_function", TimeFunction.factory),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        if context is not None and context.check:
            self.spatial_model.check(self.displacement_type, self.uncertainty_type, context)
        self.epoch = None
        self.tfunc = None

    def _zero(self):
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    def setUncertainty(self, result):
        if self.uncertainty_type == "none":
            result[H_UNCERTAINTY_ORDINATE] = self.horizontal_uncertainty
            result[V_UNCERTAINTY_ORDINATE] = self.vertical_uncertainty
        elif self.uncertainty_type == "horizontal":
            result[V_UNCERTAINTY_ORDINATE] = self.vertical_uncertainty
        elif self.uncertainty_type == "vertical":
            result[H_UNCERTAINTY_ORDINATE] = self.horizontal_uncertainty

    def timeFunctionValueAt(self, epoch):
        if epoch != self.epoch:
            self.tfunc = self.time_function.valueAt(epoch)
        return self.tfunc

    def spatialModelValueAt(self, lon, lat):
        result = self.spatial_model.valueAt(lon, lat)
        self.setUncertainty(result)
        return result

    def valueAt(self, x, y, epoch):
        if not self.extent.contains(x, y):
            result = self._zero()
        else:

            if self.tfunc == 0.0:
                result = self._zero()
            else:
                result = self.spatialModelValueAt(x, y) * self.tfunc
        return result


class Authority(DictObject):
    fields = [
        Field("name", str, regex=r".*\S"),
        Field("url", str, regex=r".*\S"),
        Field("address", str, regex=r".*\S"),
        Field("email", str, regex=r"^.*\@.*"),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class Link(DictObject):

    fields = [
        Field("href", str, regex=r"^https?\:\/\/"),
        Field("rel", str, regex=r"^(about|source|metadata|license)$"),
        Field("type", str, regex=r"^(text\/html|application\/zip|application\/xml)$"),
        Field("title", str),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class Evaluator:

    ScaledComponent = namedtuple("ScaledComponent", "scale component")

    def __init__(self, model):
        self.model = model
        self.epoch = None
        self.scaledComponents = []

    def setEpoch(self, epoch):
        if epoch == self.epoch:
            return
        self.scaleComponents = []
        self.epoch = epoch
        for component in self.model.components:
            scale = component.timeFunctionValueAt(epoch)
            if scale != 0.0:
                self.scaledComponents.append(self.ScaledComponent(scale, component))

    def valueAt(self, lon, lat):
        result = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        for sc in self.scaledComponents:
            try:
                compValue = sc.component.spatialModelValueAt(lon, lat)
                # Convert error to variance
                compValue[VarianceIndices] *= compValue[VarianceIndices]
                result += compValue
            except RangeError:
                continue
        result[VarianceIndices] = np.sqrt(result[VarianceIndices])
        return result


class DeformationModel(DictObject):
    fields = [
        Field("file_type", str, regex=r"^deformation_model_master_file$"),
        Field("format_version", str, regex=r"1.0"),
        Field("name", str),
        Field("version", str),
        Field("publication_date", TimeFunction.parseTime),
        Field("license", str),
        Field("description", str),
        Field("authority", Authority),
        Field("links", [Link], optional=True),
        Field("source_crs", str, regex=r"^EPSG\:\d+$"),
        Field("target_crs", str, regex=r"^EPSG\:\d+$"),
        Field("definition_crs", str, regex=r"^EPSG\:\d+$"),
        Field("reference_epoch", TimeFunction.decimalYear),
        Field("uncertainty_reference_epoch", TimeFunction.decimalYear),
        Field("horizontal_offset_unit", str, regex=r"^metre|degree$"),
        Field("vertical_offset_unit", str, regex=r"^metre$"),
        Field("horizontal_uncertainty_type", str, regex=r"^circular 95% confidence limit$"),
        Field("horizontal_uncertainty_unit", str, regex=r"^metre$"),
        Field("vertical_uncertainty_type", str, regex=r"^95% confidence limit$"),
        Field("vertical_uncertainty_unit", str, regex=r"^metre$"),
        Field("horizontal_offset_method", str, regex=r"^addition$"),
        Field("extent", Extent),
        Field("time_extent", TimeExtent),
        Field("components", [Component]),
    ]

    @staticmethod
    def LoadJson(sourcefile, check=False):
        with open(sourcefile, "r") as jsonf:
            value = json.load(jsonf)
        context = Context(sourcefile, check)
        model = DeformationModel(value, context)
        return model

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        if context.check:
            if self.definition_crs != self.source_crs:
                raise ValueError(
                    "Source CRS ({0}) and definition CRS ({1}) must be the same:".format(self.source_crs, self.definition_crs)
                )
        self.isgeographic = None
        self.ellipsoid_a = None
        self.ellipsoid_b = None
        self.a2 = None
        self.b2 = None
        self.evaluator = Evaluator(self)
        self.evaluator0 = None

    def valueAt(self, lon, lat, epoch):
        if not self.extent.contains(lon, lat):
            raise RangeError()
        self.evaluator.setEpoch(epoch)
        return self.evaluator.valueAt(lon, lat)

    def setCrsEllipsoid(self, a, rf=None):
        """
        Define the ellipsoid used by the source, target, and definition CRS.  
        Required to apply the offset to coordinates
        """
        self.isgeographic = True
        self.ellipsoid_a = a
        if rf is None:
            self.ellipsoid_b = a
        else:
            self.ellipsoid_b = a * (1.0 - 1.0 / rf)
        self.a2 = self.ellipsoid_a * self.ellipsoid_a
        self.b2 = self.ellipsoid_b * self.ellipsoid_b

    def setCrsIsProjection(self):
        """
        Defined the source, target, and definition CRS to be in projections coordinates.
        Units are assumed to be metres.
        """
        self.isgeographic = False

    def enScaleFactors(self, lon, lat):
        if self.isgeographic is None:
            raise RuntimeError("CRS information not available - use setCrsEllipsoid or setCrsIsProjection")
        if not self.isgeographic:
            return 1.0, 1.0
        slat = math.sin(math.radians(lat))
        slat2 = slat * slat
        clat = math.cos(math.radians(lat))
        clat2 = clat * clat
        a2c2 = self.a2 * clat2
        b2s2 = self.b2 * slat2
        bsac = b2s2 + a2c2
        dlonde = math.sqrt(bsac) / a2c2
        dlatdn = bsac / (self.a2 * self.b2)
        return dlonde, dlatdn


def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(description="Calculate deformation in a deformation model JSON file")
    parser.add_argument("json_file", help="Name of the deformation model JSON file")
    parser.add_argument("-p", "--point", nargs=3, help="lon,lat,epoch to evaluate")
    parser.add_argument("-i", "--point-file", help="Name of CSV file - columns id,lon,lat,epoch")
    parser.add_argument("-o", "--output-file", help="Name of CSV output file")
    parser.add_argument("-c", "--check", action="store_true", help="Check the deformation model")
    args = parser.parse_args()
    model = DeformationModel.LoadJson(args.json_file, check=args.check)
    if args.point:
        lon = float(args.point[0])
        lat = float(args.point[1])
        epoch = float(args.point[2])
        print(model.valueAt(lon, lat, epoch))
    if args.point_file:
        import csv

        csvi = csv.DictReader(open(args.point_file, "r"))
        for field in ("id", "lon", "lat", "epoch"):
            if field not in csvi.fieldnames:
                raise RuntimeError("Input file {0} is missing field {1}".format(args.point_file, field))
        if args.output_file:
            csvo = csv.writer(open(args.output_file, "w"))
        else:
            csvo = csv.writer(sys.stdout)
        csvo.writerow(["id", "lon", "lat", "epoch", "de", "dn", "du", "eh", "ev"])
        nrow = 1
        for row in csvi:
            nrow += 1
            try:
                id = row["id"]
                lon = float(row["lon"])
                lat = float(row["lat"])
                epochstr = row["epoch"]
                try:
                    epoch = float(epochstr)
                except:
                    epoch = TimeFunction.decimalYear(epochstr)
                result = model.valueAt(lon, lat, epoch)
                output = [id, str(lon), str(lat), epochstr] + [str(r) for r in result]
                csvo.writerow(output)
            except Exception as ex:
                print("Error in row {0}: {1}".format(nrow, ex))


if __name__ == "__main__":
    main()

