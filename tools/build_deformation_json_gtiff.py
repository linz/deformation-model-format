#!/usr/bin/python3

# Script to create the master file GeoTIFF file creation scripts to encode
# the LINZ Deformation model

from __future__ import print_function

from collections import namedtuple, OrderedDict
from datetime import datetime
from os.path import isdir
from os.path import join as joinpath
from subprocess import call
import argparse
import hashlib
import os
import os.path
import re
import json
import sys
from LINZ.DeformationModel import Model, Time


# Arguments used by the program

parser = argparse.ArgumentParser("Build a LINZDEF deformation file from a published deformation model")
parser.add_argument("model_dir", help="Base directory of model, containing model and tools directories")
parser.add_argument("target_dir", help="Target directory in which to build model")
parser.add_argument("model_name", help="Base name of final model")
parser.add_argument("-s", "--source-crs", default="EPSG:4959", help="Source/interpolation CRS for the model")
parser.add_argument("-t", "--target-crs", default="EPSG:7907", help="Target CRS for the model")
parser.add_argument(
    "-l", "--license", default="Creative Commons Attribution 4.0 International", help="License under which model is published"
)
parser.add_argument("-r", "--reference-date", default="2001-01-01", help="Reference date for model")
parser.add_argument("-u", "--default-uncertainty", type=float, default=0.01, help="Default uncertainty for model")
parser.add_argument("-c", "--compact-metadata", action="store_true", help="Reduce size of metadata in GeoTIFF directory")
parser.add_argument("-n", "--no-optimize", action="store_true", help="Don't optimized GeoTIFF for cloud")
parser.add_argument("--no-build-grids", action="store_false", dest="build_grids", help="Don't build GeoTIFF grids")
parser.add_argument("--uint16-encoding", action="store_true", help="Use uint16 storage with linear scaling/offseting")
parser.add_argument("--split-child-grids", action="store_true", help="Split child grids that overlap multiple parents")
parser.add_argument(
    "-i", "--individual-grids", action="store_true", help="Write separate grid files for each child (for testing)"
)
parser.add_argument("-a", "--all-versions", action="store_true", help="Create master files for all versions of model")
parser.add_argument(
    "-k", "--keep-gtiff-definition-files", action="store_true", help="Keep JSON definitions used to build GeoTIFF files"
)
parser.add_argument("-v", "--verbose", action="store_true", help="More diagnostic output")

# Classes used to compile model


class TimeEvent:
    def __init__(self, time, f0, time1, f1):
        self.time0 = time
        self.time1 = time1
        self.days = time.daysAfter(refdate)
        self.f0 = f0
        self.f1 = f1

    def __str__(self):
        return "Time:{0} f0:{1} {2} f1:{3}".format(self.time0, self.f0, self.time1, self.f1)

    def __lt__(self, other):
        return (self.time0, self.time1) < (other.time0, other.time1)

    def __repr__(self):
        return str(self)


class Extent:
    def __init__(self, minlon, minlat, maxlon, maxlat):
        self._minlon = minlon
        self._maxlon = maxlon
        self._minlat = minlat
        self._maxlat = maxlat

    def copy(self):
        return Extent(self._minlon, self._minlat, self._maxlon, self._maxlat)

    def unionWith(self, other):
        if other._minlon < self._minlon:
            self._minlon = other._minlon
        if other._maxlon > self._maxlon:
            self._maxlon = other._maxlon
        if other._minlat < self._minlat:
            self._minlat = other._minlat
        if other._maxlat > self._maxlat:
            self._maxlat = other._maxlat

    def spec(self):
        return [self._minlon, self._minlat, self._maxlon, self._maxlat]


TimeStep = namedtuple("TimeStep", "mtype t0 f0 t1 f1")
DefSeq = namedtuple("DefSeq", "component description gridtype versionstart versionend zerobeyond steps grids extent")
DefComp = namedtuple("DefComp", "date factor before after")

clean_description = lambda d: re.sub(r"[\r\n]+", "\n", d.strip())

timeformat = "%Y-%m-%dT00:00:00Z"

# ==============================================

args = parser.parse_args()
allversions = args.all_versions
verbose = args.verbose
source_crs = args.source_crs
target_crs = args.target_crs
component_uncertainty = [args.default_uncertainty, args.default_uncertainty]
licensetype = args.license
refdate = Time.Time(datetime.strptime(args.reference_date, "%Y-%m-%d"))

md = args.model_dir
if not isdir(md):
    raise RuntimeError("Invalid model directory: " + md)

bd = args.target_dir
if not isdir(bd):
    if os.path.exists(bd):
        raise RuntimeError("Invalid target build directory: " + bd)
    else:
        os.makedirs(bd)

if args.build_grids:
    import deformation_csv_to_gtiff

defname = args.model_name

m = Model.Model(joinpath(md, "model"))
authority = m.metadata("authority")

# Compile the sequences (components) used in the model.  These are compiled from
# the current model components, but assembling all the grids which comprise nested
# grid models.

sequences = []

for c in m.components(allversions=allversions):
    print("Adding source component:", c.name)
    component = c.submodel
    # print component
    tm = c.timeFunction
    mtype = tm.time_function
    if mtype not in ("velocity", "step", "ramp"):
        raise RuntimeError("Cannot handle temporal model type " + mtype)
    # print type(c.timeComponent.model())
    # print mtype,tm.factor0,tm.factor1,tm.time0,tm.time1

    grids = None
    gridtype = "none"
    zerobeyond = "yes"
    if c.versionRevoked != "0" and not allversions:
        continue
    versionstart = c.versionAdded
    versionend = c.versionRevoked
    extent = None
    descriptions = OrderedDict()
    for sm in c.spatialModel.models():
        if sm.spatial_model != "llgrid":
            raise RuntimeError("Cannot handle spatial model type " + sm.spatial_model)
        descriptions[sm.description] = 1
        # print sm.model().gridSpec()
        # print '    ',sm.model().gridFile()
        # print '    ',sm.columns
        gridtype = sm.displacement_type + ":" + sm.error_type
        model = sm.model()
        smextent = [sm.min_lon, sm.min_lat, sm.max_lon, sm.max_lat]
        copyright = "{0} ({1}): Released under {2}".format(authority, versionstart[:4], licensetype)
        if grids is None:
            grids = OrderedDict(
                [
                    ("crs", source_crs),
                    ("type", gridtype),
                    ("copyright", copyright),
                    ("description", sm.description),
                    ("version", versionstart),
                    ("grids", []),
                ]
            )
        grids["grids"].append(
            OrderedDict([("filename", sm.file1), ("extent", smextent), ("size", [sm.npoints1, sm.npoints2])])
        )
        if not sm.spatial_complete:
            zerobeyond = "no"
        mextent = Extent(*smextent)
        if extent is None:
            extent = mextent
        else:
            extent.unionWith(mextent)

    description = "\n".join([d for d in descriptions])
    description = description.replace("\r", "")

    step = TimeStep(mtype, tm.time0, tm.factor0, tm.time1, tm.factor1)

    found = False
    for s in sequences:
        if (
            s.component == component
            and s.zerobeyond == zerobeyond
            and s.versionstart == versionstart
            and s.versionend == versionend
            and s.grids == grids
            and s.gridtype == gridtype
        ):
            s.steps.append(step)
            s.extent.unionWith(extent)
            found = True
            break

    if not found:
        sequences.append(DefSeq(component, description, gridtype, versionstart, versionend, zerobeyond, [step], grids, extent))

# Compile components and versions

SeqComp = namedtuple("SeqComp", "time factor before after nested")
small = 0.00001
gdf_files = {}
seqcomps = {}

for sequence in sequences:
    if sequence.component not in seqcomps:
        seqcomps[sequence.component] = 0
    seqcomps[sequence.component] += 1

seqcomps = {c: "" for c in seqcomps if seqcomps[c] > 1}

components = []
gridnames = set()

for sequence in sequences:
    compname = sequence.component
    ncomp = 0
    while compname in seqcomps:
        ncomp += 1
        compname = sequence.component + "_c{0}".format(ncomp)
    seqcomps[compname] = sequence.component

    print("Building proj component:", compname)

    subsequences = []
    events = []
    timefuncs = []

    for s in sequence.steps:
        if s.mtype == "velocity":
            timefuncs.append(
                (OrderedDict([("type", "velocity"), ("parameters", {"reference_epoch": s.t0.strftime(timeformat)})]), s.t0)
            )
        elif s.mtype == "step":
            events.append(TimeEvent(s.t0, s.f0, s.t0, s.f1))
        elif s.mtype == "ramp":
            events.append(TimeEvent(s.t0, s.f0, s.t1, s.f1))
        else:
            raise RuntimeError("Cannot handle time model type " + s.mtype)

    time_model = []
    eventtime = None
    if len(events) > 0:
        time_model = []
        events.sort()
        eventtime = events[0].time0
        e0 = None
        for e in events:
            if e0 and e0.time1.daysAfter(e.time0) > 0.001:
                raise RuntimeError("Cannot handle overlapping time events in series")
            v0 = 0.0 if len(time_model) == 0 else time_model[-1][1]
            for t in time_model:
                t[1] += e.f0
            time_model.append([e.time0, e.f0 + v0])
            time_model.append([e.time1, e.f1 + v0])
        for i in reversed(range(len(time_model) - 1)):
            t0 = time_model[i]
            t1 = time_model[i + 1]
            if abs(t0[0].daysAfter(t1[0])) < 0.001 and abs(t0[1] - t1[1]) < 0.00001:
                time_model[i : i + 1] = []

    if (
        len(time_model) == 2
        and time_model[0][0] == time_model[1][0]
        and time_model[0][1] == 0.0
        and time_model[1][1] == 1.0
        and time_model[0][0]
    ):
        timefuncs.append(
            (OrderedDict([("type", "step"), ("parameters", {"step_epoch": time_model[0][0].strftime(timeformat)})]), eventtime)
        )
    elif (
        len(time_model) == 2
        and time_model[0][0] == time_model[1][0]
        and time_model[0][1] == -1.0
        and time_model[1][1] == 0.0
        and time_model[0][0]
    ):
        timefuncs.append(
            (
                OrderedDict([("type", "reverse_step"), ("parameters", {"step_epoch": time_model[0][0].strftime(timeformat)})]),
                eventtime,
            )
        )
    elif len(time_model) > 1:
        timefuncs.append(
            (
                OrderedDict(
                    [
                        ("type", "piecewise"),
                        (
                            "parameters",
                            OrderedDict(
                                [
                                    ("before_first", "zero" if time_model[0][1] == 0.0 else "constant"),
                                    ("after_last", "zero" if time_model[-1][1] == 0.0 else "constant"),
                                    (
                                        "model",
                                        [
                                            OrderedDict(
                                                [("epoch", t[0].strftime(timeformat)), ("scale_factor", round(t[1], 4))]
                                            )
                                            for t in time_model
                                        ],
                                    ),
                                ]
                            ),
                        ),
                    ]
                ),
                eventtime,
            )
        )

    # Now construct the sequences in the definition file...

    filespecs = {}
    gridspecs = None

    gkey = ":".join(sorted([g["filename"] for g in sequence.grids["grids"]]))
    if gkey in filespecs:
        gridspecs = filespecs[gkey]
        if verbose:
            print("Using existing grid {0}".format(gridspec))
    else:
        gridspecs = []
        gname = sequence.grids["grids"][0]["filename"]
        gname = os.path.dirname(gname)
        gname = os.path.basename(gname)
        gname = gname.replace("patch_", "")
        gname = gname.replace("_", "")
        gname = defname + "-" + gname
        ngname = 0
        while True:
            ngname += 1
            gridname = gname + "-grid{0:02d}".format(ngname) + ".tif"
            if gridname not in gridnames:
                break
        gridnames.add(gridname)
        if verbose:
            print("Building grid {0}".format(gridname))
        gtiffj = os.path.join(bd, gridname) + ".json"
        with open(gtiffj, "w") as gridf:
            gridf.write(json.dumps(sequence.grids, indent=2))
        if args.build_grids:
            gridfiles = deformation_csv_to_gtiff.build_deformation_gtiff(gtiffj, gridname, args, basedir=bd)
            if not args.keep_gtiff_definition_files:
                os.unlink(gtiffj)
        else:
            gridfiles = [gridname + ".json"]

        for gridfile in gridfiles:
            gtiff = os.path.join(bd, gridfile)
            if not os.path.exists(gtiff):
                raise RuntimeError("Failed to create GeoTIFF {0}".format(gtiff))
            md5 = hashlib.md5()
            with open(gtiff, "rb") as gf:
                while True:
                    buffer = gf.read(2048)
                    if len(buffer) == 0:
                        break
                    md5.update(buffer)
            gridspecs.append(
                OrderedDict(
                    [
                        ("type", "GeoTIFF"),
                        ("interpolation_method", "bilinear"),
                        ("filename", gridfile),
                        ("md5_checksum", md5.hexdigest()),
                    ]
                )
            )
        filespecs[gkey] = gridspecs

    # Add a component for each time function (should only be one, but could be multiple velocities in theory)

    disptype, unctype = sequence.gridtype.split(":")

    for nfunc, functime in enumerate(timefuncs):
        for gridspec in gridspecs:
            components.append(
                (
                    sequence,
                    functime[1],
                    OrderedDict(
                        [
                            ("description", sequence.description),
                            ("displacement_type", disptype),
                            ("uncertainty_type", unctype),
                            ("horizontal_uncertainty", component_uncertainty[0]),
                            ("vertical_uncertainty", component_uncertainty[1]),
                            ("extent", OrderedDict([("type", "bbox"), ("parameters", {"bbox": sequence.extent.spec()})])),
                            ("spatial_model", gridspec),
                            ("time_function", functime[0]),
                        ]
                    ),
                )
            )

# Links to information about the deformation model.

about = OrderedDict(
    [
        ("href", "https://www.linz.govt.nz/nzgd2000"),
        ("rel", "about"),
        ("type", "text/html"),
        ("title", "About the NZGD2000 deformation model"),
    ]
)
source = OrderedDict(
    [
        ("href", "https://www.geodesy.linz.govt.nz/download/nzgd2000_deformation_model"),
        ("rel", "source"),
        ("type", "application/zip"),
        ("title", "Authoritative source of the NZGD2000 deformation model"),
    ]
)
license = OrderedDict(
    [
        ("href", "https://creativecommons.org/licenses/by/4.0/"),
        ("rel", "license"),
        ("type", "text/html"),
        ("title", "Creative Commons Attribution 4.0 International license"),
    ]
)
metadatafunc = lambda v: OrderedDict(
    [
        (
            "href",
            "https://www.geodesy.linz.govt.nz/download/nzgd2000/metadata/nzgd2000_deformation_{version}_metadata.xml".replace(
                "{version}", v
            ),
        ),
        ("rel", "metadata"),
        ("type", "application/xml"),
        ("title", " ISO 19115 XML encoded metadata regarding the deformation model"),
    ]
)
versions = [v for v in m.versions()] if allversions else [m.currentVersion()]
for v in versions:
    vseq = [c for c in components if c[0].versionstart <= v and (c[0].versionend == "0" or c[0].versionend > v)]
    vseq.sort(key=lambda c: c[1])
    if len(vseq) == 0:
        continue
    vcomps = [s[2] for s in vseq]
    vextent = vseq[0][0].extent.copy()
    for s in vseq[1:]:
        vextent.unionWith(s[0].extent)
    modelspec = OrderedDict(
        [
            ("file_type", "deformation_model_master_file"),
            ("format_version", "1.0"),
            ("name", m.metadata("model_name")),
            ("version", v),
            ("publication_date", m.versionInfo(v).release_date.strftime(timeformat)),
            ("license", licensetype),
            ("description", m.metadata("description").replace("\r", "")),
            (
                "authority",
                OrderedDict(
                    [
                        ("name", m.metadata("authority")),
                        ("url", m.metadata("authority_website")),
                        ("address", m.metadata("authority_address")),
                        ("email", m.metadata("authority_email")),
                    ]
                ),
            ),
            ("links", [about, source, license, metadatafunc(v)]),
            ("source_crs", source_crs),
            ("target_crs", target_crs),
            ("definition_crs", source_crs),
            ("reference_epoch", refdate.strftime(timeformat)),
            # This is an arbitrary reference epoch, but more realistic than 2000-01-01!
            ("uncertainty_reference_epoch", datetime(2018, 12, 1).strftime(timeformat)),
            ("horizontal_offset_unit", "metre"),
            ("vertical_offset_unit", "metre"),
            ("horizontal_uncertainty_type", "circular 95% confidence limit"),
            ("horizontal_uncertainty_unit", "metre"),
            ("vertical_uncertainty_type", "95% confidence limit"),
            ("vertical_uncertainty_unit", "metre"),
            ("horizontal_offset_method", "addition"),
            ("extent", OrderedDict([("type", "bbox"), ("parameters", {"bbox": vextent.spec()})])),
            (
                "time_extent",
                OrderedDict(
                    [("first", datetime(1900, 1, 1).strftime(timeformat)), ("last", datetime(2050, 1, 1).strftime(timeformat))]
                ),
            ),
            ("components", vcomps),
        ]
    )
    modeljson = json.dumps(modelspec, indent=2)
    mergefunc = lambda m: m.group(1) + re.sub(r"\s+", "", m.group(2))
    modeljson = re.sub(r"(\"bbox\"\:\s)([^\]]*\])", mergefunc, modeljson)
    modeljson = re.sub(r"(\"attributes\"\:\s)([^\]]*\])", mergefunc, modeljson)
    modeljson = re.sub(r"(\"default_uncertainty\"\:\s)([^\]]*\])", mergefunc, modeljson)
    modeljson = re.sub(r"(\{)(\s*\"epoch\"[^\}]+\"scale_factor\"[^\}]+)", mergefunc, modeljson)
    deffile = defname + "-" + v + ".json"
    with open(os.path.join(bd, deffile), "w") as dfh:
        dfh.write(modeljson)
    print("Created proj deformation master file {0}".format(deffile))
