#!/usr/bin/python3

import fileunittest
import os.path
import sys
import csv

testdir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(testdir, "data")
rootdir = os.path.dirname(testdir)
libdir = os.path.join(rootdir, "LINZ", "deformation")
testpt_extension = ".testdata.csv"
sys.path.insert(0, libdir)
from grid import DeformationGridGeoTIFF


class GeotiffGridTest(fileunittest.TestCase):
    def checkGrid(self, name, gridfile, testfile):
        # print("Running test with " + testfile)
        tfile = os.path.join(datadir, testfile)
        gfile = os.path.join(datadir, gridfile)
        assert os.path.exists(tfile), "Test file {0} is missing".format(tfile)
        assert os.path.exists(gfile), "Grid file {0} is missing".format(gfile)

        try:
            grid = DeformationGridGeoTIFF(gfile, geographic=True)
            with open(tfile) as tfh:
                tfr = csv.reader(tfh)
                tfr.__next__()
                for npt, r in enumerate(tfr):
                    lon, lat = float(r[0]), float(r[1])
                    # print("testing {0} {1}".format(lon, lat))
                    testname = name + " point {0}".format(npt + 1)
                    self.checkRun(testname, lambda: grid.calcValue(lon, lat))
        except Exception as ex:
            self.reportFail("Failed to run test {0}: {1}".format(name, str(ex)))

    def test001_bilinear_interpolation(self):
        self.checkGrid("Bilinear interpolation", "grid2x2.tif", "grid2x2.test001.csv")

    def test002_geodetic_bilinear_interpolation(self):
        self.checkGrid("Geocentric bilinear interpolation", "grid2x2gb.tif", "grid2x2gb.test002.csv")

    def test003_select_grid_cells(self):
        self.checkGrid("Grid cell selection", "grid3x4.tif", "grid3x4.test003.csv")

    def test_004_range_error(self):
        self.checkGrid("Points out of range", "grid3x4.tif", "grid3x4.test004.csv")

    def test_005_polar_grid(self):
        self.checkGrid("Geocentric bilinear near pole", "polar.tif", "polar.test005.csv")

    def test_006_antimeridian(self):
        self.checkGrid("Offset 360 east", "grid2x2_east.tif", "grid2x2.test001.csv")
        self.checkGrid("Offset 360 west", "grid2x2_west.tif", "grid2x2.test001.csv")


if __name__ == "__main__":
    fileunittest.main()
