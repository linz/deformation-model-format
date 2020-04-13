#!/usr/bin/python3

import fileunittest
import os.path
import sys
import csv
from datetime import datetime, date

testdir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(testdir, "data")
rootdir = os.path.dirname(testdir)
libdir = os.path.join(rootdir, "LINZ", "deformation")
testpt_extension = ".testdata.csv"
sys.path.insert(0, libdir)
from model import TimeFunction


class TimeFunctionTest(fileunittest.TestCase):
    def checkFunction(self, definition, dates, testid=""):
        function = TimeFunction.factory(definition)
        functype = type(function).__name__
        if testid != "":
            testid = testid.strip() + " "
        for date in dates:
            testname = "{0} {1}at {2:.3f}".format(functype, testid, date)
            value = function.valueAt(date)
            self.check(testname, value)

    def test001_date_calcs(self):
        # Check parseTime function
        self.check("parseTime 1", TimeFunction.parseTime("2012-09-13T10:25:36Z"))
        self.check("parseTime 2", TimeFunction.parseTime("1997-12-16T00:00:00Z"))
        self.checkRun("parseTime 3", lambda: TimeFunction.parseTime("1997-12-16"))
        self.checkRun("parseTime 4", lambda: TimeFunction.parseTime("2012-13-09T10:25:36Z"))
        self.checkRun("parseTime 5", lambda: TimeFunction.parseTime("2012-09-13T10:25:36"))
        # Check decimalYear function
        self.check("decimalYear 1", TimeFunction.decimalYear("2012-09-13T10:25:36Z"))
        self.check("decimalYear 2", TimeFunction.decimalYear(datetime(2012, 9, 13, 10, 25, 36)))
        self.check("decimalYear 3", TimeFunction.decimalYear(2019.423))
        self.check("decimalYear 4", TimeFunction.decimalYear(1987))
        self.checkRun("decimalYear 5", lambda: TimeFunction.decimalYear(date(2019, 9, 13)))
        self.checkRun("decimalYear 6", lambda: TimeFunction.decimalYear("2012-09-13T10:25:36"))

    def test002_constant_function(self):
        self.checkFunction({"type": "constant", "parameters": {}}, [1983.0, 2018.9])

    def test003_velocity_function(self):
        self.checkFunction(
            {"type": "velocity", "parameters": {"reference_epoch": "2000-01-01T00:00:00Z"}}, [1987.0, 2000.0, 2013.4], "1"
        )
        self.checkFunction(
            {"type": "velocity", "parameters": {"reference_epoch": "2012-09-13T00:00:00Z"}}, [1987.0, 2000.0, 2013.4], "2"
        )

    def test004_step_function(self):
        self.checkFunction(
            {"type": "step", "parameters": {"step_epoch": "2012-09-13T00:00:00Z"}}, [1987.0, 2012.699, 2017.700, 2013.4]
        )

    def test005_reverse_step_function(self):
        self.checkFunction(
            {"type": "reverse_step", "parameters": {"step_epoch": "2012-09-13T00:00:00Z"}},
            [1987.0, 2012.699, 2017.700, 2013.4],
        )

    def test006_piecewise_function(self):
        dates = [2000.0, 2012.699, 2012.7, 2013.0, 2013.698, 2013.7, 2013.8, 2014.098, 2014.1, 2020.0]
        model = [
            {"epoch": "2012-09-13T00:00:00Z", "scale_factor": 1.5},
            {"epoch": "2013-09-13T00:00:00Z", "scale_factor": 2.5},
            {"epoch": "2014-02-06T00:00:00Z", "scale_factor": 2.1},
        ]
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "zero", "after_last": "zero", "model": model}}, dates, "zz"
        )
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "zero", "after_last": "constant", "model": model}},
            dates,
            "zc",
        )
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "zero", "after_last": "linear", "model": model}}, dates, "zl"
        )
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "constant", "after_last": "zero", "model": model}},
            dates,
            "cz",
        )
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "linear", "after_last": "zero", "model": model}}, dates, "lz"
        )
        self.checkFunction(
            {"type": "piecewise", "parameters": {"before_first": "linear", "after_last": "linear", "model": model[:2]}},
            dates,
            "ll2",
        )

    def test007_exponential_function(self):
        pass


if __name__ == "__main__":
    fileunittest.main()
