import sys
import os
import os.path
import re
import json
import inspect
import atexit
import unittest
import numpy as np

# Floating point test tolerance
defaultDelta = 1.0e-8

# Set True to write the results from all tests to stdout
dumpFile = None
dumpDelimiter = "{\n"


class TestCase(unittest.TestCase):
    """
    Subclass of unittest.TestCase to support generation and use of a file
    of test results.  Basic structure of a test module using this is:

        from LINZ import fileunittest

        class MyTestCase( fileunittest.TestCase ):

            def test_001( self ):
                result=calc_result()
                self.check('test 001 result',result)

            def test002( self ):
                if notgood():
                    self.reportFail('Test 2 failed')

        if __name__ == "__main__":
            fileunittest.main()

    Run with the --dump option to create a file with the same name as the 
    test module but extension .results.new  Run without --dump to test against
    file .results.
    """

    resultsFile = None
    dumpFile = None
    testResults = {}

    @classmethod
    def setUpClass(cls):
        global dumpFile

        clsfile = inspect.getfile(cls)
        TestCase.resultsFile = os.path.splitext(clsfile)[0] + ".results"
        try:
            if dumpFile is None:
                with open(TestCase.resultsFile) as rf:
                    TestCase.testResults = json.loads(rf.read())
        except Exception as ex:
            print("Failed to load results from Tescase.resultsFile: {0}".format(ex))
            raise

    def toJsonStr(self, output):
        try:
            jsonstr = json.dumps(output)
        except TypeError:
            jsonstr = json.dumps(repr(output))
        return jsonstr

    def check(self, testname, output, message=None, delta=None):
        """
        Call to check that output matches the value for testname in the .results
        file.  message is the error message if it fails.  delta is a numerical
        tolerance for floating point tests.  delta is used in the unittest.assertAlmostEqual
        call.
        """
        global dumpFile
        global dumpDelimiter
        testcode = testname.lower()
        testcode = re.sub(r"\s+", "_", testcode)
        testcode = re.sub(r"\W", "_", testcode)
        if isinstance(output, np.ndarray):
            output = output.tolist()
        if dumpFile is not None:
            if testcode in self.testResults:
                raise RuntimeError("Invalid duplicate test code: " + testcode)
            with open(dumpFile, "a") as dumph:
                dumph.write(dumpDelimiter)
                dumpDelimiter = ",\n"
                dumph.write('"{0}": {1}'.format(testcode, self.toJsonStr(output)))
            self.testResults[testcode] = output
        else:
            message = message or testname + " incorrect"
            if testcode not in self.testResults:
                self.fail("Test result {0} missing in results file".format(testcode))
            else:
                expected = self.testResults.get(testcode)
                if type(expected) == str and type(output) != str:
                    output = repr(output)
                self.checkEqual(output, expected, message, delta)

    def checkRun(self, testname, function, message=None, delta=None):
        # Run a function with no parameters and either check the output
        # or check the error generated.
        try:
            self.check(testname, function(), message, delta)
        except AssertionError as ex:
            raise
        except Exception as ex:
            self.check(testname, ex, message, delta)

    def checkEqual(self, output, expected, message, delta):
        # Equality check with tolerance on floating point numbers and handling
        # of numpy arrays
        global defaultDelta
        delta = delta or defaultDelta
        if isinstance(output, float) and isinstance(expected, float):
            message = message + " ({0} != {1})".format(output, expected)
            self.assertAlmostEqual(output, expected, msg=message, delta=delta)
        elif isinstance(output, str) and isinstance(expected, str):
            if output != expected:
                self.fail(message + "  ({0} != {1})".format(output, expected))
        elif isinstance(output, dict) and isinstance(expected, dict):
            for k in list(expected.keys()):
                self.checkEqual(output[k], expected[k], message + "[{0}]".format(k), delta)
        elif hasattr(output, "__getitem__") and hasattr(expected, "__getitem__"):
            for i, e in enumerate(expected):
                o = output[i]
                self.checkEqual(o, e, message + "[{0}]".format(i), delta)
        else:
            try:
                self.assertEqual(output, expected, msg=message)
            except Exception as ex:
                self.fail(message + ": Cannot test equality of {0} and {1}: " + str(ex))

    def reportFail(self, message):
        """
        Function to directly report a failed test
        """
        global dumpFile
        if dumpFile is None:
            self.fail(message)


def main():
    """
    Main function called from subclasses to run the tests directly
    """
    global dumpFile
    global dumpDelimiter
    if "--dump" in sys.argv:
        sys.argv.remove("--dump")
        print("**** Dumping output - not running tests!")
        dumpFile = os.path.splitext(os.path.realpath(sys.argv[0]))[0] + ".results.new"
        if os.path.exists(dumpFile):
            os.remove(dumpFile)

        def enddump():
            if dumpDelimiter != "{\n":
                with open(dumpFile, "a") as dumph:
                    dumph.write("\n}\n")

        atexit.register(enddump)
    try:
        unittest.main()
    except Exception as ex:
        print("Failed to run unit tests: {0}".format(ex))
