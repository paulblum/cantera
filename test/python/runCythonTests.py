"""
Unit tests for Cantera's Cython-based Python module.

This script gathers all the tests defined 'cantera.test' module, runs them,
and prints a report. Extra command line arguments can be used to run subsets
of the test suite, e.g.:

all tests from 'test_thermo.py' and 'test_kinetics.py':

    python runCythonTests.py thermo kinetics

all tests from the 'test_reactor.TesTIdealGasReactor' class:

    python runCythonTests.py reactor.TestIdealGasReactor

a single test:

    python runCythonTests.py onedim.TestDiffusionFlame.test_mixture_averaged
"""

import sys
import os
from pathlib import Path

cantera_root = os.path.relpath(__file__).split(os.sep)[:-1] + ['..', '..']
os.chdir(os.sep.join(cantera_root + ['test', 'work']))

import unittest
try:
    import pytest
except ImportError:
    pytest = None
import cantera
import cantera.test

class TestResult(unittest.TextTestResult):
    def __init__(self, *args, **kwargs):
        unittest.TextTestResult.__init__(self, *args, **kwargs)
        self.outName = 'python-results.txt'
        with open(self.outName, 'w') as f:
            pass # just create an empty output file

    def reformat(self, test_string):
        name, cls = test_string.split()
        cls = cls.replace('(cantera.test.', '').replace(')','')
        return '%s.%s' % (cls, name)

    def addSuccess(self, test):
        with open(self.outName, 'a') as f:
            f.write('PASS: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addSuccess(self, test)

    def addFailure(self, test, err):
        with open(self.outName, 'a') as f:
            f.write('FAIL: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)

    def addError(self, test, err):
        with open(self.outName, 'a') as f:
            f.write('ERROR: %s\n' % self.reformat(str(test)))
        unittest.TextTestResult.addFailure(self, test, err)


if __name__ == '__main__':
    print('\n* INFO: using Cantera module found at this location:')
    print('*     ', repr(cantera.__file__))
    print('* INFO: Cantera version:', cantera.__version__)
    print('* INFO: Git commit:', cantera.__git_commit__, '\n')
    sys.stdout.flush()

    if len(sys.argv) > 1 and sys.argv[1] == "fast_fail":
        fast_fail = True
        subset_start = 2
    else:
        fast_fail = False
        subset_start = 1

    if pytest is not None:
        base = Path(cantera.__file__).parent.joinpath('test')
        subsets = []
        for name in sys.argv[subset_start:]:
            subsets.append(str(base.joinpath(f"test_{name}.py")))

        if not subsets:
            subsets.append(str(base))

        pytest_args = ["-raP", "--durations=50", "--junitxml=pytest.xml"]
        if fast_fail:
            pytest_args.insert(0, "-x")

        ret_code = pytest.main(pytest_args + subsets)
        sys.exit(ret_code)
    else:
        loader = unittest.TestLoader()
        runner = unittest.TextTestRunner(
            verbosity=2, resultclass=TestResult, failfast=fast_fail
        )
        suite = unittest.TestSuite()
        subsets = []
        for name in sys.argv[subset_start:]:
            subsets.append('cantera.test.test_' + name)

        if not subsets:
            subsets.append('cantera.test')

        suite = loader.loadTestsFromNames(subsets)

        results = runner.run(suite)
        sys.exit(len(results.errors) + len(results.failures))
