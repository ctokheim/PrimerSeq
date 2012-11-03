import matplotlib
from distutils.core import setup
from py2exe import *
import sys

sys.argv.append('py2exe')
includes = ["matplotlib.backends.backend_tkagg"]
excludes = []
setup(
    windows = [{'script': 'PrimerApp.py'}],
    options={"py2exe": {"includes": includes, "excludes": excludes}},
    data_files=matplotlib.get_py2exe_datafiles(),
)