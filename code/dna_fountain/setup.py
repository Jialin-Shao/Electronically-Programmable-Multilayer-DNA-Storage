"""
Copyright (C) 2016 Yaniv Erlich
License: GPLv3-or-later. See COPYING file for details.
"""

from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'MyProject',
  ext_modules = cythonize(["*.pyx"], gdb_debug=True),
)
