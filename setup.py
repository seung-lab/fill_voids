#!/usr/bin/env python
import setuptools

class NumpyImport:
  def __repr__(self):
    import numpy as np

    return np.get_include()

  __fspath__ = __repr__

# NOTE: If fill_voids.cpp does not exist, you must run
# cython -3 --cplus fill_voids.pyx

setuptools.setup(
  setup_requires=['pbr', 'numpy', 'cython'],
  extras_require={
     'test': ['scipy'],
  },
  ext_modules=[
    setuptools.Extension(
      'fill_voids.fill_voids',
      sources=[ 'fill_voids/fill_voids.pyx' ],
      language='c++',
      include_dirs=[ str(NumpyImport()) ],
      extra_compile_args=[
        '-std=c++11', '-O3'
      ]
    ),
  ],
  packages=setuptools.find_packages(),
  long_description_content_type='text/markdown',
  pbr=True,
)

