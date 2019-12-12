#!/usr/bin/env python
import setuptools

import numpy as np

# NOTE: If fill_voids.cpp does not exist, you must run
# cython -3 --cplus fill_voids.pyx

setuptools.setup(
  setup_requires=['pbr', 'numpy'],
  extras_require={
     ':python_version == "2.7"': ['futures'],
     'test': ['scipy'],
  },
  ext_modules=[
    setuptools.Extension(
      'fill_voids',
      sources=[ 'fill_voids.cpp' ],
      language='c++',
      include_dirs=[ np.get_include() ],
      extra_compile_args=[
        '-std=c++11', '-O3', '-ffast-math'
      ]
    ),
  ],
  long_description_content_type='text/markdown',
  pbr=True,
)

