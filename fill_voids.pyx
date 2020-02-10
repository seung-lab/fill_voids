"""
Fill in holes in 3D binary images. 

Similar to scipy.morphology.binary_fill_holes
but designed to be more performant.

The power of this library derives from the emptiness of space and the
latent energies of the void.

Author: William Silversmith
Affiliation: Seung Lab, Princeton Neuroscience Institute
Date: December 2019

*****************************************************************
This file is part of fill_voids.

fill_voids is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fill_voids is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fill_voids.  If not, see <https://www.gnu.org/licenses/>.
*****************************************************************
"""
cimport cython
from libc.stdlib cimport calloc, free
from libc.stdint cimport (
  int8_t, int16_t, int32_t, int64_t,
  uint8_t, uint16_t, uint32_t, uint64_t
)
from libcpp cimport bool
from cpython cimport array 
import array
import sys

from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map
from libcpp.utility cimport pair as cpp_pair

cimport numpy as cnp
import numpy as np

import fastremap

cdef extern from "math.h":
  float INFINITY

ctypedef fused NUMBER: 
  int8_t
  int16_t
  int32_t
  int64_t
  uint8_t
  uint16_t
  uint32_t
  uint64_t
  unsigned char
  float 
  double

cdef extern from "fill_voids.hpp" namespace "fill_voids":
  cdef void binary_fill_holes[T](
    T* labels, 
    size_t sx, size_t sy, size_t sz
  )

def fill(labels, in_place=False):
  """
  fill(cnp.ndarray[NUMBER, cast=True, ndim=3] labels, in_place=False)

  labels: a binary valued numpy array of any common 
    integer or floating dtype

  in_place: bool, Allow modification of the input array (saves memory)

  Return: a void filled binary image of the same dtype
  """
  if labels.size == 0:
    return labels
  return _fill(labels, in_place)

def _fill(cnp.ndarray[NUMBER, cast=True, ndim=3] labels, in_place=False):
  if not in_place:
    labels = np.copy(labels, order='F')
  else:
    label = fastremap.asfortranarray(labels)

  dtype = labels.dtype

  if dtype in (np.uint8, np.int8, np.bool):
    binary_fill_holes[uint8_t](<uint8_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint16, np.int16):
    binary_fill_holes[uint16_t](<uint16_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint32, np.int32):
    binary_fill_holes[uint32_t](<uint32_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint64, np.int64):
    binary_fill_holes[uint64_t](<uint64_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype == np.float32:
    binary_fill_holes[float](<float*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype == np.float64:
    binary_fill_holes[double](<double*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  else:
    raise TypeError("Type {} not supported.".format(dtype))

  return labels

def void_shard():
  """??? what's this ???"""
  print("Play Starcraft 2!")

