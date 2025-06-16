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
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fill_voids is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with fill_voids.  If not, see <https://www.gnu.org/licenses/>.
*****************************************************************
"""
cimport cython
from libc.stdlib cimport calloc, free
from libc.stdint cimport (
  int8_t, int16_t, int32_t, int64_t,
  uint8_t, uint16_t, uint32_t, uint64_t
)
from libcpp cimport bool as native_bool
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
  cdef size_t binary_fill_holes2d[T](
    T* labels, 
    size_t sx, size_t sy
  )
  cdef size_t binary_fill_holes3d[T](
    T* labels, 
    size_t sx, size_t sy, size_t sz
  )


class DimensionError(Exception):
  pass


@cython.binding(True)
def fill(labels, in_place=False, return_fill_count=False):
  """
  Fills holes in a 1D, 2D, or 3D binary image.

  labels: a binary valued numpy array of any common 
    integer or floating dtype

  in_place: bool, Allow modification of the input array (saves memory)
  return_fill_count: Also return the number of voxels that were filled in.

  Let IMG = a void filled binary image of the same dtype as labels

  if return_fill_count:
    Return: (IMG, number of filled in background voxels)
  else:
    Return: IMG
  """
  ndim = labels.ndim 
  shape = labels.shape 

  if labels.ndim < 2:
    labels = labels[..., np.newaxis]
  while labels.ndim > 3:
    if labels.shape[-1] == 1:
      labels = labels[..., 0]
    else:
      raise DimensionError("The input volume must be (effectively) a 1D, 2D or 3D image: " + str(shape))

  dtype = labels.dtype
  if labels.dtype == bool:
    labels = labels.view(np.uint8)

  if labels.size == 0:
    num_filled = 0
  elif labels.ndim == 2:
    (labels, num_filled) = _fill2d(labels, in_place)
  elif labels.ndim == 3:
    (labels, num_filled) = _fill3d(labels, in_place)
  else:
    raise DimensionError("fill_voids only handles 1D, 2D, and 3D data. Got: " + str(shape))

  while labels.ndim > ndim:
    labels = labels[..., 0]
  while labels.ndim < ndim:
    labels = labels[..., np.newaxis]

  labels = labels.view(dtype)

  if return_fill_count:
    return (labels, num_filled)
  else:
    return labels

def _fill3d(cnp.ndarray[NUMBER, cast=True, ndim=3] labels, in_place=False):
  if not in_place:
    labels = np.copy(labels, order='F')
  else:
    labels = fastremap.asfortranarray(labels)

  dtype = labels.dtype

  cdef size_t num_filled = 0

  if dtype in (np.uint8, np.int8, bool):
    num_filled = binary_fill_holes3d[uint8_t](<uint8_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint16, np.int16):
    num_filled = binary_fill_holes3d[uint16_t](<uint16_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint32, np.int32):
    num_filled = binary_fill_holes3d[uint32_t](<uint32_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype in (np.uint64, np.int64):
    num_filled = binary_fill_holes3d[uint64_t](<uint64_t*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype == np.float32:
    num_filled = binary_fill_holes3d[float](<float*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  elif dtype == np.float64:
    num_filled = binary_fill_holes3d[double](<double*>&labels[0,0,0], labels.shape[0], labels.shape[1], labels.shape[2])
  else:
    raise TypeError("Type {} not supported.".format(dtype))

  return (labels, num_filled)

def _fill2d(cnp.ndarray[NUMBER, cast=True, ndim=2] labels, in_place=False):
  if not in_place:
    labels = np.copy(labels, order='F')
  else:
    labels = fastremap.asfortranarray(labels)

  dtype = labels.dtype

  cdef size_t num_filled = 0

  if dtype in (np.uint8, np.int8, bool):
    num_filled = binary_fill_holes2d[uint8_t](<uint8_t*>&labels[0,0], labels.shape[0], labels.shape[1])
  elif dtype in (np.uint16, np.int16):
    num_filled = binary_fill_holes2d[uint16_t](<uint16_t*>&labels[0,0], labels.shape[0], labels.shape[1])
  elif dtype in (np.uint32, np.int32):
    num_filled = binary_fill_holes2d[uint32_t](<uint32_t*>&labels[0,0], labels.shape[0], labels.shape[1])
  elif dtype in (np.uint64, np.int64):
    num_filled = binary_fill_holes2d[uint64_t](<uint64_t*>&labels[0,0], labels.shape[0], labels.shape[1])
  elif dtype == np.float32:
    num_filled = binary_fill_holes2d[float](<float*>&labels[0,0], labels.shape[0], labels.shape[1])
  elif dtype == np.float64:
    num_filled = binary_fill_holes2d[double](<double*>&labels[0,0], labels.shape[0], labels.shape[1])
  else:
    raise TypeError("Type {} not supported.".format(dtype))

  return (labels, num_filled)

def void_shard():
  """??? what's this ???"""
  print("Play Starcraft 2!")

