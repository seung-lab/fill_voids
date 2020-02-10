import pytest 

import fill_voids
import scipy.ndimage
from scipy.ndimage.morphology import binary_fill_holes

from tqdm import tqdm

import numpy as np 

img = np.load('test_data.npy')
SEGIDS = np.unique(img)[1:]

def test_scipy_comparison():
  segids = np.copy(SEGIDS)
  np.random.shuffle(segids)

  for segid in tqdm(segids[:10]):
    print(segid)
    binimg = img == segid
    slices = scipy.ndimage.find_objects(binimg)[0]
    binimg = binimg[slices]

    orig_binimg = np.copy(binimg, order='F')
    fv = fill_voids.fill(binimg, in_place=False)
    fvip = fill_voids.fill(binimg, in_place=True)

    assert np.all(fv == fvip)

    spy = binary_fill_holes(binimg)

    assert np.all(fv == spy)

DTYPES = (
  np.bool, np.int8, np.uint8, np.uint16, np.int16, 
  np.int32, np.uint32, np.int64, np.uint64, 
  np.float32, np.float64
)

@pytest.mark.parametrize("dtype", DTYPES)
def test_dtypes(dtype):
  binimg = img == SEGIDS[0]

  comparison = fill_voids.fill(binimg, in_place=False)
  res = fill_voids.fill(binimg.astype(dtype), in_place=False)
  assert np.all(comparison == res)

def test_zero_array():
  labels = np.zeros((0,), dtype=np.uint8)
  # just don't throw an exception
  fill_voids.fill(labels, in_place=False)
  fill_voids.fill(labels, in_place=True)

  labels = np.zeros((128,128,128), dtype=np.uint8)
  fill_voids.fill(labels, in_place=True)
  assert not np.any(labels)

def test_return_count():
  labels = np.ones((10, 10, 10), dtype=bool)
  labels[3:6,3:6,3:6] = False

  filled = fill_voids.fill(labels)
  assert np.all(filled == 1)

  filled, ct = fill_voids.fill(labels, return_fill_count=True)
  assert np.any(labels == False)
  assert ct == 27

@pytest.mark.parametrize("dimension", [1,2,3,4,5,6])
def test_dimensions(dimension):
  size = [5] * dimension
  for i in range(3, dimension):
    size[i] = 1

  labels = np.ones(size, dtype=np.uint8)
  labels = fill_voids.fill(labels)
  assert labels.ndim == dimension

  if dimension <= 3:
    return
    
  size[dimension - 1] = 2
  labels = np.ones(size, dtype=np.uint8)
  try:
    labels = fill_voids.fill(labels)
    assert False 
  except fill_voids.DimensionError:
    pass
