import fill_voids
from scipy.ndimage.morphology import binary_fill_holes

from tqdm import tqdm

import numpy as np 

img = np.load('test_data.npy')
SEGIDS = np.unique(img)[1:]

def test_scipy_comparison():
  segids = np.copy(SEGIDS)
  np.random.shuffle(segids)

  for segid in tqdm(segids[:5]):
    binimg = img == segid
    fv = fill_voids.fill(binimg, in_place=True)
    spy = binary_fill_holes(binimg)

    assert np.all(fv == spy)

def test_dtypes():
  DTYPES = (
    np.bool, np.int8, np.uint8, np.uint16, np.int16, 
    np.int32, np.uint32, np.int64, np.uint64, 
    np.float32, np.float64
  )

  binimg = img == SEGIDS[0]

  comparison = fill_voids.fill(binimg, in_place=False)

  for dtype in DTYPES:
    res = fill_voids.fill(binimg.astype(dtype), in_place=False)
    assert np.all(comparison == res)






