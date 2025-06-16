import typing
from typing import Literal, Union, overload

import numpy as np
from numpy.typing import NDArray

_T = typing.TypeVar("_T", bound=np.generic)

class DimensionError(Exception): ...

@overload
def fill(
    labels: NDArray[_T],
    in_place: bool = False,
    *,
    return_fill_count: Literal[False] = False,
) -> NDArray[_T]: ...
@overload
def fill(
    labels: NDArray[_T],
    in_place: bool,
    return_fill_count: Literal[False] = False,
) -> NDArray[_T]: ...
@overload
def fill(
    labels: NDArray[_T],
    in_place: bool = False,
    *,
    return_fill_count: Literal[True],
) -> tuple[NDArray[_T], int]: ...
@overload
def fill(
    labels: NDArray[_T],
    in_place: bool,
    return_fill_count: Literal[True],
) -> tuple[NDArray[_T], int]: ...
def fill(  # type: ignore[misc]
    labels: NDArray[_T], in_place: bool = False, return_fill_count: bool = False
) -> Union[NDArray[_T], tuple[NDArray[_T], int]]:
    """Fills holes in a 1D, 2D, or 3D binary image.

    Args:
        labels: a binary valued numpy array of any common
            integer or floating dtype
        in_place: bool, Allow modification of the input array (saves memory)
        return_fill_count: Also return the number of voxels that were filled in.

    Returns:
        A void filled binary image of the same dtype as labels with the number
        of filled in background voxels if return_fill_count is True.
    """

def void_shard() -> None: ...
