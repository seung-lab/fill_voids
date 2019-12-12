# fill_voids: Faster binary_fill_holes
```python
import fill_voids

img = ... # 3d binary image 
filled_image = fill_voids.fill(img, in_place=False) # in_place allows editing of original image

```

3D void filling algorithm, similar function as scipy.ndimage.morphology.binary_fill_holes

The purpose of this repo is to make it convenient to improve upon the scipy hole filling 
algorithm which only applies to binary images and uses slow serial dilations. 

The current version of the improvements is to use a flood fill of the background and 
label everything not filled as foreground. This is already significantly faster, but 
we can do better by reducing memory usage, increasing speed further, and supporting multiple 
labels.

### Current version

1. Allocate a uint8 zeroed 3D image that is 2 voxels larger than the binary input image along each dimension.
2. Paint the input image into the new image such that there is a black border around it setting foreground as 2.
3. Flood fill (six connected) from a spot on the border (arbitarily taken to be (0,0,0)) and label the flood fill as 1.
4. Write out a binary image the same size as the input from the temporary buffer where foreground is set as buffer != 1 (i.e. 0 or 2).

We improve performance significantly by using libdivide to make computing x,y,z coordinates from array index faster, and also by testing whether a neighbor has been visited or contains foreground already before putting it on the stack.

### Binary Version Improvements  

It would be possible to skip the allocating and painting steps by walking along each side of the image and adding newly encountered voids as seed points to the stack. This reduces the memory usage to near zero.

### Multi-Label Improvements 

Similarly to the connected-components-3d and euclidean-distance-3d projects, in connectomics, it can be common to want to apply void filling algorithms to all labels within a densely packed volume. A multi-label algorithm can be much faster than even the fastest serial application of a binary algorithm. Here's how this might go given an input image I:


1. Compute M = max(I)
2. Perform the fill as in the binary algorithm labeling the surrounding void as M+1. This means all voids are now either legitimate and can be filled or holes in-between labels.
3. Raster scan through the volume. If a new void is encountered, we must determine if it is fillable or an in-between one which will not be filled.
4. On encountering the void, record the last label seen and contour trace around it. If only that label is encountered during contour tracing, it is fillable. If another label is encountered, it is not fillable. 
5. During the contour trace, mark the trace using an integer not already used, such as M+2. If that label is encountered in the future, you'll know what to fill between it and the next label encountered based on the fillable determination. This phase stops when either the twin of the first M+2 label is encountered or when futher contour tracing isn't possible (in the case of single voxel gaps).
6. (Inner Labels) If another label is encountered in the middle of a void, contour trace around it and mark the boundary with the same M+2 label that started the current fill.

### SciPy Comparison 

<p style="font-style: italics;" align="center">
<img height=384 src="https://raw.githubusercontent.com/seung-lab/fill_voids/master/comparison.png" alt="Filling five labels using SciPy binary_fill_holes vs fill_voids from a 512x512x512 densely labeled connectomics segmentation. (black) fill_voids 0.1 (blue) scipy 1.3.3" /><br>
Fig. 1: Filling five labels using SciPy binary_fill_holes vs fill_voids from a 512x512x512 densely labeled connectomics segmentation. (black) fill_voids 0.1 (blue) scipy 1.3.3
</p>

In this test, fill_voids is significantly faster than scipy at lower memory. 