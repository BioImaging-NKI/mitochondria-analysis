# mitochondria-analysis
### Morphological analysis of the mitochondria network in cultured cells

A [Fiji]((https://imagej.net/software/fiji/)) macro to automatically segment and measure the mitochondria network morphological properties and intensity in 3D.

Please note that before running this macro, cell segmentation should already be performed, and saved as 2D ImageJ ROI files, named `[File name without extension]_ROIs.zip` in the same folder as the input images.

### Workflow
- The intensity values for each cell are normalized (between 25th and 98th percentile) to account for differences in illumination over the field of view.
- Mitochondria are segmented using the Voronoi Threshold labeler method, as implemented in the [BioVoxxel 3D Toolbox](https://biovoxxel.github.io/bv3dbox/) (currently with hardcoded 0.5 pixels Gaussian filter, Moments threshold, and no further separation of connected labels).
- Neighboring labels with a certain maximum distance from each other are merged using [CLIJ2 libraries](https://clij.github.io/clij2/), after which the number of separate mitochondria structures (labels) in the cell were counted.
- For quantifying other morphological features the label image is binarized and skeletonized, and the skeleton analyzed using the [Analyze Skeleton plugin](https://imagej.net/plugins/analyze-skeleton/).

Images are saved with segmented (merged) mitochondria network structures, and the skeleton as overlays, as well as a table with morphological and intensity features.

![Animated](https://github.com/user-attachments/assets/844b7a37-2ab9-4675-93e0-762f5bed5ea7)

![image](https://github.com/user-attachments/assets/4f42e89b-0f56-4790-a694-e99d1340ff14)
