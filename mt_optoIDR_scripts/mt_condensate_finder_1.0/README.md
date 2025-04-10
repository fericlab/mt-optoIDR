![Static Badge](https://img.shields.io/badge/mt_Droplets-1.0-red)
# How to use mt_condensate_finder_1.0 to segment mitochondrial droplets

## Requiremnts 
## Install necessary dependencies
    pip install opencv-python
    pip install czifile
    pip install trackpy

## Running 
 Given an input file, size data and centroid estimates for an inputted threshold are returned. The program previews the image with the thresholded punctae and allows for modifications of the threshold on the fly.   
## Inputs:   
    # imagefile     - [directory] complete file directory of the input image. takes .czi files
    # channel       - [int] channel index of .czi file
    # p_to_mu       - [float] pixel to micron conversion ratio
    # I_thr         - [float] initial threshold value. increase to make more specific and 
    #                 decrease to gather more puncta
    # back_i        - [float] estimate of the cytoplasmic background signal
    # p_size        - [int] approximate radius (in pixels) of the puncta. Overestimate
    # wndw          - [int] window size from the center to the size taken around each puncta, in pixels
    # interactive   - [1/0] 1 means user output is not suppressed and 0 suppresses all
##                 outputs
    # Outputs:
    # centroids     - returns the x and y coordinates of the centroid
    # pixelAxis     - returns the major and minor axis lengths in pixels
    # punctaAxis    - returns the major and minor axis lengths in microns
    # aspectRatio   - returns the aspect ratio of the punctae
    # threshold     - returns the final threshold value used for the image
    
