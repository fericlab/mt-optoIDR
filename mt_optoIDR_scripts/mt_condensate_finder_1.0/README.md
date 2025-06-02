![Static Badge](https://img.shields.io/badge/mt%20condensate%20finder%20-%201.0-red)
# Google Colab Notebook
## This Colab notebook can be used online for single-file use cases. 
Colab notebook: https://colab.research.google.com/drive/1ri0wl73W9sENBOHmWwQ6B1Lrs1YzBbGV?usp=sharing
# How to use mt_condensate_finder_1.0 to segment mitochondrial droplets

## Requiremnts 
## Install necessary dependencies
    pip install opencv-python
    pip install czifile
    pip install trackpy

## Running 
Given an input file, size data, and centroid estimates for an inputted threshold are returned. The program previews the image with the thresholded puncta and allows for modifications of the threshold on the fly.   
## Inputs:   
    # folder        - [directory] complete directory of the input folder containing files to be run
    # p_to_mu       - [float] pixel to micron conversion ratio
## Optional Inputs:
    # All arrays can have either one value for each channel, or one element which will be applied to all channels
    # timesteps     - [array of ints] an array containing the first timestep and the final timestep to be analyzed delimited by commas
    # zstacks       - [array of ints] an array containing the first zstack and the final zstack to be analyzed delimited by commas
    # channels      - [array of ints] an array the channel indices of the czifile to be analyzed
    # intensities   - [array of ints] an array of intensities per each channel
    # backints      - [array of floats] estimate background noise outside of the mitochondrial network
    # masksizes     - [array of ints] approximate radius (in pixels) of the puncta. Overestimate
    # windowsizes   - [array of ints] window size from the center to the size taken around each puncta
    # crop          - [xmin, xmax, ymin, ymax] array that crops the original image to the arrays values
    # plot          - 'Suppress' if you want to suppress the display of the segmentation for each image trace
## Outputs:
    # totallist     - returns a dictionary of dictionaries that contain the puncta found in each image trace
    # arraylistt    - returns a dictionary of dictionaries that contain the processed image of each puncta in each image trace
    # imagelistt    - returns a dictionary of dictionaries that contain the raw image of each puncta in each image trace
    # dframe        - returns a dataframe containing parameters for each puncta, only if the puncta were fit nicely 
    # rawdframe     - returns a dataframe containing parameters for each puncta
    
