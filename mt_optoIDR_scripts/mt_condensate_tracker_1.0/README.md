![Static Badge](https://img.shields.io/badge/mt%20condensate%20tracker-%201.0-green)
# Google Colab Notebook
## This Colab notebook can be used online for single-file use cases. 
Colab notebook: https://drive.google.com/file/d/1eyep1YvroG-o79zs1A5RXpEsLlXviWPn/view?usp=sharing 

# How to use mt-condensate-tracker for tracking droplets in mitochondria

## Requirements 
## Install necessary dependencies 
    pip install czifile
    pip install trackpy

## Guide
    1. Download the mt_droplet_tracker.py and tracker_runing.py   
    2. Keep your .czi or .tiff files in the same folder   
    3. In tracker_runing.py, change the file path to the same directory 
    4. Adjust other tracking parameters based on your imaging data.  

## trackdroplets Parameters 
        
    Parameters
    ----------
    file :              -[file or file like object] czi file or a TIF file
    filetype :          -[str] file type 'CZI' or 'TIF' 
    start_frame :       -[int] first frame number, first frame is 0
    end_frame :         -[int] end frame number
    channel :           -[int] channel number, first channel is 0        
    p_threshold :       -[float] precentile threshold to isolate puncta by trackpy. From 0-100%. Higher values are more specific
    search_range:       -[float or tuple] search range for linking   
    memory :            -[int] memory is the frames that can be skipped. 0 = no skipping     
    imaging_interval :  -[int] imaging intervals in secounds     
    annotate :          -[bool] True if you want to see the droplets identified      
    filter_threshold :  -[int] only considers tracks that have been detected over a certain amount of frames   
    correct_drift :     -[bool] True if drift correct on
        

    Returns
    -------
    DataFrame
        Data frame with linked puncta
## individualTracks Parameters 
    Parameters
    ----------
    df_linked :     -[DataFrame] linked dataframe from trackdroplets function
    
            
    Returns
    -------
    tracks_list :   -[DataFrame] A list of individual tracks 

## RotateTracks Parameters
    Parameters
    ----------
    tracks :         -[list] isolated track list from individualTracks       
    plotting :       -[bool] True if you want to see the rotated track
               
    Returns
    -------
    xAlined_tracks :      -[DataFrame] dataframe with xAlined coordinates
## lagTData Parameters  
    Parameters
    ----------
    tracks_list :     -[list] list with individual tracks rotated (xAlined)
    xAlined :         -[bool] if True, it will use xAlined x and y coordinates. if False, it will use the original coordinates
    mpp :             -[float] microns per pixel

    Returns
    -------
    lagTData :        -[DataFrame] conatin dx, dy, and corresponsing lagtime

## singleVanHovePlots Parameters
    Parameters
    ----------
    lagTData :      -[DataFrame] data frame from lagTData function
    trackID :       -[str] ID of the track ex. track '00' (00 = particle number + file number)
    lagT :          -[int] possible lag time in seconds
    nbins :         -[int] number of bins 
   
    Returns
    -------
    None. Only for visulaization purposes. 

## dataPool_vanHoverPlot Parameters
    Parameters
    ----------
    lagdata :       -[DataFrame] data frame from lagTData function
    lagT :          -[int] possible lag time in seconds
    dimention :     -[str] 'x' , 'y', or 'xy'
    nbins :         -[int] number of bins 
       

    Returns
    -------
    d_total :       -[list] raw distance values  
    bin_centers :   -[array] bin centers of the histogram (density = True, with normalized probability in y)
    counts :        -[array] normalized probability since density = True.


## MSD_cal Parameters 
    Parameters
    ----------
    lagttracks :        -[list] list of calculated lag ts from lagtData function 
    msd_dimention :     -[str] 'x' , 'y', or 'xy'
    ref_graph :         -[bool] if True plots slope = 1 graph
        
    
    Returns
    -------
    MSD_data 



# Output
    trackingData.csv :  contain calculated dx and dy data and appropriate lag times based on other assigned parameters (ex. xAlined is True or False)   
    1st_track_data_xAlined.csv :  an example data set for the first isolated track   

All the other data can be extracted from saved lists.
