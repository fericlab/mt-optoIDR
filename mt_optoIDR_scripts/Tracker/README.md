# mt_drop_Tracker for tracking droplets in mitochondria 

# How to use Tracker.py
1. Download the Tracker.py and tracker_runing.py   
2. Keep your .czi to .tiff file in the same folder   
3. In tracker_runing.py, change the file path to the same directory 
4. Adjust other tracking parameters based on your imaging data.  

## trackdroplets function Parameters 
        
    Parameters
    ----------
    file :              -[file or file like object] czi file or a TIF file
    filetype :          -[str] file type 'CZI' or 'TIF' 
    start_frame :       -[int] first frame number, first frame is 0
    end_frame :         -[int] end frame number
    channel :           -[int] channel number, first channel is 0        
    p_threshold :       -[float] precentile threshold to isolate puncta by trackpy 0-100%, higher the value more specific
    search_range:       -[float or tuple] search range for linking   
    memory :            -[int] memory is the frames that can be skipped 0=no skipping     
    imaging_interval :  -[int] imaging intervals in secounds     
    annotate :          -[bool] True if you want to see the droplets identified      
    filter_threshold :  -[int] filter if it can not be tracked at least this amount of frames    
    correct_drift :     -[bool] True if drift correct on
        

    Returns
    -------
    DataFrame
        Data frame with linked puncta
## individualTracks parameters 
    Parameters
    ----------
    df_linked :     -[DataFrame] linked dataframe from trackdroplets function
    
            
    Returns
    -------
    tracks_list :   -[DataFrame] A list of individual tracks 

## RotateTracks parameters
    Parameters
    ----------
    tracks :         -[list] isolated track list from individualTracks       
    plotting :       -[bool] True if you want to see the rotated track
               
    Returns
    -------
    xAlined_tracks :      -[DataFrame] dataframe with xAlined coordinates
## lagTData parameters  
    Parameters
    ----------
    tracks_list :     -[list] list with individual tracks rotated (xAlined)
    xAlined :         -[bool] if True , it will use xAlined x and y coordinates, False, use original coordinates
    mpp :             -[float] microns per pixes

    Returns
    -------
    lagTData :        -[DataFrame] conatin dx, dy and corresponsing ladtime

## singleVanHovePlots parameters
    Parameters
    ----------
    lagTData :      -[DataFrame] data frame from lagTData function
    trackID :       -[str] ID of the track ex. track '00' (00 = particle number + file number)
    lagT :          -[int] possible lag time in sec
    nbins :         -[int] number of bins 
   
    Returns
    -------
    None. Only for visulaization propose. 

## dataPool_vanHoverPlot parameters
    Parameters
    ----------
    lagdata :       -[DataFrame] data frame from lagTData function
    lagT :          -[int] possible lag time in sec
    dimention :     -[str] 'x' , 'y', or 'xy'
    nbins :         -[int] number of bins 
       

    Returns
    -------
    d_total :       -[list] raw distance values  
    bin_centers :   -[array] bin centers of the histogram (density = True, with normalized probability in y)
    counts :        -[array] normalized probability since density = True.


## MSD_cal parameters 
    Parameters
    ----------
    lagttracks :        -[list] list of calculated lag ts from lagtData function 
    msd_dimention :     -[str] 'x' , 'y', or 'xy'
    ref_graph :         -[bool] if True plots slope = 1 graph
        
    
    Returns
    -------
    MSD_data 



# Output
    trackingData.csv :  contain calculated dx and dy data and appropriate lag times based on other parameters assigned (ex. xAlined is True or False)   
    1st_track_data_xAlined.csv :  an example data set for the first track isolated   

All the other data can be extracted from saved lists.
