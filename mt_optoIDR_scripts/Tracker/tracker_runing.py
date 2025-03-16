import Tracker as tr
import os
import pandas as pd
import itertools

#STEP 1: TRACKING
###### For each file run track droplets and individual tracks ##########
All_tracks =[]
path=r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\Data_Analysis\Dynamics\Data\Drp1k38a'
directory = os.listdir(path)
files=[]
for i in range(len(directory)):
    if directory[i].endswith('.czi'): # removes a hidden file 
        files.append(path + '/' + directory[i])
    else:
        pass
        

for i in range(len(files)):
    file=files[i]
    Linked_tracks = tr.tackdroplets(file=rf'{file}',
                                 start_frame=0, 
                                 end_frame=200, 
                                 channel=1, # 0 is the first channel
                                 p_threshold=95, #to identify the puncta by trackpy 1-100%
                                 search_range=10, #for trackpy linking fuction
                                 memory=5, #memory is the frames that can be skipped 0=no skipping
                                 imaging_interval=2, #in sec
                                 annotate=False,
                                 filter_threshold=10, #filter if it can not be tracked at least this amount of frames
                                 correct_drift=False) #correct drifts

    #if Linked retuts none skip
    if Linked_tracks is None :
        continue
    #add the file number
    Linked_tracks['file number']=i

    tracks_list = tr.individualTracks(df_linked=Linked_tracks) #Input data frame containing linked tracks 

    All_tracks.append(tracks_list)


#STEP 2: GET ALL THE TRACKS TO A SINGLE LSIT 
# Flatten using itertools.chain , combine all tracks to a single list
Combined_tracks = list(itertools.chain.from_iterable(All_tracks))   

#STEP 3: ALINE TO THE X AXIS
xAlined = tr.RotateTracks(tracks=Combined_tracks, plotting=False)

#STEP 4: CALCULATE, dx AND dy FOR ALL POSSIBLE lag t S
All_lagtdata = tr.lagTData(tracks_list=xAlined, xAlined=True, mpp=0.0488706004242)



#STEP 5: WRITE DATA INTO CSV
####write data to a csv###
combined = pd.concat(All_lagtdata, ignore_index=True)
combined.to_csv('data_Drp1k38a_xAlined.csv', index=False)





