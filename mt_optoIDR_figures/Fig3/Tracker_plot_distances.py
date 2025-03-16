### Tracking droplets Sanjaya AG 20241106 2.51 AM :) ###
import czifile as czifile
import matplotlib.pyplot as plt
import trackpy as tp
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import linregress
import tifffile
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as ticker

#plot design parameters
font_properties = {
    'family': 'Arial',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

def tackdroplets(file, start_frame, end_frame, channel, p_threshold, search_range, memory, imaging_interval, annotate, filter_threshold, filetype):
    '''
    Parameters
    ----------
    file : TYPE
        DESCRIPTION. file path to czi file
    start_frame : TYPE
        DESCRIPTION.
    end_frame : TYPE
        DESCRIPTION.
    channel : TYPE
        DESCRIPTION.
    p_threshold : TYPE
        DESCRIPTION. % Threshold for trackpy locate function
    search_range : TYPE
        DESCRIPTION.search range for linking function 
    memory : TYPE
        DESCRIPTION.How many frames can be skkiped to link
    imaging_interval : TYPE
        DESCRIPTION. Imaging intervals in sec
    annotate : TYPE
        DESCRIPTION.
    filter_threshold : TYPE
        DESCRIPTION. Remove if this numbers of frames can not be tracked 

    Returns
    -------
    linked : TYPE
        DESCRIPTION. Linking data frame

    '''
    if filetype == 'TIF':
        tif_image = tifffile.imread(file)
        frames = tif_image[start_frame:end_frame,channel,:,:]
        
    elif filetype == 'CZI':   
        #import the czi file 
        cziimage = czifile.imread(file)
        frames = cziimage[start_frame:end_frame,channel,0,:,:,0]
    else:
        print('please specify the file type : TIF or CZI')
        
    
    #create the dataframe to store data from locate function
    tp_df= pd.DataFrame()
    
    #loop the locate function through all the frames , do not use batch so you loose the precentile thresholding feature 
    for i in range (len(frames)):
        thr_value = tp.find.percentile_threshold(frames[i], p_threshold)
        f = tp.locate(frames[i], 
                      diameter=(15, 13), 
                      threshold=thr_value)
        if len(f)==0: #if droplets were not delected go to next frame
            continue
        else:
            #add the frame number to the DataFrame
            f['frame']=int(i) 
            #add time to the dataframe
            f['time'] = i*imaging_interval
            
            tp_df = pd.concat([tp_df,f])
            
            if annotate==True:
                fig, ax = plt.subplots()
                tp.annotate(f, frames[i], ax=ax, plot_style={'markersize': 5})
                ax.set_title(f'Frame {i}')
                plt.show()
    
    ### linking tracks######
    linked = tp.link(tp_df, search_range=search_range, memory=memory)
    fig1, (ax1, ax2) = plt.subplots(1,2)
    tp.plot_traj(linked, superimpose=frames[0], colorby='particle', ax=ax1, label=False)
    ax1.set_title('All tracked')
    
    ##filtering -OPTIONAL
    linked = tp.filter_stubs(linked, threshold=filter_threshold) #threshold = minimun # of frames to survive 
    tp.plot_traj(linked, superimpose=frames[0], ax=ax2)
    ax2.set_title(f'Filtered: Filter={filter_threshold} ')
    plt.tight_layout()
    #fig1.savefig('All_track_on_the_image.svg')
    plt.show()
    
    #only filtered image
    fig2, ax3 = plt.subplots()
    ax3.imshow(frames[13], cmap='inferno', vmax=3000, vmin=0)
    tp.plot_traj(linked, ax=ax3)
    
    
    
    # Apply axis formatter to convert pixels to micrometers
    pixel_to_um =0.049
    ax3.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{x * pixel_to_um:.1f}'))
    ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{y * pixel_to_um:.1f}'))
    
    

    ax3.tick_params(axis='x', labelsize=14)
    ax3.tick_params(axis='y', labelsize=14)
    num_ticks = 4
    ax3.xaxis.set_major_locator(MaxNLocator(num_ticks))
    ax3.yaxis.set_major_locator(MaxNLocator(num_ticks))
    
    
    ax3.set_xlabel('X position (µm)', **font_properties)
    ax3.set_ylabel('Y position (µm)', **font_properties)

    plt.tight_layout()
    fig2.savefig('track_on_the_image.pdf')
    plt.show()
    
    
    return linked


## steps --> Isolate individual trcks --> get a linear plot --> rotate parallel to the x axis --> calculate the x and y distance individually 
#STEP 1: Isolate individual tracks 
def individualTracks(df_linked):
    '''
    Parameters
    ----------
    df_linked : TYPE
        DESCRIPTION. Data frame retures from 'trackdroplet' function.

    Returns
    -------
    tracks_list : Dataframe
        DESCRIPTION. A list of individual tracks 

    '''
    droplets_list = df_linked['particle'].unique() #get unique IDs
    
    tracks_list =[] #store all the individual tracks
    for i in droplets_list:
        track_i = df_linked[df_linked['particle']==i] #filter out tracks 
        tracks_list.append(track_i)
    
    return tracks_list


# Rotate each track parallel to the x axis 
def RotateTracks(tracks, plotting):
    print('Runing xAline....')
    xAlined_tracks=[]
    
    for i in tracks:
        #Plot and fit a linear fit
        y=np.array(i['y'])
        x=np.array(i['x'])
        
        #step2.1
        # Fit a linear model 
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        r2=r_value**2
        #print(f'R^2 of the fitting = {r2}')
        
        #Get the angle between x and the filled line , angle is the arctan
        # Calculate the angle in radians
        angle_radians = np.arctan(slope)
        angle_deg = np.degrees(angle_radians)
        #print(angle_radians)

        # Generate the fitted line
        y_fit = slope * x + intercept
        
        if plotting==True:
            plt.scatter(x, y, label='Data points')
            plt.plot(x, y_fit, color='red', label='Fitted line')
            plt.legend()
            plt.title(f'angle to x = {angle_deg}')
            plt.gca().invert_yaxis() #invert is neccesory 
            plt.show()
        
        #step2.2
        ##########rotate coordinates############

        coordinates = i[['x', 'y']].to_numpy()
        
        #find the centroid
        centroid = np.mean(coordinates, axis=0)
        #Translate points to center of gravity
        Trans_coordinates = coordinates - centroid
        

        # Define the angle of rotation in degrees
        theta_rad = angle_radians

        # Define the rotation matrix
        rotation_matrix = np.array([
            [np.cos(theta_rad), -np.sin(theta_rad)],
            [np.sin(theta_rad),  np.cos(theta_rad)]
        ])

        # Rotate all points
        rotated_xy = np.dot(Trans_coordinates , rotation_matrix) #np.dot used to do matrix multiplication here not dot product
        # Translate the points back to their original position
        rotated_xy = rotated_xy + centroid
        
        
        if plotting==True:
            plt.figure(figsize=(4, 4))
            plt.plot(coordinates[:,0], coordinates[:,1], label='original')
            plt.plot(rotated_xy[:,0], rotated_xy[:,1], label=f'rotated by {angle_deg}')
            plt.legend()
            #plt.xlim(-10,60)
            #plt.ylim(0,60)
            plt.gca().invert_yaxis()
            #plt.savefig('rotated_tack.svg')
            plt.show()
       
        #put this to a the dataframe
        #create a copy of the dataframe so I do not modify the original
        i_copy = i.copy()
        
        i_copy.loc[:,'x_alined']=rotated_xy[:,0]
        i_copy.loc[:, 'y_alined']=rotated_xy[:,1]
        xAlined_tracks.append(i_copy)
    
    print('DONE')
    return xAlined_tracks



################ Extracting dx, dy and lagT data ############
#Need a dataframe with lag time, dx and dy
def lagTData(tracks_list, xAlined, mpp):
    #THIS FUNCTION GETS dX, dY FOR EACH POSSIBLE lag_time (dt). 
    #RETURN A LIST FOR EACH TRACK 
    print('Runing lagT...')
    lagTData = [] #list to store lagt data
    
    for i in tracks_list:
        lagTDataFrame=pd.DataFrame()
        
        
        if xAlined==True:
            y=np.array(i['y_alined']) #get to a np array for easy numerical manipulations 
            x=np.array(i['x_alined'])
        elif xAlined == False:
            y=np.array(i['y']) #get to a np array for easy numerical manipulations 
            x=np.array(i['x'])
        else:
            print('specify the aligment')
            
        t=0
        t_list = np.array(i['time'])
        particles = np.array(i['particle'])
        
        for k in range(len(i)):
            t+=1 #lagtime iteratable
            if t==len(i):
                break
           
            for j in range(len(i)): #j is the varaible to get each value
                if j==(len(i)-t):
                    break
                
                x1=x[j]
                x2=x[j+t]
                
                y1=y[j]
                y2=y[j+t]
                
                dx=(x2-x1)*mpp #convert it o microns 
                dy=(y2-y1)*mpp 
                
                dt=(t_list[j+t])-(t_list[j])
                
                particle=particles[j]
                
                new_row = pd.DataFrame({'dy':[dy], 'dx': [dx], 'dt':[dt], 'particle':[particle]})
                lagTDataFrame = pd.concat([lagTDataFrame, new_row])              
                
        lagTData.append(lagTDataFrame)
    
    print('DONE')
    return lagTData
                



    
#######################################################  Running for plotting ###############

linked = tackdroplets(file=r'C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/Fig3/Tracking_images/WT/20241004_cell12_crop.tif', 
             start_frame=0, 
             end_frame=45, 
             channel=1, 
             p_threshold=85, 
             search_range=20, 
             memory=2, 
             imaging_interval=2, 
             annotate=False, 
             filter_threshold=20,
             filetype='TIF')

tracks_list = individualTracks(df_linked=linked) #Get individual tracks


xAlined = RotateTracks(tracks=tracks_list, plotting=True) #Aline to X-axis if want


########################## PLOTTING #############################

############# Plotting x and y change on the same graph #######################
mpp = 0.049 #microns per pixels
normalize = False
#get the track
track = tracks_list[0] ###Give the tracklist as you want (XAlined or not)
#track = xAlined[0]

#get x y cordinates 
x_cord = np.array(track.x)
y_cord = np.array(track.y)

#get the difference
x_diff = np.diff(x_cord)
y_diff = np.diff(y_cord)

#convert from pixels to um 
x_diff = x_diff*mpp
y_diff = y_diff*mpp

#get absolute values
x_diff = np.abs(x_diff)
y_diff = np.abs(y_diff)
#get the time
time = np.array(track.time[:-1])



time_refd = time - time[0] #shift the data so first become zero


#max distance is it in x or y?
max_dis = x_diff.max()
print('max distance: ', max_dis)

#Normalize the movement
if normalize == True:
    x_diff = x_diff/max_dis
    y_diff = y_diff/max_dis

plt.figure(figsize=(5, 4))


#plot
strt =0
endd = 45
plt.plot(time_refd[strt:endd], x_diff[strt:endd], label='ΔX', color='#7E287E', alpha=1, marker='s')
plt.plot(time_refd[strt:endd], y_diff[strt:endd], label='ΔY', color='gray', alpha=1, marker='s')


##export data 
distance_data = pd.DataFrame({'time(s)': time_refd,
                              'x_diff': x_diff,
                              'y_diff': y_diff})


distance_data.to_csv('Distance_data_WT.csv')

# Draw a vertical line at the mean
x_mean= x_diff.mean()
y_mean = y_diff.mean()

#plt.axhline(x_mean, color='#7E287E', linestyle='--')
#plt.axhline(y_mean, color='gray', linestyle='--')


#plt.axhline(x_mean, color='#2CA02C', linestyle='--', label=f'Mean X = {x_mean:.2f}')
#plt.axhline(y_mean, color='gray', linestyle='--', label=f'Mean Y = {y_mean:.2f}')


plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('Time (s)', **font_properties)

if normalize == True:
    plt.ylabel('Normalized |Δd| (a.u)', **font_properties)
    plt.ylim(0,1.1)
#plt.title('Temporal change in x and y')
else:
    plt.ylabel('|Δd| (µm)', **font_properties)
    
#plt.xlim(-1, 45)
plt.ylim(0, 2.5)
plt.legend(frameon=False, fontsize=14)
plt.tight_layout()
plt.savefig('xy_change.pdf')



##export the tracking data
tracks_list[0].to_csv('track_data_WT.csv')


















