### Tracking droplets Sanjaya AG 20241106 
import czifile as czifile
import matplotlib.pyplot as plt
import trackpy as tp
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import linregress
import tifffile

def tackdroplets(file, filetype, start_frame, end_frame, channel, p_threshold, search_range, memory, imaging_interval, annotate, filter_threshold, correct_drift):
    
    
    #THIS FUNCTION TRACK DROPLETS AND LINK THEM THEN CORRECT FOR ANY DRIFTS 
    
    #import the the file 
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
    
    print(len(frames))
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
                tp.annotate(f, frames[i], ax=ax, plot_style={'markersize': 10})
                ax.set_title(f'Frame {i}')
                plt.show()
    
    #add the file name #with the path to the dataframe
    tp_df['file path']=file
    
    ### linking tracks######
    linked = tp.link(tp_df, search_range=search_range, memory=memory)
    fig1, (ax1, ax2) = plt.subplots(1,2)
    tp.plot_traj(linked, superimpose=frames[0], colorby='particle', ax=ax1, label=False)
    ax1.set_title('All tracked')
    
    
    ##filtering -OPTIONAL
    linked = tp.filter_stubs(linked, threshold=filter_threshold) #threshold = minimun # of frames to survive 
    
    
    #check whether linked and filtered steps resulted any trajectories if yes go ..
    if len(linked)>0: #plot trajectories if only it exist
        tp.plot_traj(linked, superimpose=frames[0], ax=ax2)
        ax2.set_title(f'Filtered: Filter={filter_threshold} ')
        plt.tight_layout()
        plt.show()
        
        #correct any drifts if needed
        if correct_drift==True:
            #compute drift
            drift = tp.compute_drift(linked, smoothing=4)
            fig, ax3=plt.subplots()
            drift.plot(ax=ax3)
            ax3.set_title('Drift')
            
            #substract drift
            drift_corrected = tp.subtract_drift(linked, drift)
            drift_corrected = drift_corrected.reset_index(drop=True) #delete duplicated indexes 
            
            fig2, (ax4,ax5) = plt.subplots(1,2)
            tp.plot_traj(linked, superimpose=frames[0], ax=ax4)
            tp.plot_traj(drift_corrected, superimpose=frames[0], ax=ax5)
            
            ax4.set_title('Filtered')
            ax5.set_title('Drift corrected')
            plt.show()
            
            return drift_corrected
        
        else:
            return linked
        
    else:
        print('--------Zero tracjectories after filtering step-------')
   
    
    
   
    
    


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
            plt.plot(coordinates[:,0], coordinates[:,1], label='original')
            plt.plot(rotated_xy[:,0], rotated_xy[:,1], label=f'rotated by {angle_deg}')
            plt.title('xAlined track')
            plt.legend()
            plt.gca().invert_yaxis()
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
        file_numbers = np.array(i['file number'])
        
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
                file_number = file_numbers[j]
                
                new_row = pd.DataFrame({'dy':[dy], 'dx': [dx], 'dt':[dt], 'particle':[particle], 'file number':[file_number]})
                lagTDataFrame = pd.concat([lagTDataFrame, new_row])              
                
        lagTData.append(lagTDataFrame)
    
    print('DONE')
    return lagTData
                
                 


################################################### van Hover Plots ###############################

def singleVanHovePlots(lagTData, trackID, lagT, nbins):
    #extract data from the lagTData funciton 
    data = lagTData[trackID]
    #isolate data for a given lagt
    data_for_lagt = data[data['dt']==lagT]
    
    #get dx data 
    dx=np.array(data_for_lagt['dx'])
    
    #this is in terms of histogram
    
    plt.hist(dx, bins=nbins, density = True) #density=True --> normalized probability, area under the cruve is =1
    plt.yscale('log')
    plt.xlabel('Displacement of x (\u03BCm)')
    plt.ylabel('log(probability)')
    plt.title(f'Hist: van Hove plot for lagt = {lagT}, trackID = {trackID}, nbins={nbins}')
    plt.show()
    
    #with data
    # Calculate histogram data with density=True to normalize
    counts, bin_edges = np.histogram(dx, bins=nbins, density=True)

    # Calculate bin centers for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 #(lastvalue-firstvalue)/2
    
    
    #try to fit a Gaussian distribution to data 
    mean, std = norm.fit(dx) #Naormal fit (gaussian)--> get mean and std
    #plot the fucntion 
    x_values = np.linspace(min(dx), max(dx), 1000)
    gaussian_curve = norm.pdf(x_values, mean, std)
    
    fig, ax = plt.subplots()
    #counts=np.log10(counts) #make each value log
    ax.plot(bin_centers, counts, '-', marker='o', color='gray', label='Data')  # 'o' creates data points only
    ax.plot(x_values, gaussian_curve, label='Gaussian Fit', color='red')
    ax.set_yscale('log')  
    ax.set_xlabel('Displacement of x (\u03BCm)')
    ax.set_ylabel('log (Probability)')
    ax.set_title(f'van Hove plot for lagt = {lagT}, trackID = {trackID}, nbins={nbins}')
    ax.legend()
    plt.show()

    
    return None
    

    
############################## Polling all the lag T data together ##############
def dataPool_vanHoverPlot(lagdata, lagT, dimention, nbins):
    
    # this function pool all the data from the dimention you want
    
    d_total = []
    for i in lagdata:
        #isolate data for a given lagt
        data_for_lagt = i[i['dt']==lagT]
        
        if dimention=='x':
            d_total.extend(np.array(data_for_lagt['dx']))
        elif dimention=='y':
            d_total.extend(np.array(data_for_lagt['dy']))
        elif dimention=='xy':
            d_total.extend(np.array(data_for_lagt['dy']))
            d_total.extend(np.array(data_for_lagt['dx']))
        else:
            print('Please specify the dimention x, y or xy')
    #val Hover plot
    #get dx data 
    d=d_total
    
    #this is in terms of histogram
    
    plt.hist(d, bins=nbins, density = True) #density=True --> normalized probability, area under the cruve is =1
    plt.yscale('log')
    plt.xlabel('Displacement of x (\u03BCm)')
    plt.ylabel('log(probability)')
    N=len(d)
    plt.title(f'Hist: Pooled {dimention}: van Hove plot for lagt = {lagT}, nbins={nbins}, N={N}')
    plt.show()
    
    #with data
    # Calculate histogram data with density=True to normalize
    counts, bin_edges = np.histogram(d, bins=nbins, density=True)

    # Calculate bin centers for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 #(lastvalue-firstvalue)/2
    
    
    #try to fit a Gaussian distribution to data 
    mean, std = norm.fit(d) #Naormal fit (gaussian)--> get mean and std
    #plot the fucntion 
    x_values = np.linspace(min(d), max(d), 1000)
    gaussian_curve = norm.pdf(x_values, mean, std)
    
    fig, ax = plt.subplots()
    #counts=np.log10(counts) #make each value log
    ax.plot(bin_centers, counts, '-', marker='o', color='gray', label='Data')  # 'o' creates data points only
    ax.plot(x_values, gaussian_curve, label='Gaussian Fit', color='red')
    ax.set_yscale('log')  
    ax.set_xlabel(f'Displacement of {dimention} (\u03BCm)')
    ax.set_ylabel('log (Probability)')
    N=len(d)
    ax.set_title(f' Pooled {dimention} van Hove plot for lagt = {lagT}, nbins={nbins}, N={N}')
    ax.legend()
    plt.show()
    

    return d_total, bin_centers, counts
        
        


########### Calculating MSD values #############################

def MSD_cal(lagttracks, msd_dimention, ref_graph): 
    #get all the lag t data to a one dataframe
    combined_lagtdata = pd.concat(lagttracks, ignore_index=True)
    
    MSD_data = (
        combined_lagtdata.groupby('dt')[['dy', 'dx']] #get each unique dt, and extract dx and dy
        .apply(lambda group: (group**2).mean()) #apply square elemet wise and get the mean
        .reset_index() #put dt to a column
    )
    #Modify column names
    MSD_data.columns = ['dt', 'mean_squared_dy', 'mean_squared_dx']
    
    #MSD(τ) = <Δr(τ ))2> = <(ΔX(τ ))2> + <(ΔY (τ ))2> + <(ΔZ(τ ))2>
    if msd_dimention == 'x':
        MSD_data['msd']=MSD_data['mean_squared_dx'] 
    elif msd_dimention == 'y':
        MSD_data['msd']=MSD_data['mean_squared_dy']
    elif msd_dimention == 'xy':
        MSD_data['msd']=MSD_data['mean_squared_dy'] + MSD_data['mean_squared_dx'] 
    else:
        print('specify the dimention x,y or xy')
    
    #plotting 
    plt.plot(MSD_data['dt'], MSD_data['msd'], label='MDS')
    plt.xlabel('lag time (s)')
    plt.ylabel('MSD (\u03bcm$^2$)')
    plt.title(f'log-log plot of MSD in {msd_dimention}')
    
   
    
    #REFERENCE plot slope=1
    ref_x = np.linspace(MSD_data.dt.min(), MSD_data.dt.max(), 1000)
    y=ref_x
    
    if ref_graph==True:
        plt.plot(ref_x, y, label="Slope=1")
        
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.show()
    
    
    return MSD_data



    

    
    
