import czifile as czifile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import trackpy as tp
import os
import cv2
from scipy.optimize import curve_fit


def droplet_density(file, mpp, plotting, x_bleach, y_bleach, see_droplets,see_masked_droplets, p_thre_for_network, p_thre_for_droplet_back, p_thre_to_locate_droplets):
    #read the czi file 
    frames = czifile.imread(file)
    
    data_frame = pd.DataFrame(columns=['frame number', 'number of droplets', 'network area', 'average network int', 'average droplet int'])
    whole_dataf = pd.DataFrame()
    
    individual_droplet_data = []
    
    for i in range(len(frames)):
        
        if i<3: # skip first 4 frames
            continue
        
        #IN TRACKPY --> IDENTIFY A DROPLET
        image = frames[i,1,0,:,:,0]
        
        
        ##select a defined area from the activation spot and crop this part
        x=x_bleach
        y=y_bleach
        r=200 #half width of the extraction window
        image = image[int(y)-int(r):int(y)+int(r), int(x)-int(r):int(x)+int(r)] 
        
        #----------------------------------------------------------------------- part 1
        #identify puncta 
        threshold = tp.find.percentile_threshold(image, p_thre_to_locate_droplets)
        
        blobs_thresholded = tp.locate(image, (5, 13), threshold=threshold)
        #locate blob like features (droplets here)
        
        if len(blobs_thresholded) == 0: #if droplets were not detected skip this frame
            continue
       
       
            
        
        #----------------------------------------------------------------------- part 2
        #get the mitochondrial area
        net_threshold = tp.find.percentile_threshold(image, p_thre_for_network) #find a threshold using adaptive threshold 
        network = image > net_threshold
        network = network * image #not binary 
        
        #non - zero pixels 
        network_non_zero = network[network>0]
        net_area = mpp * mpp * len(network_non_zero) #squared micrometer 
        
        #--------------------------------------------------------------------- mask the droplets to get the network intensity
        droplets_mask = np.zeros(network.shape[:2], dtype=np.uint8)
        
        for _, row in blobs_thresholded.iterrows(): #apply the mask to the image
            x, y = int(row['x']), int(row['y'])
            cv2.circle(droplets_mask, (x,y), radius=4, color=255, thickness=-1)
        
        #get the average intensity of the masked droplets
        droplet_intensity = cv2.mean(network, mask=droplets_mask)[0]
        
        masked_network = network.copy()
        masked_network[droplets_mask==255]=0 #make 0 where the mask is applied 
        
        #ave network int
        network_masked_non_zero = masked_network[masked_network>0]
        ave_net_int = np.mean(network_masked_non_zero) #average network intensity
        #network area 
        
        
        
        #----------------------------------------------------------------------- vizualization 
        if plotting == True and i==4:
            plt.imshow(network)
            plt.title(f'Backg substraction threshold = {p_thre_for_network}')
            plt.show()
        if plotting ==True:
            #annotate located droplets
            tp.annotate(blobs_thresholded, image,  plot_style = {'markersize':6})
            plt.show()
        if see_masked_droplets == True:
            plt.imshow(masked_network)
            plt.title('droplets masked-no droplets')
            plt.show()
        
        
        #-------------------------------------------------------------------------- part 3 - analyze individual droplets if needs    
        ##Isolate each droplet and calculate the inside concentration and background concentration
        #isolation of droplets
        if see_droplets == True:
            for j in range(len(blobs_thresholded)):
                
                drop_x_loc = blobs_thresholded.x[j]
                drop_y_loc = blobs_thresholded.y[j]
                dw=8 #droplet_window
                droplet =  image[int(drop_y_loc)-int(dw): int(drop_y_loc)+int(dw), int(drop_x_loc)-int(dw): int(drop_x_loc)+int(dw)]
                
                #check whether the droplet is cropped properly (when close to the boundary cropping does not work properly)
                if droplet.shape[0] == 0 or droplet.shape[1] == 0:
                    print(f"Invalid crop for droplet {j}. Skipping.")
                    continue
                
                #background substraction
                net_threshold = tp.find.percentile_threshold(image, p_thre_for_droplet_back)
                droplet_thresholded = droplet > net_threshold #binatry 
                droplet = droplet_thresholded * droplet
                
                
                
                #DEFINE the inside and outside
                # Create a mask
                dop_center = (dw, dw)
                m_radius =4 #mask radius
                mask = np.zeros_like(droplet, dtype=np.uint8)
                cv2.circle(mask, dop_center, m_radius, 255, -1)  # circle as mask
                
                # Darken the masked region (set to black), keep the outside region unchanged
                masked_droplet = droplet.copy()
                masked_droplet[mask == 255] = 0  # Set masked area to dark (0 intensity)
               
                
                # Calculate average intensity inside and outside
                inside_intensity = cv2.mean(droplet, mask=mask)[0]
                #print(inside_intensity)
                
                #for the outside intensity consider only non zero pixels 
                outside_non_zero_px = masked_droplet[masked_droplet>0]
                outside_intensity = np.mean(outside_non_zero_px)
                #print(outside_intensity)
                
                if see_droplets==True:
                    #visualize 
                    fig, ax1 = plt.subplots(1,2)
                    ax1[0].imshow(droplet)
                    ax1[1].imshow(masked_droplet, vmin=droplet.min(), vmax=droplet.max())
                    ax1[0].set_title(f'droplet {j}')
                    ax1[1].set_title('masked droplet')
                    plt.show()
                
                #store these inside a dataframe for further use 
                individual_droplet_data.append({'droplet':j, 
                                                'x': drop_x_loc,
                                                'y': drop_y_loc,
                                                'inside intensity': inside_intensity,
                                                'outside intensity': outside_intensity,
                                                'frame': i})
        
        #---------------------------------------------------------------------------- export  data
        #record the number of droplets and other data to the data frame
        num_droplets = len(blobs_thresholded)
        data_frame.loc[i]=[i, num_droplets, net_area, ave_net_int, droplet_intensity]
        
        #write the whole dataframe
        blobs_thresholded['frame number']=i #add a column with the frame number
        whole_dataf=pd.concat([whole_dataf,blobs_thresholded], axis=0, ignore_index=True)
        
        #convert to a dataframe
        individual_drop_data_df = pd.DataFrame(individual_droplet_data)
    
    return data_frame, whole_dataf, individual_drop_data_df



###################### RUN the SCRIPT ###########################################





file=r'C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/S7/Activation/20250207_DDX28-mCh-Cry2olig_LocalAct_AL_CellNo4-01_Airyscan_processed.czi'

data, whole_df, individual_drop_data = droplet_density(file=file,
                             mpp=0.049,
                             plotting=True,
                             see_droplets=False, #if ture it isolate single droplets and get the inside and outside average intensity 
                             see_masked_droplets=False, 
                             p_thre_for_network=80, 
                             p_thre_for_droplet_back=70,
                             x_bleach= 327,
                             y_bleach= 360,
                             p_thre_to_locate_droplets=94)





############################################## Pllotting and observe data #########################################

#plot design parameters
font_properties = {
    'family': 'Arial',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42
tick_font_size = 14





#calculate droplet density
data['droplet density'] = data['number of droplets']/data['network area']
plt.plot(data['frame number'], data['droplet density'], alpha=0.5, color='gray')

plt.xlabel('frame number')
plt.ylabel('Droplet Density (droplets/um^2)')
#plt.legend()
#plt.ylim(0,1)
plt.show()


#------------------------------------------------------------------------------------------------ Fit to a model to get nucleation rate , J
########fit to the nucleation model refer the paper Shimobayashi, S.F., Ronceray, P., Sanders, D.W. et al. Nucleation landscape of biomolecular condensates. Nature 599, 503–506 (2021). https://doi.org/10.1038/s41586-021-03905-5

time = data['frame number'] *5 #5 = time between two frames is 5 secounds 
rho_data = data['droplet density']

time = time -20 #as activation stats at 5th frame 

#define the model
def nucleation_equation(t, rho_0, t_0, tau):
    return rho_0 * (1 - np.exp(-(t - t_0) / tau))


params, covariance = curve_fit(nucleation_equation, time, rho_data, p0=[0.01,20,2]) #can provide initial guesses too p0=[.., ..., ..., ]
rho_0_fit, t_0_fit, tau_fit = params

# Calculate fitted density and nucleation rate
rho_fit = nucleation_equation(time, rho_0_fit, t_0_fit, tau_fit)

## Nucleation rate J = dRho/dt = rho_0/tau * e**-(t-t0)/tau
#at t=0, J = rho_0/tau * e**t0/t
J = rho_0_fit/tau_fit * np.exp(t_0_fit/tau_fit)

#plot
plt.figure(figsize=(5,4))
plt.scatter(time, rho_data, color='gray', label='Data', alpha=0.6)
plt.plot(time, rho_fit, label='fitr', color='black')
plt.xlabel('Time(s)', **font_properties)
plt.ylabel('ρ (droplets/µm2)', **font_properties)
plt.annotate(f'Nucleation rate (J(t0)) = {J:.5f}', xy=(0.1,0.3), xycoords='axes fraction')

plt.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
plt.tick_params(axis='y', labelsize=tick_font_size)

plt.ylim(0,0.3)
plt.xlim(0, 495)

plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('Nucleation rate.pdf')
plt.show()

#### Export data ##

Nucleation_rate_Data = pd.DataFrame({'time(s)': time,
                                     'av_drop_density': rho_data,
                                     'fit': rho_fit
                                        })


Nucleation_rate_Data.to_csv('Nucleation_rates_data.csv')









############################################################################################

#save data
file_name = os.path.basename(file)
data.to_csv(file_name + '.csv', index=False)



############### from individual droplet data
#plt.plot(individual_drop_data.frame.unique(), individual_drop_data.groupby('frame')['inside intensity'].mean(), label='inside av int')

#plt.plot(individual_drop_data.frame.unique(), individual_drop_data.groupby('frame')['outside intensity'].mean(), label='outside av int')
#plt.legend()


plt.plot(data['frame number'], data['average droplet int'], label='ave droplet int')
plt.plot(data['frame number'], data['average network int'], label='ave network int')
plt.xlabel('frame number')
plt.ylabel('network average intensity (a.u)')
plt.legend()
plt.show()














