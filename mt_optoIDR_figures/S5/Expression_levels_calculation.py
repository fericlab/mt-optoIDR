import czifile as czifile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import trackpy as tp
import os

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




def expression_levels (file, p_thre_for_network, mpp):
    #read the czi file 
    frames = czifile.imread(file)
    
    data_frame = pd.DataFrame(columns=['frame number', 'network area', 'average network int'])
    
    for i in range(len(frames)):
        
        if i>3: # stop at frame 4, as activation happens after 
            continue
        
        #IN TRACKPY --> IDENTIFY A DROPLET
        image = frames[i,1,0,:,:,0]
        plt.imshow(image)
        plt.title('Original image')
        plt.show()
        
        #----------------------------------------------------------------------- part 2
        #get the mitochondrial area
        net_threshold = tp.find.percentile_threshold(image, p_thre_for_network) #find a threshold using adaptive threshold 
        network = image > net_threshold
        network = network * image #not binary 
        plt.imshow(network)
        plt.title(f'Thresholded at {net_threshold}')
        plt.show()
        
        #non - zero pixels 
        network_non_zero = network[network>0]
        net_area = mpp * mpp * len(network_non_zero) #squared micrometer 
        
        ave_net_int = np.mean(network_non_zero) #average network intensity
        data_frame.loc[i] = [i,net_area, ave_net_int]
       
        
       
    return data_frame


######################################## Running this ###################################

path_Drp1=r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\Data_Analysis\Expression_levels_comparision\data\Drp1k38a'
Path_WT = r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\Data_Analysis\Expression_levels_comparision\data\WT'

paths = [Path_WT, path_Drp1]
threds =[96, 60] #thresholds for WT and Drp1 bulbus 

for j in range (len(paths)):
    path = paths[j]
    directory = os.listdir(path)
    files=[]
    for i in range(len(directory)):
        if directory[i].endswith('.czi'): # removes a hidden file 
            files.append(path + '/' + directory[i])
        else:
            pass
            
    
    All_data = []
    
    for file_no in range(len(files)):
        file_name=files[file_no]
        data = expression_levels(file=file_name,
                          p_thre_for_network=threds[j], 
                          mpp=0.048)
        data['file number']=file_no #add a new column with the file number
        All_data.append(data)
    
    
    
    ########## Getting the average of each file###########
    All_data_combined = pd.concat(All_data)
    
    Average_network_ints = All_data_combined.groupby('file number')['average network int'].mean()
    
    print('Average network intensity:', Average_network_ints.mean())
    
    #save
    All_data_combined.to_csv(f'intensity_data_{j}.csv', index=False)
    plt.hist(Average_network_ints, bins=8)
    plt.show()
    

############## get the data again -- so easy to see ############

wt = pd.read_csv('intensity_data_0.csv')
drp1 = pd.read_csv('intensity_data_1.csv')

wt_means = wt.groupby('file number')['average network int'].mean()
wt_mean = wt_means.mean()
print('Wild type network average intensity:', wt_mean)


drp1_means = drp1.groupby('file number')['average network int'].mean()
drp1_mean = drp1_means.mean()
print('Drp1K38A network average intensity:', drp1_mean)


data_list = [wt_means, drp1_means]
############# violine plot #############
vp = plt.violinplot(dataset=data_list, showmeans=True, showmedians=False, showextrema=True)
plt.ylabel(" Average intensity (a.u)", **font_properties)
plt.xticks(ticks=[1, 2], labels=['WT', 'Drp1K38A',])
plt.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
plt.tick_params(axis='y', labelsize=tick_font_size)



colors = ['#564F8E', 'gray']
for i, body in enumerate(vp['bodies']):  
    body.set_facecolor(colors[i])  # Fill color
    body.set_edgecolor('black')    # Outline color
    body.set_alpha(0.7)  

vp['cmeans'].set_color('black')
vp['cbars'].set_color('black')  # Change whisker line color
vp['cbars'].set_linewidth(1)  # Adjust thickness

# Change the top and bottom caps (if needed)
vp['cmins'].set_color('black')  # Bottom cap
vp['cmaxes'].set_color('black')  # Top cap



plt.title('Average intensities of the mitochondria')
plt.annotate(f'wt_mean={wt_mean:.2f},\ndrp1_mean={drp1_mean:.2f}', xy=(0.1,0.5), xycoords='axes fraction')
plt.tight_layout()
plt.savefig('network_int_violin_plot.pdf')
plt.show()

