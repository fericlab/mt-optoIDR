import mt_condensate_tracker as tr
import os
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import czifile as czifile
import seaborn as sns

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

#STEP 1: TRACKING
###### For each file run track droplets and individual tracks ##########
path=r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper#1\Commun. Biol revisions\Revision1\FIGURES\S5\Same_cell_sizes\Size'


All_tracks =[]

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
                                 filetype='CZI', #file type CZI or TIF
                                 start_frame=0, 
                                 end_frame=20, 
                                 channel=1, # 0 is the first channel
                                 p_threshold=85, #to identify the puncta by trackpy 1-100%
                                 search_range=10, #for trackpy linking fuction
                                 memory=5, #memory is the frames that can be skipped 0=no skipping
                                 imaging_interval=2, #in sec
                                 annotate=False,
                                 filter_threshold=10, #filter if it can not be tracked at least this amount of frames
                                 correct_drift=False) #correct drifts !Warning be careful with the dift correction with dynamic systems!

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



#plot tracjctories with labels and mannualy pick tubular ones
#import the czi file 
cziimage = czifile.imread(files[0])
frame = cziimage[1,1,0,:,:,0]


#max particle number 
Max_particle_num = All_tracks[0][-1]['particle'].astype(int).max()


for t in range (Max_particle_num): #define the length , go to the last track and see the particle number
    for i in xAlined:
        particleNum=i['particle'].iloc[0]
        if t == particleNum:
            fig, ax = plt.subplots()
            ax.plot(i.x,i.y)
            ax.text(i['x'].iloc[-1], i['y'].iloc[-1], str(particleNum),
                   fontsize=6, ha='left', va='center', color='red')
            ax.set_title(particleNum)
    
            plt.imshow(frame, cmap='gray', alpha=0.8)
            plt.show()






############################################################################################ SELECTION ###########################

Tubular_particles = [5, 9, 11, 21, 40, 50, 129, 141]  
#Bulbous_particles = [ ]


################ Devide xAlined to two 1.tutular, 2.bulbous 
xAlined_tubular = []
xAlined_bulb = []

for i in xAlined:
    particleNum = i['particle'].iloc[0]   
    if particleNum in Tubular_particles:
        xAlined_tubular.append(i)
    else:
        #pass
        xAlined_bulb.append(i)
        

'''
#bulb slection

for j in xAlined:
    particleNum = j['particle'].iloc[0]   
    if particleNum in Bulbous_particles:
        xAlined_bulb.append(j)
    else:
        pass
        
'''





######################## NOw, other analysis for separate groups ######################
                    ############## TUBULER #################
############################################################################################## 

# 1. SIZE (in terms of Radius of gyration)
#isolate average radius of gyration in x and y 
Diameter_list_TB = []
for i in xAlined_tubular:
    rg = (i['size_x'] + i['size_y'])/2      #effective radius - mean
    rg = rg.median()   
    Diameter = rg*2*0.049 #microns per pixel =0.049                        

    #rg = (i['size_x'].mean() * i['size_y'].mean())**0.5 #effective radius as geometric mean
    Diameter_list_TB.append(Diameter)
    
    

Diameter_list_BB = []
for i in xAlined_bulb:
    rg = (i['size_x'] + i['size_y'])/2      #effective radius - mean
    rg = rg.median()   
    Diameter = rg*2*0.049 #microns per pixel =0.049                        

    #rg = (i['size_x'].mean() * i['size_y'].mean())**0.5 #effective radius as geometric mean
    Diameter_list_BB.append(Diameter)



#plotting 

sizes_tubular = Diameter_list_TB
sizes_bulb = Diameter_list_BB

df = pd.DataFrame({
    "size": sizes_tubular + sizes_bulb,
    "group": ["Tubular"] * len(sizes_tubular) + ["Bulb"] * len(sizes_bulb)
})

plt.figure(figsize=(5,6))

# Boxplot
sns.boxplot(
    data=df,
    x="group",
    y="size",
    width=0.3,
    linewidth=1.5,
    fliersize=0,
    boxprops=dict(edgecolor='black', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='black', linewidth=1.2),
    whiskerprops=dict(color='black', linewidth=1.2),
    capprops=dict(color='black', linewidth=1.2)
)

# Overlay datapoints
sns.stripplot(
    data=df,
    x="group",
    y="size",
    dodge=False,
    jitter=True,
    alpha=0.6,
    size=6,
    color='black'   
)

plt.ylabel('Effective Diameter (μm)')
plt.xlabel("")
plt.tight_layout()
plt.savefig('Size.pdf', dpi=600)
plt.show()




######################################################################################################## van hover plots


R_tracks_blb = tr.lagTData(xAlined_bulb , xAlined=True, mpp=0.049)
R_tracks_tube = tr.lagTData(xAlined_tubular , xAlined=True, mpp=0.049)


#########################################################################################
############# one plot ###########

lgt=10 # prefered lag time

lagt=lgt
nbins=5
#WT
d, x_wt_x, y_wt_x = tr.dataPool_vanHoverPlot(R_tracks_tube, lagT=lgt, dimention='x', nbins=nbins)

d, x_wt_y, y_wt_y = tr.dataPool_vanHoverPlot(R_tracks_tube, lagT=lgt, dimention='y', nbins=nbins)

#Drp1k38a

d, x_drp1_x, y_drp1_x = tr.dataPool_vanHoverPlot(R_tracks_blb, lagT=lgt, dimention='x', nbins=nbins)

d, x_drp1_y, y_drp1_y = tr.dataPool_vanHoverPlot(R_tracks_blb, lagT=lgt, dimention='y', nbins=nbins)

plt.figure(figsize=(5, 4))
plt.plot(x_wt_x, y_wt_x, label='WT$_x$', alpha=1, marker='s', color='#7E287E', linewidth=1)
plt.plot(x_wt_y, y_wt_y, label='WT$_y$',  alpha=1, marker='s', color='gray', linewidth=1)

plt.plot(x_drp1_x, y_drp1_x, label='Drp1$^{K38A}_x$',  alpha=1, marker='o', color='#7E287E', linewidth=1)
plt.plot(x_drp1_y, y_drp1_y, label='Drp1$^{K38A}_y$',  alpha=1, marker='o', color = 'gray', linewidth=1)

plt.xlabel('Displacement (\u03BCm)')
plt.ylabel('Probability')
plt.yscale('log')
plt.title(f'Overlayed van Hove Plots lagt:{lagt}, nbins={nbins}')
#plt.xlim(0,1.8)
#plt.ylim(0.04,10)
plt.legend(frameon=False ) 
plt.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig('van_Hove_plot.pdf')
plt.show()

