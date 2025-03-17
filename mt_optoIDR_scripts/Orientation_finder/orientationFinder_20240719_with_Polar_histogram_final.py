import numpy as np
import czifile
import cv2
import droplet_v1_1_8July24 as drop
import skimage as ski
import sklearn as skl
import matplotlib.pyplot as plt
import os
from math import pi 
import math
import pandas as pd


def mitoOrientation(folder, tdf):
    # creates a list of all of the files in a given folder directory
    f = os.listdir(folder)
    files = []
    images = []
    for i in range(len(f)):
        if f[i].endswith('.czi'): # may be unnecesary, removes a hidden file #sanjaya
            files.append(folder + '/' + f[i])
        else:
            pass
    for j in range(len(files)):
        # reads the image and pulls out a frame and a channel
        filedata = czifile.imread(files[j])
        images.append(filedata)
        
    ftot = tdf['file number'].unique() 
    ctot = tdf['channel number'].unique() 
    frtot = tdf['frame number'].unique() 
    stot = tdf['stack number'].unique() 
    
    fskele = {}
    fim = {}
    for i in ftot:
        cskele = {}
        cim = {}
        for j in ctot:
            frskele = {}
            frim = {}
            for k in frtot:
                sskele = {}
                sim = {}
                for l in stot:
                    i2 = int(i.split()[1])
                    j2 = int(j.split()[1])
                    k2 = k #int(k.split()[1])
                    l2 = int(l.split()[1])
                    filedata = images[i2-1]
                    # creates the proper format for each frame of the image
                    if len(filedata.shape) == 6:
                        image = filedata[k2-1,j2-1,0,:,:,0]
                    elif len(filedata.shape) == 8:
                        image = filedata[k2-1,0,j2-1,0,l2-1,:,:,0]   
                    # linearizes the image to work with OpenCV's adaptive thresholding
                    minim = min(image.ravel())
                    maxim = max(image.ravel())
                    linear_image = np.multiply((image - minim), 255/(maxim-minim))
                    # creates a new image that masks the puncta based on the adaptive thresholding scheme
                    processed = cv2.GaussianBlur(linear_image, (3,3), 0)
                    adaptive = cv2.adaptiveThreshold(np.uint8(processed), 255, 
                                                     cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                                                     cv2.THRESH_BINARY, 11, -0.3)
                    # finds mean of this number
                    backImage = (1-adaptive/255) * image
                    backImageGuess = 2*sum(sum(backImage))/len(np.nonzero(backImage)[0][:])
                    back_i = backImageGuess
                    
                    for z in range(3):
                        adaptive = cv2.GaussianBlur(adaptive, (15,15), 0)
                        adaptive = cv2.adaptiveThreshold(np.uint8(adaptive), 255, 
                                                         cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                                                         cv2.THRESH_BINARY, 11, -0.3)
                        
                    adaptive = ski.morphology.skeletonize(adaptive)
                    adaptive = (adaptive*image>back_i)
                    
                    sskele[l] = adaptive
                    sim[l] = image
                frskele[k] = sskele
                frim[k] = sim
            cskele[j] = frskele
            cim[j] = frim
        fskele[i] = cskele
        fim[i] = cim

    r, c = tdf.shape
    
    slopes = []
    tslopes = []
    alignment = []
    talignment = []
    rsquared = []
    for i in range(r):
        
        workingImage = fskele[tdf['file number'][i]][tdf['channel number'][i]][tdf['frame number'][i]][tdf['stack number'][i]]
        y = np.uint16(tdf['y-centroids'][i])
        x = np.uint16(tdf['x-centroids'][i])
        workingImage = workingImage[y-7:y+8, x-7:x+8]

        a = np.nonzero(workingImage)
        m, b = np.polyfit(a[1], a[0], 1)
        r2 = skl.metrics.r2_score(a[0], m*a[1]+b)
        rsquared.append(r2)
        
        oDiff = np.abs(((tdf['orientation'][i] % pi) - ((np.arctan2(m, 1) +2*pi) % pi)))
        
        if oDiff > pi/2:
            oDiff = pi - oDiff
            
        tslopes.append(m)
        talignment.append(oDiff)
        
        if (r2 >= 0.70) and (tdf['aspect ratio'][i]>1):  #############  Filtering 
            slopes.append(m)
            alignment.append(oDiff)
        else:
            slopes.append(np.nan)
            alignment.append(np.nan)
        
    return slopes, alignment, tslopes, talignment, rsquared, fskele, fim





################ Plotting puncta - Surya##################################


def punctaPlotter(i, fim, tdf, tslopes, talignment, r2):
    
        plt.figure()
        xset = np.linspace(0,12,200)
        yset = (xset-7)*np.tan(tdf['orientation'][i])+7
        yset2 = (xset-7)*tslopes[i]+7
        
        image = fim[tdf['file number'][i]][tdf['channel number'][i]][tdf['frame number'][i]][tdf['stack number'][i]]
        y = np.uint16(tdf['y-centroids'][i])
        x = np.uint16(tdf['x-centroids'][i])
        
        # finds the indices of the brightest puncta
        maxValue = np.max(image[y-7:y+8, x-7:x+8])
        maxR, maxC = np.where(image[y-7:y+8, x-7:x+8] == maxValue)
    
        # adjusts puncta coords accordingly
        y = (y + maxR[0] - 7)
        x = (x + maxC[0] - 7)

        image = image[y-7:y+8, x-7:x+8]
        
        plt.ylim(12,0)
        plt.xlim(0,12)
        plt.xlabel('index = ' + str(i) + '   |   '+'r^2 = ' + str(r2[i]))
        plt.title('blue = droplet, orange = mito')
        plt.annotate('    Angle between = ' + str(np.round(talignment[i], 3))  + ' rad',(12,7))
        
        plt.imshow(image)
        plt.plot(xset, yset,label="puncta")
        plt.plot(xset, yset2, label='mitochondria')
        plt.legend()
        plt.show()
        
        return


############################ Runing above functions ##############################

pwd=r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\Data_Analysis\Size_analysis\Size_partitonCoeff_AspectRatio\All_cells\data'

o, n, m, l, q = drop.condensate(pwd, 
                                0.0488706004242 , 
                                intensities = ['Percentile', 95], 
                                timesteps = [1], 
                                zstacks = [1], 
                                channels = [2], 
                                windowsizes = [5], 
                                masksizes = [4], 
                                backints = 'Calculated')
                                



slopes, alignment, tslopes, talignment, rsquared, fskele, images = mitoOrientation(pwd, l)



#######################look at some droplets####################### 
for i in range(20):
    punctaPlotter(i, images, l, tslopes, alignment, rsquared)


plt.hist(alignment, bins=50)
plt.title('angle variation')
plt.xlabel('Angle (rad)')
plt.ylabel('Count')
plt.show()

#plotting puncta size vs angle
plt.scatter(l['puncta major'], alignment, color='green', alpha=0.5)
plt.xlabel('Punctr major (um)')
plt.ylabel('Angle (rad)')
plt.title('Puncta major vs Angle')
plt.show()

#plotting puncta size vs angle
plt.scatter(l['puncta minor'], alignment, alpha=0.7)
plt.xlabel('Punctr minor (um)')
plt.ylabel('Angle (rad)')
plt.title('Puncta minor vs Angle')
plt.show()





################ Plotting polar histogram - Sanjaya############################

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




#from alignment data remove NaN values 
data_noNaN = [angle for angle in alignment if not math.isnan(angle)]
data = np.array(data_noNaN)


data_df = pd.DataFrame(data_noNaN)
data_df.to_csv('angles.csv')

#get histogram data
bins_num = 20
counts , bins = np.histogram(data, bins=bins_num)
width = bins[1]-bins[0]
theta = bins[:-1]+width/2 #center the bar to each bin

#create the polar figure 
ax = plt.subplot(111, projection='polar')

ax.bar(x=theta, height=counts, width=width, bottom=0, alpha=0.5, color='green')


# Set the theta limit to show only 0 to 90 degrees
ax.set_thetamax(90)
ax.set_thetamin(0)

# Set the zero location to the top - N", "NW", "W", "SW", "S", "SE", "E", or "NE"
#ax.set_theta_zero_location(90)
ax.set_theta_direction(-1)
ax.set_theta_offset((pi/2))

N=len(data)
precent = round((len(data)/len(alignment))*100 , 1)
ax.set_title(f'Polar histogram of angles \n {precent}% of the time/data in this plot, N={N}')

ax.tick_params(axis='both', which='both', labelsize=14)



plt.tight_layout()
plt.savefig('polar_histogram.pdf')
plt.show()

#plot a histogram to verify the plotting 
plt.hist(data, bins=bins_num)
plt.title('angle variation to compare with polar histogram')
plt.xlabel('Angle (rad)')
plt.ylabel('Count')
plt.show()




