import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r'C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/S3/Images/Act/S3_Int_data_crop2.csv')

channels = ['DNA', 'DDX28', 'GRSF1']

colors=['violet', 'green', 'red']

for i in range (len(channels)):
    channel = channels[i]
    intensities = df[channel]
    
    x=df.iloc[0:len(intensities), 0] #get the distance from the first column untill the last raw of i th channel
    
    norm_intensities = intensities/intensities.max()
    
    clr=colors[i]
    plt.plot(x, norm_intensities, label=channel, alpha=0.6, linewidth=2, color=clr)
    plt.xlabel('Distance (µm)', fontsize=18)
    plt.ylabel('Normalized intensity (a.u)', fontsize=18)
   
    


# Increase the font size of x and y axis ticks
plt.tick_params(axis='both', which='major', labelsize=20)
plt.legend(frameon=False)
#plt.title('Normalized intensity')
plt.ylim(0,1.2)
plt.savefig(r'intensity_profile_S3_crop2.svg')
plt.show()





