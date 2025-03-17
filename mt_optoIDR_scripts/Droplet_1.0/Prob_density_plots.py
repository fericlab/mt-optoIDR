import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import gaussian_kde
import numpy as np

df=pd.read_csv(r'C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/Data_Analysis/PlayGround/droplets2.0/data_droplet2.0.csv')

#save or no, if save type 'save'
save='save'
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

def size_Denhist_major():
    
    
    #import data
    data=df['puncta major']
    data= data*1000 #convert to nm
    
    #Gaussian kernel 
    bw=0.3 #bandwidth for the Gaussian kernel
    kde = gaussian_kde(data)
    kde.set_bandwidth(bw)

    #get the data space
    data_space = np.linspace(50, data.max(), 1000)
    pdf_values = kde.evaluate(data_space) # estimate probability densities 

    bins=50
    
    fig, ax=plt.subplots(figsize=(6,4,)) #units in inches, size 
    
    c='#404042'
    
    ax.hist(data, density=True, bins=bins, alpha=0.4, color=c, label=f'histogram, bins={bins}', histtype='bar')
    ax.plot(data_space, pdf_values, label=f'Gaussian KDE bw={bw}', color=c)
    
    #fill the kde
    ax.fill_between(data_space, pdf_values, alpha=0.2, color=c)
    ax.legend()
    
    #limits
    #plt.xlim(0.500)
    ax.set_xlabel('Droplet size (nm)', **font_properties)
    ax.set_ylabel('Probability density', **font_properties)
    
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)  # Font size for y-axis tick labels


    #Mode based on the kde
    mode=data_space[np.argmax(pdf_values)]
    ax.axvline(mode)
    mode=round(mode,3)
    
    N=len(data)
    
    #std
    std=round(data.std(), 3)
    ax.annotate(f'STD={std} \nN={N}', xy=(50,0.004))
    
    
    #plt.legend()
    size_mean = round((df['puncta major'].mean() * 1000), 3)
    ax.set_title(f'Size distribution_major_ax Mean={size_mean} nm, mode:{mode}')
    
     
    #plt.grid(True, linestyle='--', linewidth=0.2, color='#d3d3d3')
    plt.show()
    #save
    if save=='save':
        fig.savefig(fname='Size_Den_distibution_major_axis.svg', dpi=600, bbox_inches='tight')
   

def size_Denhist_minor():
    #import data
    data=df['puncta minor']
    data= data*1000 #convert to nm
    
    #Gaussian kernel 
    bw=0.3 #bandwidth for the Gaussian kernel
    kde = gaussian_kde(data)
    kde.set_bandwidth(bw)

    #get the data space
    data_space = np.linspace(100, data.max(), 1000)
    pdf_values = kde.evaluate(data_space) # estimate probability densities 

    bins=50
    
    fig, ax=plt.subplots(figsize=(6,4,)) #units in inches, size 
    
    
    c='#404042'
    ax.hist(data, density=True, bins=bins, alpha=0.4, color=c, label=f'histogram, bins={bins}', histtype='bar')
    ax.plot(data_space, pdf_values, label=f'Gaussian KDE bw={bw}', color=c)

    #fill the kde
    ax.fill_between(data_space, pdf_values, alpha=0.2, color=c)
    ax.legend()
    
    #limits
    #plt.xlim(0.500)
    ax.set_xlabel('Droplet size (nm)', **font_properties)
    ax.set_ylabel('Probability density', **font_properties)
    
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)  # Font size for y-axis tick labels

    
    #Mode based on the kde
    mode=data_space[np.argmax(pdf_values)]
    ax.axvline(mode)
    
    N=len(data)
    
    #std
    std=round(data.std(), 3)
    ax.annotate(f'STD={std} \nN={N}', xy=(110,0.01))
    
    #plt.legend()
    size_mean = round((df['puncta minor'].mean() * 1000), 2)
    ax.set_title(f'Size Den_distribution_minor_ax Mean:{size_mean} nm, mode:{mode}')
    
    #plt.grid(True, linestyle='--', linewidth=0.5, color='#d3d3d3')
    plt.show()
    #save
    if save=='save':
        fig.savefig(fname='Size_Dendistibution_minor_axis.svg', #file name
                dpi=600, bbox_inches='tight')
    


def size_Denhist_major_minor():
    #import data, #convert to nm
    data_major=df['puncta major'] *1000
    data_minor=df['puncta minor'] *1000
     
    
    #Gaussian kernel 
    bw=0.3 #bandwidth for the Gaussian kernel
    kde_major = gaussian_kde(data_major)
    kde_minor = gaussian_kde(data_minor)
    
    kde_major.set_bandwidth(bw)
    kde_minor.set_bandwidth(bw)
    
    #get the data space
    data_space = np.linspace(100, data_major.max(), 1000) #data space of droplet size
    
    # estimate probability densities 
    pdf_values_major = kde_major.evaluate(data_space)
    pdf_values_minor = kde_minor.evaluate(data_space)

    bins=50
    
    
    fig, ax=plt.subplots(figsize=(5,4)) #units in inches, size 
    
    #plot_major axis
    c_major='#564F8E'
    ax.hist(data_major, density=True, bins=bins, alpha=0.4, color=c_major, label=f'histogram_major, bins={bins}', histtype='bar', range=(100,data_major.max()))
    ax.plot(data_space, pdf_values_major, label=f'Gaussian KDE bw={bw}', color=c_major)
    ax.fill_between(data_space, pdf_values_major, alpha=0.2, color=c_major) #fill the kdes
    
    #plot_minor axis
    c_minor='#161616'
    ax.hist(data_minor, density=True, bins=bins, alpha=0.4, color=c_minor, label=f'histogram_minor, bins={bins}', histtype='bar', range=(100,data_major.max()))
    ax.plot(data_space, pdf_values_minor, label=f'Gaussian KDE bw={bw}', color=c_minor)
    ax.fill_between(data_space, pdf_values_minor, alpha=0.2, color=c_minor) #fill the kdes
    
    
    ax.legend(frameon=False)
    
    #limits
    #plt.xlim(0.500)
    ax.set_xlabel('Droplet size (nm)', **font_properties)
    ax.set_ylabel('Probability density', **font_properties)
    
    #Specify the number of ticks 
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)  # Font size for y-axis tick labels

    #add some detalis
    size_mean_major = round((df['puncta major'].mean() * 1000), 2)
    size_mean_minor = round((df['puncta minor'].mean() * 1000), 2)
    
    ax.set_title(f'Size distribution \nMajor Mean: {size_mean_major} nm \nMinor Mean: {size_mean_minor} nm')
    
    #plt.grid(True, linestyle='--', linewidth=0.5, color='#d3d3d3')
    plt.show()
    #save
    if save=='save':
        fig.savefig(fname='Size_Den_distibution_major_minor.pdf', #file name
                dpi=600, bbox_inches='tight')
    


def aspectR_Denhist():
    
    #import data
    data=df['aspect ratio']
    
    #Gaussian kernel 
    bw=0.3 #bandwidth for the Gaussian kernel
    kde = gaussian_kde(data)
    kde.set_bandwidth(bw)

    #get the data space
    data_space = np.linspace(0, data.max(), 1000)
    pdf_values = kde.evaluate(data_space) # estimate probability densities 

    bins=50
    
    
    #plot aspect ratio
    fig, ax=plt.subplots(figsize=(5,4)) 
   
    c='#564F8E' #color
    ax.hist(data, density=True, bins=bins, alpha=0.4, color=c, label=f'histogram, bins={bins}', histtype='bar')
    ax.plot(data_space, pdf_values, label=f'Gaussian KDE bw={bw}', color=c)

    #fill the kde
    ax.fill_between(data_space, pdf_values, alpha=0.2, color=c)
    ax.legend(frameon=False)
    
    #Specify the number of ticks 
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)

    ax.set_xlabel('Aspect ratio', **font_properties)
    ax.set_ylabel('Probability density', **font_properties)
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)  # Font size for y-axis tick labels
    ax.set_ylim(0,2.5)
    ax.set_xlim(0.5,3.3)
    #plt.grid(True, linestyle='--', linewidth=0.5, color='#d3d3d3')
    
    #Mode based on the kde
    mode=data_space[np.argmax(pdf_values)]
    ax.axvline(mode)
    
    N=len(data)
    
    #std
    std=round(data.std(), 3)
    ax.annotate(f'STD={std} \nN={N}', xy=(2.5,1.25))
    
    aspectR_mean=round(df['aspect ratio'].mean(), 2)
    ax.set_title(f'Aspect ratio distribution, Mean:{aspectR_mean}, Mode:{mode}')
    plt.show()
    
    if save=='save':
        fig.savefig(fname='Aspect_ration_Den_hist.pdf', dpi=600, bbox_inches='tight')
    
    
   

def partition_coeff_Den():
    #import data
    data=df['partition coeff']
    
    #Gaussian kernel 
    bw=0.3 #bandwidth for the Gaussian kernel
    kde = gaussian_kde(data)
    kde.set_bandwidth(bw)

    #get the data space
    data_space = np.linspace(0, data.max(), 1000)
    pdf_values = kde.evaluate(data_space) # estimate probability densities 

    bins=50
    
    
    #plot aspect ratio
    fig, ax=plt.subplots(figsize=(5,4)) 
   
    c='#564F8E' #color
    ax.hist(data, density=True, bins=bins, alpha=0.4, color=c, label=f'histogram, bins={bins}', histtype='bar')
    ax.plot(data_space, pdf_values, label=f'Gaussian KDE bw={bw}', color=c, alpha=1)

    #fill the kde
    ax.fill_between(data_space, pdf_values, alpha=0.2, color=c)
    ax.legend()
    
    #Specify the number of ticks 
    ax.locator_params(axis='x', nbins=6)
    ax.locator_params(axis='y', nbins=6)
    
    
    ax.set_xlabel('Partition coefficients', **font_properties)
    ax.set_ylabel('Probability density', **font_properties)
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)  # Font size for y-axis tick labels
    ax.set_ylim(0,0.75)
    ax.set_xlim(1,9.5)
    #Mode based on the kde
    mode=data_space[np.argmax(pdf_values)]
    ax.axvline(mode)
    
    N=len(data)
    
    #std
    std=round(data.std(), 3)
    ax.annotate(f'STD={std} \nN={N}', xy=(7,0.3))
    
    partcoef_mean=round(data.mean(), 2)
    ax.set_title(f'Partition coeff distribution, mean:{partcoef_mean}, Mode:{mode}')
    
    #plt.grid(True, linestyle='--', linewidth=0.5, color='#d3d3d3')
    plt.show()
    #save
    if save=='save':
        fig.savefig(fname='Partition_coeff_Den_distribution.svg', #file name
                dpi=600, bbox_inches='tight')
    



size_Denhist_major()
size_Denhist_minor()
size_Denhist_major_minor()
aspectR_Denhist()
partition_coeff_Den()


######################################
print('Sample size (Droplet number) is : ', len(df['puncta major']))



################################Tests for the normal distribution################


