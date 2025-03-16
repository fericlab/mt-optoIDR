#dynamics analysis plotting
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.stats import norm



#Loading data
data_WT = pd.read_csv('C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/S6/van-Hove-plots/data_WT_xAlined.csv')
data_Drp1 = pd.read_csv('C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/S6/van-Hove-plots/data_Drp1k38a_xAlined_mannually.csv')

#van Hover plots
def vanHoverPlot(lagdata, lagT, dimention, nbins, Absolute_values, gaussian, plotting):
    
    # this function pool all the data from the dimention you want
    
    d_total = []
   
    #isolate data for a given lagt
    data_for_lagt = lagdata[lagdata['dt']==lagT]
    
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
    #get dx,y data 
    d=d_total
    
    if Absolute_values == True:
        d=np.abs(d) #get absolute values
    
   
    #with data
    # Calculate histogram data with density=True to normalize
    counts, bin_edges = np.histogram(d, bins=nbins, density=True)

    # Calculate bin centers for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2 #(lastvalue-firstvalue)/2
    
    ##### Plotting#######
    
    if plotting == True:
        #try to fit a Gaussian distribution to data 
        mean, std = norm.fit(d) #Naormal fit (gaussian)--> get mean and std
        #plot the fucntion 
        x_values = np.linspace(min(d), max(d), 1000)
        gaussian_curve = norm.pdf(x_values, mean, std)
        
        #this is in terms of histogram
        
        plt.hist(d, bins=nbins, density = True) #density=True --> normalized probability, area under the cruve is =1
        plt.yscale('log')
        plt.xlabel('Displacement of x (\u03BCm)')
        plt.ylabel('log(probability)')
        N=len(d)
        plt.title(f'Hist: Pooled {dimention}: van Hove plot for lagt = {lagT}, nbins={nbins}, N={N}')
        plt.show()
        
        fig, ax = plt.subplots()
        if gaussian==True:
            ax.plot(x_values, gaussian_curve, label='Gaussian Fit', color='red')
        #counts=np.log10(counts) #make each value log
        ax.plot(bin_centers, counts, '-', marker='o', color='gray', label='Data')  # 'o' creates data points only
        ax.set_yscale('log')  
        ax.set_xlabel(f'Displacement of {dimention} (\u03BCm)')
        ax.set_ylabel('log (Probability)')
        N=len(d)
        ax.set_title(f' Pooled {dimention} van Hove plot for lagt = {lagT}, nbins={nbins}, N={N}')
        ax.legend()
        #plt.xlim(0,1.2)
        #plt.ylim(0,10)
        plt.show()
    
    

    return d_total, bin_centers, counts



#################### Plotting All graphs on the same canvas#################

#plot design parameters
font_properties = {
    'family': 'Arial',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

lgt=2
for i in range (50):
    lagt=lgt
    Absolute_values = True
    gausian = False
    plott = False
    nbins=20
    
    #WT
    d, x_wt_x, y_wt_x = vanHoverPlot(lagdata=data_WT, 
             lagT=lagt, 
             dimention='x', 
             nbins=nbins,
             Absolute_values=Absolute_values,
             gaussian=gausian,
             plotting=plott)
    
    d, x_wt_y, y_wt_y = vanHoverPlot(lagdata=data_WT, 
             lagT=lagt, 
             dimention='y', 
             nbins=nbins,
             Absolute_values=Absolute_values,
             gaussian=gausian,
             plotting=plott)
    
    #Drp1k38a  
    d, x_drp1_x, y_drp1_x = vanHoverPlot(lagdata=data_Drp1, 
             lagT=lagt, 
             dimention='x', 
             nbins=nbins,
             Absolute_values=Absolute_values,
             gaussian=gausian,
             plotting=plott)
    
    d, x_drp1_y, y_drp1_y = vanHoverPlot(lagdata=data_Drp1, 
             lagT=lagt, 
             dimention='y', 
             nbins=nbins,
             Absolute_values=Absolute_values,
             gaussian=gausian,
             plotting=plott)
    
    
    plt.figure(figsize=(5, 4))
    plt.plot(x_wt_x, y_wt_x, label='WT$_x$', alpha=1, marker='s', color='#7E287E', linewidth=1)
    plt.plot(x_wt_y, y_wt_y, label='WT$_y$',  alpha=1, marker='s', color='gray', linewidth=1)
    
    plt.plot(x_drp1_x, y_drp1_x, label='Drp1$^{K38A}_x$',  alpha=1, marker='o', color='#7E287E', linewidth=1)
    plt.plot(x_drp1_y, y_drp1_y, label='Drp1$^{K38A}_y$',  alpha=1, marker='o', color = 'gray', linewidth=1)
    
    plt.xlabel('Displacement (\u03BCm)', **font_properties)
    plt.ylabel('Probability', **font_properties)
    plt.yscale('log')
    plt.title(f'Overlayed van Hove Plots lagt:{lagt}, nbins={nbins}')
    plt.xlim(0,2.1)
    plt.ylim(0.04,14)
    plt.legend(frameon=False, prop = font_properties ) 
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.annotate(f'Lag time τ = {lagt}s', xy=(0.05,0.9), xycoords="axes fraction", fontsize=14)
    plt.tight_layout()
    
    #saving nessecory images 
    if lagt in [2, 10, 20, 50, 100]: #if lag t is 2, 10 or 20 it will save the pdf
        plt.savefig(f'van_Hove_plot_lagt{lagt}.pdf')
        
        
        ### Export data 
        van_Hove_plot_data = pd.DataFrame({'Displacement_WT_x': x_wt_x,
                                           'Prob_WT_x': y_wt_x,
                                           
                                           'Displacement_WT_y': x_wt_y,
                                           'Prob_WT_y': y_wt_y,
                                           
                                           'Displacement_Drp1_x': x_drp1_x,
                                           'Prob_Drp1_x': y_drp1_x,
                                           
                                           'Displacement_Drp1_y': x_drp1_y,
                                           'Prob_Drp1_y': y_drp1_y
                                           })

        van_Hove_plot_data.to_csv(f'van_Hove_plot_Data_LagT{lagt}s.csv')
        
    plt.show()
    

    lgt+=2











