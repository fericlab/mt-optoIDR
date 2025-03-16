# -*- coding: utf-8 -*-
"""
Created on Thu May  2 10:44:58 2024

@author: sanja
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42
tick_font_size = 16


#load data

df = pd.read_csv('Int_values_normalized.csv')


#get the x and y axis data
x=df['Distance_uM']

y1=df['Normalized']
#y2=df['Normalized_PostAct']
#y3=df['FUS']

label_list=['Post-Activation'] #give your lables

#Labels 
x_label='Distance (\u00B5m)'
y_label='Normalized Intensity (a.u)'

graph_title='Distance vs intensity'

#save as
file_name='intensity_plots_droplet.svg' #jeg, png, svg etc
dpi_value=600

#plot design parameters
font_properties = {
    'family': 'Arial',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}


#plot the figure
plt.figure(figsize=(5,4)) #units in inches, size = A7
plt.style.use('tableau-colorblind10')

#plot multiple graphs
plt.plot(x, list(zip(y1)), #x and y axis
         linewidth=2, marker='', 
         linestyle='-',  
         alpha=0.5,
         label=label_list, color='#99AF82') #give your lables


#x and Y axis limits
plt.xlim(0, 1)
plt.ylim(0.05, 1.2)

#specify ticklabel space ,Ex: plt.xticks(np.arange(0, 1, 0.2)) from 0-1, space 0.2
plt.xticks(np.arange(0, 1.2, 0.2))
#grid



#axis


#x and y labels
plt.xlabel(x_label, **font_properties)
plt.ylabel(y_label, **font_properties)


plt.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
plt.tick_params(axis='y', labelsize=tick_font_size)

#title
#plt.title(graph_title, **font_properties)

#save
plt.savefig('intensity_profile.pdf')
plt.show()

