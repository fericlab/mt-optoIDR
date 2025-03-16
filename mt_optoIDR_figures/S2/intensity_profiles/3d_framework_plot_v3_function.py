# -*- coding: utf-8 -*-
"""
Created on Thu May  2 09:39:57 2024

@author: sanja
"""


import matplotlib.pyplot as plt
import imageio
import cv2
import numpy as np




#this makes text editable
plt.rcParams['pdf.fonttype'] = 42
tick_font_size = 14




def three_d_wireframe_plot (file_name, clr_hex, font_size, save_name):
    
    """
    Input:
    file_name = droplet file (should be a png for this)
    clr_hex = color for the wireframe, hex code 
    font_size = font size for tick labels 
    save_name = name to save the file 
    
    Return:
        two images , original image, and wireframe image with background (image)
    
    Example 
        three_d_wireframe_plot('droplet_0.png', '#32CD32' , 8, '3d_surface_of_the_droplet')
    """
    
    #import the image as an png --> written for png
    image = imageio.v2.imread(str(file_name))
    #convert into gryascale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    plt.imshow(image)
    plt.tight_layout()
    plt.savefig(save_name+'Original_image', dpi=600)
    
    #get the Y and X
    height, width = gray_image.shape
    
    # Create meshgrid , multiply this grid by the pixel to micron conversion factor to get the grid in micron
    X, Y = np.meshgrid(np.arange(width)*0.048, np.arange(height)*0.048)
    Z = np.zeros_like(X)
    float_image = gray_image.astype(float) /(np.max(gray_image)) #Normalization for the surface plot as facecolors need 0-1 
    
    #plot design parameters
    font_properties = {
        'family': 'Arial',
        'size': 15,
        'weight': 'normal',
        'style': 'normal'
    }

    
    # Plot image as , 1. wireframe, 2.surface (background)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    
    # make the panes transparent or change color (OPTIONAL)
    # Example: ax.xaxis.set_pane_color('red') or use RGBA values 
    #ax.xaxis.set_pane_color('lightgray')
    #ax.yaxis.set_pane_color('lightgray')
    #ax.zaxis.set_pane_color('lightgray')
    
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0))
    
    # make the grid lines transparent or change color (OPTIONAL)
    ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
     
    # Add wireframe
    ax.plot_wireframe(X, Y, float_image, color=str(clr_hex), alpha=0.6) #plot wireframe for the normalized intensities 
    
    #plot the surface
    ax.plot_surface(X, Y, Z, facecolors=plt.cm.gray(float_image), alpha=1) #for facecolor need data in RGBA format , convert if not
    
    #x and y labels 
    ax.set_xlabel('Distance (\u00B5m)', **font_properties)
    ax.set_ylabel('Distance (\u00B5m)', **font_properties)
    ax.set_zlabel('Normalized Intensity (a.u)', **font_properties)
    
    # Rotate tick labels and set font size
    ax.xaxis.set_tick_params(rotation=0, labelsize=font_size)
    ax.yaxis.set_tick_params(rotation=-20, labelsize=font_size)
    ax.zaxis.set_tick_params(labelsize=font_size)
    

   #chnage the z axis position to left
    ax.zaxis.set_ticks_position('lower')
    ax.zaxis.set_label_position('lower')
    
    ax.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
    ax.tick_params(axis='y', labelsize=tick_font_size)
    ax.tick_params(axis='z', labelsize=tick_font_size)

    #save the image 
    plt.tight_layout()
    plt.savefig(str(save_name+'.pdf'))
    plt.show()
    
    return X, Y, float_image


x, y, floatImage = three_d_wireframe_plot('droplet_0.png', '#9AB083' , 8, '3d_surface_of_the_droplet')


