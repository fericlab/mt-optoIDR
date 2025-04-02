import czifile as czifile
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import pandas as pd
import trackpy as tp
import cv2
from math import pi
import warnings
import os

def condensate (folder, p_to_mu, **static_conditions):
    
    # Given an input file, returns size data and centroid estimates for an
    # inputted threshold. The program previews the image with the thresholded
    # punctae and allows for modifications of the threshold on the fly. 
    #
    # Inputs:
    # imagefile     - [directory] complete file directory of the input image. takes .czi files
    # channel       - [int] channel index of .czi file
    # p_to_mu       - [float] pixel to micron conversion ratio
    # I_thr         - [float] initial threshold value. increase to make more specific and 
    #                 decrease to gather more puncta
    # back_i        - [float] estimate of the cytoplasmic background signal
    # p_size        - [int] approximate radius (in pixels) of the puncta. Overestimate
    # wndw          - [int] window size from the center to the size taken around each puncta, in pixels
    # interactive   - [1/0] 1 means user output is not suppressed and 0 suppresses all
    #                 outputs
    # Outputs:
    # centroids     - returns the x and y coordinates of the centroid
    # pixelAxis     - returns the major and minor axis lengths in pixels
    # punctaAxis    - returns the major and minor axis lengths in microns
    # aspectRatio   - returns the aspect ratio of the punctae
    # threshold     - returns the final threshold value used for the image
    
    # creates a list of all of the files in a given folder directory
    f = os.listdir(folder)

    totaldicts = []
    rawtotaldicts = []
    
    files = []
    
    timesteps = static_conditions.pop('timesteps', [])
    zstacks = static_conditions.pop('zstacks', [])
    channel = static_conditions.pop('channels', [])
    iArray = static_conditions.pop('intensities', [])
    pArray = static_conditions.pop('masksizes', [])
    bArray = static_conditions.pop('backints', [])
    wArray = static_conditions.pop('windowsizes', [])
    cropcoords = static_conditions.pop('crop',[])
    suppression = static_conditions.pop('plot', [])
    
    for i in range(len(f)):
        if f[i].endswith('.czi'): # may be unnecesary, removes a hidden file #sanjaya
            files.append(folder + '/' + f[i])
        else:
            pass
            
    # gaussian curve
    def gaussian2D(data, amp, x0, y0, a, b, c, d):
        x, y = data
        x0 = float(x0)
        y0 = float(y0)
        g = d + amp*np.exp( - (a*((x-x0)**2) + b*(x-x0)*(y-y0) + c*((y-y0)**2)))                                   
        return g.ravel()

    # creates dictionary for files
    totallist = {}
    arraylistt = {}
    imagelistt = {}
    # for each file in the folder, the program runs
    for j in range(len(files)):
        
        # gets user inputs on relevant parameters that need to be studies
        print("\nThe current file is "+ files[j] + '\n')
                  
        if timesteps == []:
            f1 = int(input("Enter beginning frame #: \n" ))
            fn = int(input("\nEnter final frame #: \n" ))
        else:
            f1 = timesteps[0]
            fn = timesteps[-1]
        if zstacks == []:
            z1 = int(input("\nEnter beginning z-stack #: \n" ))
            zn = int(input("\nEnter final z-stack #: \n" ))
        else:
            z1 = zstacks[0]
            zn = zstacks[-1]
        if channel == []:
            c = input("\nEnter relevant channel numbers (seperated by commas): \n")
            c = c.split(',')
            channels = np.uint16(np.array(c)) - 1
        else:
            channels = np.array(channel) - 1
            
        # creates arrays which will eventually be used iterate the program
        frames = np.uint16(np.linspace(f1-1, fn-1, fn-f1+1))
        zstack = np.uint16(np.linspace(z1-1, zn-1, zn-z1+1)) 
        
        # creates dictionaries for channels
        channellist = {}
        arraylistc = {}
        imagelistc = {}
        for k in channels:
            
            # reads the image and pulls out a frame and a channel
            filedata = czifile.imread(files[j])
            
            # creates the proper format for each frame of the image
            if len(filedata.shape) == 6:
                image = filedata[f1-1,int(k),0,:,:,0]
            elif len(filedata.shape) == 8:
                image = filedata[f1-1,0,int(k),0,int(z1-1),:,:,0]
                
            plt.title('Original Image')
            plt.imshow(image)
            plt.pause(0.5)
            
            if cropcoords == []:
                pass
            elif len(cropcoords) == 4:
                ymin = np.uint64(cropcoords[2])
                ymax = np.uint64(cropcoords[3])
                xmin = np.uint64(cropcoords[0])
                xmax = np.uint64(cropcoords[1])
                croppedimage = image[ymin:ymax,xmin:xmax]
                
                plt.title('Cropped Image')
                plt.imshow(croppedimage)
                plt.pause(0.5)
                
            else:
                ymin = np.uint64(cropcoords[1]-cropcoords[2])
                ymax = np.uint64(cropcoords[1]+cropcoords[2]+1)
                xmin = np.uint64(cropcoords[0]-cropcoords[2])
                xmax = np.uint64(cropcoords[0]+cropcoords[2]+1)
                croppedimage = image[ymin:ymax,xmin:xmax]
                
                plt.title('Cropped Image')
                plt.imshow(croppedimage)
                plt.pause(0.5)
                
            if iArray == []:
                
                I_thr = float(input("\nInitial Threshold Guess is: \n"))
                # sets the sentinel condition to ensure that the while loop runs at least once
                sen_condition = 0
                # initial segmentation, and user input
                while sen_condition != 1: 
                    
                    # creates the proper format for each frame of the image
                    if len(filedata.shape) == 6:
                        image = filedata[f1-1,int(k),0,:,:,0]
                    elif len(filedata.shape) == 8:
                        image = filedata[f1-1,0,int(k),0,int(z1-1),:,:,0]
                        
                    # finds the shape of the image
                    h, w = image.shape
                    
                    # uses trackpy to locate the puncta within the image matrix
                    seg = tp.locate(image, 11, threshold = I_thr)
                    seg1 = np.zeros((len(seg['x']), 2))
                    seg1[:,0] = np.uint16(seg['y'])
                    seg1[:,1] = np.uint16(seg['x'])
                    
                    if cropcoords == []:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<np.uint64(h/30), :] = 0
                        seg1[seg1[:,1]<np.uint64(h/30), :] = 0
                        seg1[seg1[:,0]>h-np.uint64(h/30), :] = 0
                        seg1[seg1[:,1]>h-np.uint64(h/30), :] = 0
                    elif len(cropcoords) == 4:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<cropcoords[2], :] = 0
                        seg1[seg1[:,1]<cropcoords[0], :] = 0
                        seg1[seg1[:,0]>cropcoords[3], :] = 0
                        seg1[seg1[:,1]>cropcoords[1], :] = 0
                    else:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<cropcoords[1]-cropcoords[2], :] = 0
                        seg1[seg1[:,1]<cropcoords[0]-cropcoords[2], :] = 0
                        seg1[seg1[:,0]>cropcoords[1]+cropcoords[2]+1, :] = 0
                        seg1[seg1[:,1]>cropcoords[0]+cropcoords[2]+1, :] = 0

                    
                    seg1 = seg1[seg1[:,0]!=0]
                    
                    # preallocates punctae centroid estimates
                    seg2 = np.zeros((len(seg1),2))
                    
                    # centers the image based on the brightest spot 
                    for i in range(len(seg1)):
                        
                        # find max and min windows
                        xmin = np.uint16(seg1[i,0])-5
                        xmax = np.uint16(seg1[i,0])+5
                        ymin = np.uint16(seg1[i,1])-5
                        ymax = np.uint16(seg1[i,1])+5
                        
                        # finds the indices of the brightest puncta
                        maxValue = np.max(image[xmin:xmax,ymin:ymax])
                        maxR, maxC = np.where(image[xmin:xmax,ymin:ymax] == maxValue)
                    
                        # adjusts puncta coords accordingly
                        seg2[i, 0] = np.float64(seg1[i, 0] + maxR[0] - 5)
                        seg2[i, 1] = np.float64(seg1[i, 1] + maxC[0] - 5)
                
                    # removes nonzero elements (from punctae that were too close to edge)
                    seg2 = seg2[np.matrix.nonzero(seg2[:,1])]
            
                    # shows the image matrix, such that the user can change the threshold
                    plt.title('Segmentation')
                    plt.imshow(image)
                    plt.scatter(seg2[:,1],seg2[:,0], facecolors='none', edgecolors='r')
                    plt.pause(0.01)
                
                    # determines whether the loop will continue
                    x = input("\nIf the threshold condition is good enter 1: \n")
                    sen_condition = int(x)
                    
                    # asks for new threshold
                    if sen_condition != 1:
                        I_thr = int(input("\nEnter new threshold: \n"))
            elif iArray[0] == "Percentile":
                I_thr = tp.percentile_threshold(image, iArray[1])
                # finds the shape of the image
                h, w = image.shape
                
                # uses trackpy to locate the puncta within the image matrix
                seg = tp.locate(image, 11, threshold = I_thr)
                seg1 = np.zeros((len(seg['x']), 2))
                seg1[:,0] = np.uint16(seg['y'])
                seg1[:,1] = np.uint16(seg['x'])
                
                if cropcoords == []:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,0]>h-np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]>h-np.uint64(h/30), :] = 0
                elif len(cropcoords) == 4:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0], :] = 0
                    seg1[seg1[:,0]>cropcoords[3], :] = 0
                    seg1[seg1[:,1]>cropcoords[1], :] = 0
                else:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[1]-cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0]-cropcoords[2], :] = 0
                    seg1[seg1[:,0]>cropcoords[1]+cropcoords[2]+1, :] = 0
                    seg1[seg1[:,1]>cropcoords[0]+cropcoords[2]+1, :] = 0
            
                seg1 = seg1[seg1[:,0]!=0]
                
                # preallocates punctae centroid estimates
                seg2 = np.zeros((len(seg1),2))
                
                # centers the image based on the brightest spot 
                for i in range(len(seg1)):
                    
                    # find max and min windows
                    xmin = np.uint16(seg1[i,0])-5
                    xmax = np.uint16(seg1[i,0])+5
                    ymin = np.uint16(seg1[i,1])-5
                    ymax = np.uint16(seg1[i,1])+5
                    
                    # finds the indices of the brightest puncta
                    maxValue = np.max(image[xmin:xmax,ymin:ymax])
                    maxR, maxC = np.where(image[xmin:xmax,ymin:ymax] == maxValue)
                
                    # adjusts puncta coords accordingly
                    seg2[i, 0] = np.float64(seg1[i, 0] + maxR[0] - 5)
                    seg2[i, 1] = np.float64(seg1[i, 1] + maxC[0] - 5)
            
                # removes zero elements (from punctae that were too close to edge)
                seg2 = seg2[np.matrix.nonzero(seg2[:,1])]
            
            elif len(iArray) == 1:
                I_thr = np.float64(iArray[0])
                # finds the shape of the image
                h, w = image.shape
                
                # uses trackpy to locate the puncta within the image matrix
                seg = tp.locate(image, 11, threshold = I_thr)
                seg1 = np.zeros((len(seg['x']), 2))
                seg1[:,0] = np.uint16(seg['y'])
                seg1[:,1] = np.uint16(seg['x'])
                
                if cropcoords == []:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,0]>h-np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]>h-np.uint64(h/30), :] = 0
                elif len(cropcoords) == 4:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0], :] = 0
                    seg1[seg1[:,0]>cropcoords[3], :] = 0
                    seg1[seg1[:,1]>cropcoords[1], :] = 0
                else:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[1]-cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0]-cropcoords[2], :] = 0
                    seg1[seg1[:,0]>cropcoords[1]+cropcoords[2]+1, :] = 0
                    seg1[seg1[:,1]>cropcoords[0]+cropcoords[2]+1, :] = 0
                
                seg1 = seg1[seg1[:,0]!=0]
                    
                # preallocates punctae centroid estimates
                seg2 = np.zeros((len(seg1),2))
                
                # centers the image based on the brightest spot 
                for i in range(len(seg1)):
                    
                    # find max and min windows
                    xmin = np.uint16(seg1[i,0])-5
                    xmax = np.uint16(seg1[i,0])+5
                    ymin = np.uint16(seg1[i,1])-5
                    ymax = np.uint16(seg1[i,1])+5
                    
                    # finds the indices of the brightest puncta
                    maxValue = np.max(image[xmin:xmax,ymin:ymax])
                    maxR, maxC = np.where(image[xmin:xmax,ymin:ymax] == maxValue)
                
                    # adjusts puncta coords accordingly
                    seg2[i, 0] = np.float64(seg1[i, 0] + maxR[0] - 5)
                    seg2[i, 1] = np.float64(seg1[i, 1] + maxC[0] - 5)
            
                # removes zero elements (from punctae that were too close to edge)
                seg2 = seg2[np.matrix.nonzero(seg2[:,1])]
                
            else:
                I_thr = iArray[k]
                # finds the shape of the image
                h, w = image.shape
                
                # uses trackpy to locate the puncta within the image matrix
                seg = tp.locate(image, 11, threshold = I_thr)
                seg1 = np.zeros((len(seg['x']), 2))
                seg1[:,0] = np.uint16(seg['y'])
                seg1[:,1] = np.uint16(seg['x'])
                
                if cropcoords == []:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]<np.uint64(h/30), :] = 0
                    seg1[seg1[:,0]>h-np.uint64(h/30), :] = 0
                    seg1[seg1[:,1]>h-np.uint64(h/30), :] = 0
                elif len(cropcoords) == 4:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0], :] = 0
                    seg1[seg1[:,0]>cropcoords[3], :] = 0
                    seg1[seg1[:,1]>cropcoords[1], :] = 0
                else:
                    # removes punctae that are too close to the edge
                    seg1[seg1[:,0]<cropcoords[1]-cropcoords[2], :] = 0
                    seg1[seg1[:,1]<cropcoords[0]-cropcoords[2], :] = 0
                    seg1[seg1[:,0]>cropcoords[1]+cropcoords[2]+1, :] = 0
                    seg1[seg1[:,1]>cropcoords[0]+cropcoords[2]+1, :] = 0
            
                seg1 = seg1[seg1[:,0]!=0]
                
                # preallocates punctae centroid estimates
                seg2 = np.zeros((len(seg1),2))
                
                # centers the image based on the brightest spot 
                for i in range(len(seg1)):
                    
                    # find max and min windows
                    xmin = np.uint16(seg1[i,0])-5
                    xmax = np.uint16(seg1[i,0])+5
                    ymin = np.uint16(seg1[i,1])-5
                    ymax = np.uint16(seg1[i,1])+5
                    
                    # finds the indices of the brightest puncta
                    maxValue = np.max(image[xmin:xmax,ymin:ymax])
                    maxR, maxC = np.where(image[xmin:xmax,ymin:ymax] == maxValue)
                
                    # adjusts puncta coords accordingly
                    seg2[i, 0] = np.float64(seg1[i, 0] + maxR[0] - 5)
                    seg2[i, 1] = np.float64(seg1[i, 1] + maxC[0] - 5)
            
                # removes nonzero elements (from punctae that were too close to edge)
                seg2 = seg2[np.matrix.nonzero(seg2[:,1])]
            
            if len(seg2[:,0]) != 0:
                prandom = np.random.randint(0,len(seg2[:,0]), 3)
            else:
                pass
            
            if wArray == []:
                wndw = int(input("\nInitial Window Size Guess is: \n"))
                # window size
                sen_condition = 0
                while sen_condition != 1:
                    
                    # describes the small window
                    xmin1 = np.uint64(seg2[prandom[0],0]-wndw)
                    xmax1 = np.uint64(seg2[prandom[0],0]+wndw+1)
                    ymin1 = np.uint64(seg2[prandom[0],1]-wndw)
                    ymax1 = np.uint64(seg2[prandom[0],1]+wndw+1)
                    
                    # describes the small window
                    xmin2 = np.uint64(seg2[prandom[1],0]-wndw)
                    xmax2 = np.uint64(seg2[prandom[1],0]+wndw+1)
                    ymin2 = np.uint64(seg2[prandom[1],1]-wndw)
                    ymax2 = np.uint64(seg2[prandom[1],1]+wndw+1)
                    
                    # describes the small window
                    xmin3 = np.uint64(seg2[prandom[2],0]-wndw)
                    xmax3 = np.uint64(seg2[prandom[2],0]+wndw+1)
                    ymin3 = np.uint64(seg2[prandom[2],1]-wndw)
                    ymax3 = np.uint64(seg2[prandom[2],1]+wndw+1)
                    
                    # finds a small area around puncta to fit a distribution
                    punct1 = image[xmin1:xmax1,ymin1:ymax1]
                    punct2 = image[xmin2:xmax2,ymin2:ymax2]
                    punct3 = image[xmin3:xmax3,ymin3:ymax3]
                    
                    fig,ax = plt.subplots(1, 3)
                    fig.suptitle('Window Size')
                    ax[0].axis('off')
                    ax[1].axis('off')
                    ax[2].axis('off')
                    
                    # shows the image matrix, such that the user can change the window sie
                    ax[0].imshow(punct1)
                    ax[1].imshow(punct2)
                    ax[2].imshow(punct3)
        
                    plt.pause(0.01)
                
                    # determines whether the loop will continue
                    x = input("\nIf the window size is good enter 1: \n")
                    sen_condition = int(x)
                    
                    # asks for new window size
                    if sen_condition != 1:
                        wndw = int(input("\nEnter new window size: \n")) 
            elif len(wArray) == 1:
                wndw = np.int16(wArray[0])
            else:
                wndw = wArray[k]
            
            if bArray == []:
                # INTIAL BACKGROUND ESTIMATE HERE
                # linearizes the image to work with OpenCV's adaptive thresholding
                minim = min(image.ravel())
                maxim = max(image.ravel())
                linear_image = np.multiply((image - minim), 255/(maxim-minim))
    
                # creates a new image that masks the puncta based on the adaptive thresholding scheme
                processed = cv2.GaussianBlur(linear_image, (13,13), 0)
                adaptive = cv2.adaptiveThreshold(np.uint8(processed), 255, 
                                                 cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                                                 cv2.THRESH_BINARY, 11, 0)
    
                # finds mean of this number
                backImage = (1-adaptive/255) * image
                backImageGuess = 2*sum(sum(backImage))/len(np.nonzero(backImage)[0][:])
                print("\nInitial Background Intensity Guess is" + str(backImageGuess) + "\n")
                back_i = backImageGuess
                # cytoplasmic background intesity 
                sen_condition = 0
                while sen_condition != 1:
                    
                    # describes the small window
                    xmin1 = np.uint64(seg2[prandom[0],0]-wndw)
                    xmax1 = np.uint64(seg2[prandom[0],0]+wndw+1)
                    ymin1 = np.uint64(seg2[prandom[0],1]-wndw)
                    ymax1 = np.uint64(seg2[prandom[0],1]+wndw+1)
                    
                    # describes the small window
                    xmin2 = np.uint64(seg2[prandom[1],0]-wndw)
                    xmax2 = np.uint64(seg2[prandom[1],0]+wndw+1)
                    ymin2 = np.uint64(seg2[prandom[1],1]-wndw)
                    ymax2 = np.uint64(seg2[prandom[1],1]+wndw+1)
                    
                    # describes the small window
                    xmin3 = np.uint64(seg2[prandom[2],0]-wndw)
                    xmax3 = np.uint64(seg2[prandom[2],0]+wndw+1)
                    ymin3 = np.uint64(seg2[prandom[2],1]-wndw)
                    ymax3 = np.uint64(seg2[prandom[2],1]+wndw+1)
                    
                    # finds a small area around puncta to fit a distribution
                    punct1 = image[xmin1:xmax1,ymin1:ymax1]
                    punct2 = image[xmin2:xmax2,ymin2:ymax2]
                    punct3 = image[xmin3:xmax3,ymin3:ymax3]
                    
                    # isolates mitochondria
                    punct1 = punct1*(punct1>back_i)
                    punct2 = punct2*(punct2>back_i)
                    punct3 = punct3*(punct3>back_i)
                    
                    fig,ax = plt.subplots(1, 3)
                    fig.suptitle('Mitochondrial Network Isolation')
                    ax[0].axis('off')
                    ax[1].axis('off')
                    ax[2].axis('off')
                    
                    # shows the image matrix, such that the user can change the window sie
                    ax[0].imshow(punct1)
                    ax[1].imshow(punct2)
                    ax[2].imshow(punct3)
        
                    plt.pause(0.01)
                
                    # determines whether the loop will continue
                    x = input("\nIf the absolute background isolation is good enter 1: \n")
                    sen_condition = int(x)
                    
                    # asks for new background intesity
                    if sen_condition != 1:
                        back_i = int(input("\nEnter extramitochondrial space intensity: \n")) 
            elif bArray == 'Calculated':
                # INTIAL BACKGROUND ESTIMATE HERE
                # linearizes the image to work with OpenCV's adaptive thresholding
                minim = min(image.ravel())
                maxim = max(image.ravel())
                linear_image = np.multiply((image - minim), 255/(maxim-minim))
    
                # creates a new image that masks the puncta based on the adaptive thresholding scheme
                processed = cv2.GaussianBlur(linear_image, (13,13), 0)
                adaptive = cv2.adaptiveThreshold(np.uint8(processed), 255, 
                                                 cv2.ADAPTIVE_THRESH_GAUSSIAN_C, 
                                                 cv2.THRESH_BINARY, 11, 0)
    
                # finds mean of this number
                backImage = (1-adaptive/255) * image
                backImageGuess = 2*sum(sum(backImage))/len(np.nonzero(backImage)[0][:])
                back_i = backImageGuess
                
            elif len(bArray) == 1:
                back_i = np.float64(bArray[0])
            else:
                back_i = bArray[k]
                    
                # defines the image to exclude the background intensity 
                image = image*(image>back_i)
                
            # defines variables needed to make the mask
            c = np.int16((2*wndw+1)/2)
            x = np.linspace(0,2*wndw,2*wndw+1)
            y = np.linspace(0,2*wndw,2*wndw+1)
            x, y = np.meshgrid(x,y)
            sphere = np.sqrt((wndw+3)**2-(x-c)**2-(y-c)**2)
            if pArray == []:
                p_size = int(input("\nInitial Mask Guess is: \n"))
                sen_condition = 0
                while sen_condition != 1:
                    
                    # describes the small window
                    xmin1 = np.uint64(seg2[prandom[0],0]-wndw)
                    xmax1 = np.uint64(seg2[prandom[0],0]+wndw+1)
                    ymin1 = np.uint64(seg2[prandom[0],1]-wndw)
                    ymax1 = np.uint64(seg2[prandom[0],1]+wndw+1)
                    
                    # describes the small window
                    xmin2 = np.uint64(seg2[prandom[1],0]-wndw)
                    xmax2 = np.uint64(seg2[prandom[1],0]+wndw+1)
                    ymin2 = np.uint64(seg2[prandom[1],1]-wndw)
                    ymax2 = np.uint64(seg2[prandom[1],1]+wndw+1)
                    
                    # describes the small window
                    xmin3 = np.uint64(seg2[prandom[2],0]-wndw)
                    xmax3 = np.uint64(seg2[prandom[2],0]+wndw+1)
                    ymin3 = np.uint64(seg2[prandom[2],1]-wndw)
                    ymax3 = np.uint64(seg2[prandom[2],1]+wndw+1)
                    
                    # finds a small area around puncta to fit a distribution
                    punct1 = image[xmin1:xmax1,ymin1:ymax1]
                    punct2 = image[xmin2:xmax2,ymin2:ymax2]
                    punct3 = image[xmin3:xmax3,ymin3:ymax3]
                    
                    punctBack1 = punct1*(sphere<np.sqrt((wndw+3)**2-(p_size-1)**2))
                    punctBack2 = punct2*(sphere<np.sqrt((wndw+3)**2-(p_size-1)**2))
                    punctBack3 = punct3*(sphere<np.sqrt((wndw+3)**2-(p_size-1)**2))
                    
                    fig,ax = plt.subplots(1, 3)
                    fig.suptitle('Puncta Mask Size')
                    ax[0].axis('off')
                    ax[1].axis('off')
                    ax[2].axis('off')
                    
                    # shows the image matrix, such that the user can change the threshold
                    ax[0].imshow(punctBack1)
                    ax[1].imshow(punctBack2)
                    ax[2].imshow(punctBack3)
        
                    plt.pause(0.01)
                
                    # determines whether the loop will continue
                    x = input("\nIf the mask is good enter 1: \n")
                    sen_condition = int(x)
                    
                    # asks for new bounds
                    if sen_condition != 1:
                        p_size = int(input("\nEnter new puncta size: \n"))
            elif len(pArray) == 1:
                p_size = np.int16(pArray[0])
            else:
                p_size = pArray[k]
            # creates dictionaries for frames
            framelist = {}
            arraylistf = {}
            imagelistf = {}
            for l in frames:
                # creates dictionaries for zstacks
                stacklist = {}
                arraylistz = {}
                imagelistz = {}
                
                if iArray =='Percentile':
                    if len(filedata.shape) == 6:
                        image = filedata[int(l),int(k),0,:,:,0]
                    elif len(filedata.shape) == 8:
                        image = filedata[int(l),0,int(k),0,(zn-z1)/2,:,:,0]
                    I_thr = tp.percentile_threshold(image, iArray[1])
                for m in zstack:                
                    # creates the proper format for each frame of the image
                    if len(filedata.shape) == 6:
                        image = filedata[int(l),int(k),0,:,:,0]
                    elif len(filedata.shape) == 8:
                        image = filedata[int(l),0,int(k),0,int(m),:,:,0]
                        
                    rawImage = image
                    
                    # finds the shape of the image
                    h, w = image.shape
                    
                    # uses trackpy to locate the puncta within the image matrix
                    seg = tp.locate(image, 11, threshold = I_thr)
                    seg1 = np.zeros((len(seg['x']), 2))
                    seg1[:,0] = np.uint16(seg['y'])
                    seg1[:,1] = np.uint16(seg['x'])
                    
                    if len(seg) == 0:
                        punctaRemover = [];
                        continue
                    
                    if cropcoords == []:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<np.uint64(h/30), :] = 0
                        seg1[seg1[:,1]<np.uint64(h/30), :] = 0
                        seg1[seg1[:,0]>h-np.uint64(h/30), :] = 0
                        seg1[seg1[:,1]>h-np.uint64(h/30), :] = 0
                    elif len(cropcoords) == 4:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<cropcoords[2], :] = 0
                        seg1[seg1[:,1]<cropcoords[0], :] = 0
                        seg1[seg1[:,0]>cropcoords[3], :] = 0
                        seg1[seg1[:,1]>cropcoords[1], :] = 0
                    else:
                        # removes punctae that are too close to the edge
                        seg1[seg1[:,0]<cropcoords[1]-cropcoords[2], :] = 0
                        seg1[seg1[:,1]<cropcoords[0]-cropcoords[2], :] = 0
                        seg1[seg1[:,0]>cropcoords[1]+cropcoords[2]+1, :] = 0
                        seg1[seg1[:,1]>cropcoords[0]+cropcoords[2]+1, :] = 0
                    
                    seg1 = seg1[seg1[:,0]!=0]
                    
                    # preallocates punctae centroid estimates
                    seg2 = np.zeros((len(seg1),2))
                    
                    # centers the image based on the brightest spot 
                    for i in range(len(seg1)):
                        
                        # find max and min windows
                        xmin = np.uint16(seg1[i,0])-5
                        xmax = np.uint16(seg1[i,0])+5
                        ymin = np.uint16(seg1[i,1])-5
                        ymax = np.uint16(seg1[i,1])+5
                        
                        # finds the indices of the brightest puncta
                        maxValue = np.max(image[xmin:xmax,ymin:ymax])
                        maxR, maxC = np.where(image[xmin:xmax,ymin:ymax] == maxValue)
                    
                        # adjusts puncta coords accordingly
                        seg2[i, 0] = np.float64(seg1[i, 0] + maxR[0] - 5)
                        seg2[i, 1] = np.float64(seg1[i, 1] + maxC[0] - 5)
                
                    # removes nonzero elements (from punctae that were too close to edge)
                    seg2 = seg2[np.matrix.nonzero(seg2[:,1])]
                    
                    # eventual centroid output matrix
                    centroids = seg2
                    
                    # preallocates vectors
                    punctaAxis  = np.zeros((2,len(seg2)))
                    punctaRemover = []
                    punctaArray = np.zeros((2*np.uint16(wndw)+1,2*np.uint16(wndw)+1,len(seg2)))
                    imageArray = np.zeros((2*np.uint16(wndw)+3,2*np.uint16(wndw)+3,len(seg2)))
                    backInt = np.zeros((1,len(seg2)))
                    punctInt = np.zeros((1,len(seg2)))
                    orient = np.zeros((len(seg2)))
                    pcovArray = np.zeros((len(seg2), 7))
                    
                    if suppression == "Suppress":
                        pass
                    else:
                        plt.title('Segmentation')
                        plt.imshow(image)
                        #tp.annotate(seg, image)
                        plt.scatter(seg2[:,1],seg2[:,0], facecolors='none', edgecolors='r')
                        plt.pause(0.01)
                    
                    # defines variables in case the interactive mode is not used
                    image = image*(image>back_i)
                    c = np.int16((2*wndw+1)/2)
                    x = np.linspace(0,2*wndw,2*wndw+1)
                    y = np.linspace(0,2*wndw,2*wndw+1)
                    x, y = np.meshgrid(x,y)
                    sphere = np.sqrt((wndw+3)**2-(x-c)**2-(y-c)**2)
                
                    currentdict = []
                    # curve fitting section and curve analysis
                    for i in range(len(seg2)):
                        
                        # describes the small window
                        xmin = np.uint64(seg2[i,0]-wndw)
                        xmax = np.uint64(seg2[i,0]+wndw+1)
                        ymin = np.uint64(seg2[i,1]-wndw)
                        ymax = np.uint64(seg2[i,1]+wndw+1)
                        
                        # finds a small area around puncta to fit a distribution
                        punct = image[xmin:xmax,ymin:ymax]
                        
                        # describes the small window
                        xminl = np.uint64(seg2[i,0]-wndw-1)
                        xmaxl = np.uint64(seg2[i,0]+wndw+2)
                        yminl = np.uint64(seg2[i,1]-wndw-1)
                        ymaxl = np.uint64(seg2[i,1]+wndw+2)
                        
                        imageArray[:,:,i] = rawImage[xminl:xmaxl,yminl:ymaxl]
                        
                        punctInt[0,i] = np.max(punct.ravel())
                        
                        punctBack = punct*(sphere<np.sqrt((wndw+3)**2-(p_size-1)**2))
                        
                        # finds the mitrochondrial background avg intensity and isolates puncta
                        backI = (sum(sum(punctBack)))/len(np.nonzero(punctBack)[0])
                        punct = np.float64(punct) - (backI)
                        punct[punct<0] = 0
                        
                        backInt[0,i] = backI
                        
                        # assigns the puncta to a multidimensional matrix for future reference
                        punctaArray[:,:,i] = punct
                        
                        # finds shape of puncta
                        h, w = punct.shape
                        
                        # creates a meshgrid of values as x and y values for 2d gaussian 
                        x = np.linspace(0,w-1,w)
                        y = np.linspace(0,h-1,h)
                        x, y = np.meshgrid(x,y)
                        
                        #necessary unravelling to work with scikit's curve fit
                        puncta = punct.ravel()
                        
                        try:
                        
                            # finds the gaussian 
                            with warnings.catch_warnings():
                                warnings.simplefilter('ignore')
                                popt, pcov = sci.optimize.curve_fit(gaussian2D, (x, y), puncta,
                                                                    p0 = [1500,5,5,5,0,5,0]) # make interactive
                            
                            pcovArray[i,:] = pcov.diagonal()**0.5
                            
                            if (pcovArray[i,1]>0.1) or (pcovArray[i,2]>0.1):
                                print(['Puncta ' + str(i) + ' file' + str(j+1)+ ' channel ' 
                                       + str(k+1)+ ' timestep ' + str(l+1)+ ' stack ' + str(m+1) +
                                       ' is being removed due to issues with its covariance'])
                                
                            else:
                                punctaRemover.append(i)
                            
                            # Find coefficents to solve for major and minor axis      
                            A = popt[3]/(np.log(2))
                            B = popt[4]/(np.log(2))
                            C = popt[5]/(np.log(2))
                            
                            
                            # find eigenvalues of ellipse
                            eigenMinor = (A+C)/2 + np.sqrt(((A+C)**2/4)-A*C+B**2/4)
                            eigenMajor = (A+C)/2 - np.sqrt(((A+C)**2/4)-A*C+B**2/4)
                            
                            if popt[3]>popt[5]:
                                orient[i] = pi/2 - np.arcsin(-popt[4]/(eigenMinor - eigenMajor))
                            elif popt[3]<=popt[5]:
                                orient[i] = np.arcsin(-popt[4]/(eigenMinor - eigenMajor))
              
                            # find major and minor axis
                            major = 2/np.sqrt(eigenMajor)
                            minor = 2/np.sqrt(eigenMinor)
                    
                            # places the major and minor axis in the matrix
                            punctaAxis[0,i]  = np.float64(major)
                            punctaAxis[1,i]  = np.float64(minor)
                            
                            #finds centroids [1,i] - x coord, [0,i] - y coord
                            centroids[i,1] = centroids[i,1] + popt[1] - 5
                            centroids[i,0] = centroids[i,0] + popt[2] - 5
                            
                        except RuntimeError:
                            print(['Puncta ' + str(i) + ' file' + str(j)+ ' channel ' 
                                   + str(k)+ ' timestep ' + str(l)+ ' stack ' + str(m) +
                                   ' is being removed due to issues with optimization'])
                            punctaAxis[0,i]  = np.nan
                            punctaAxis[1,i]  = np.nan
                            centroids[i,1] = np.nan
                            centroids[i,0] = np.nan
                            backInt[0,i] = np.nan
                            punctInt[0,i] = np.nan
                            
                    
                    
                    # transposes the major and minor axes to be more readable
                    pixelAxis = np.transpose(punctaAxis) 
                    punctaAxis = np.transpose(punctaAxis) * p_to_mu
                    
                    # calculates the aspect ratio
                    aspectRatio = pixelAxis[:,0]/pixelAxis[:,1]
                    partCoef = punctInt/backInt
                      
                    for i in punctaRemover:
                        d = {'x-centroids': centroids[i,1], 'y-centroids': centroids[i,0], 
                             'pixel major': pixelAxis[i,0], 'pixel minor': pixelAxis[i,1],
                             'puncta major': punctaAxis[i,0], 'puncta minor': punctaAxis[i,1],
                             'aspect ratio': aspectRatio[i], 'orientation': orient[i],
                             'average background': backInt[0,i],
                             'puncta intensity': punctInt[0,i], 'partition coeff': partCoef[0,i]}
                        currentdict.append(d)
                    df = pd.DataFrame(currentdict)
                    
                    for i in punctaRemover:
                        d2 = {'x-centroids': centroids[i,1], 'y-centroids': centroids[i,0], 
                              'pixel major': pixelAxis[i,0], 'pixel minor': pixelAxis[i,1],
                              'puncta major': punctaAxis[i,0], 'puncta minor': punctaAxis[i,1],
                              'aspect ratio': aspectRatio[i], 'orientation': orient[i],
                              'average background': backInt[0,i],
                              'puncta intensity': punctInt[0,i], 'partition coeff': partCoef[0,i],
                              'file number': 'file ' + str(j+1),'channel number': 'channel ' + str(k+1),
                              'frame number': l+1, 'stack number': 'zstack ' + str(m+1),
                              'puncta index': str(i), 'covariances': pcovArray[i,:]}
                        totaldicts.append(d2)
                    
                    for i in range(len(centroids[:,1])):
                        d3 = {'x-centroids': centroids[i,1], 'y-centroids': centroids[i,0], 
                              'pixel major': pixelAxis[i,0], 'pixel minor': pixelAxis[i,1],
                              'puncta major': punctaAxis[i,0], 'puncta minor': punctaAxis[i,1],
                              'aspect ratio': aspectRatio[i], 'orientation': orient[i],
                              'average background': backInt[0,i],
                              'puncta intensity': punctInt[0,i], 'partition coeff': partCoef[0,i],
                              'file number': 'file ' + str(j+1),'channel number': 'channel ' + str(k+1),
                              'frame number': l, 'stack number': 'zstack ' + str(m+1),
                              'puncta index': str(i), 'covariances': pcovArray[i,:]}
                        rawtotaldicts.append(d3)

                    # creates a nested dictionary that has all data
                    # files -> channels -> frames -> stacks -> dataframe
                    stacklist['zstack ' + str(m+1)] = df
                    arraylistz['zstack ' + str(m+1)] = punctaArray
                    imagelistz['zstack ' + str(m+1)] = imageArray
                framelist['frame ' + str(l+1)] = stacklist
                arraylistf['frame ' + str(l+1)] = arraylistz
                imagelistf['frame ' + str(l+1)] = imagelistz
            channellist['channel ' + str(k+1)] = framelist
            arraylistc['channel ' + str(k+1)] = arraylistf
            imagelistc['channel ' + str(k+1)] = imagelistf
        totallist['file ' + str(j+1)] = channellist
        arraylistt['file ' + str(j+1)] = arraylistc
        imagelistt ['file ' + str(j+1)] = imagelistc
    dframe = pd.DataFrame(totaldicts)
    rawdframe = pd.DataFrame(rawtotaldicts)
    print('\nDone')
    return totallist, arraylistt, imagelistt, dframe, rawdframe,

def unpacker(dictionary, file, channels, frames, stacks):
    
    # Unpacks the data returned by the function condensate.
    
    # Inputs 
    # totallist     - [var] The dictionary returned by the function condensate
    # file          - [int] The corresponding file number in the folder
    # channels      - [int] The desired channel number that must be pulled out
    # frames        - [int] The desired frame number that must be pulled out
    # stacks        - [int] The desired zstack number that must be pulled out
    
    # output        - [dataframe] The dataframe corresponding to parameters

    return dictionary['file ' + str(file)]['channel ' + str(channels)]['frame ' + str(frames)]['zstack ' + str(stacks)]

def fileWriter(dataframe, directory):
    dataframe.to_csv(directory)
    return directory

