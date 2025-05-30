# mt_optoIDR_scripts
## mt_condensate_finder 1.0
Main scripts used to segment mt-opto-condensates from the mitochondrial network. 
This script calculates,   
* Major and minor axis of segmented droplets (full width at half max)   
* Aspect ratio of droplets   
* Partition coefficients    
and includes all the x and y location details and other data calculated from the trackpy locate function 


## Orientation finder
* This script calculates the orientation between mitochondrial axial axis and droplet major axis

## mt_condensate_tracker 1.0
This script contains  
* tracking functions,  
* a function to isolate individual tracks,   
* a rotate track function,   
* a function to calculate van-Hove plots for a given dimension,   
* a function to calculate MSD for a given dimension.   
