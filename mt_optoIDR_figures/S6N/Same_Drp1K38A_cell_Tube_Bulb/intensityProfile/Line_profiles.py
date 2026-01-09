import pandas as pd
import matplotlib.pyplot as plt
import tifffile as tiff
# Load data
df = pd.read_csv("tubeInt.csv")

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

#get the max int of the ffigure 
image = tiff.imread(r"C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper#1\Commun. Biol revisions\Revision1\Analysis\Same_cell_different_sizes\Figure\20241004_DDX28-mCh-Cry2olig_Drp1k38a_PostImaging_CellNo05-01_Airyscan_processed-1.tif")
chan2 = image[1,:,:]
plt.imshow(chan2)
plt.show()


# First column = distance
x = df.iloc[:, 0]

# Choose a colormap
cmap = plt.get_cmap("Accent")   # 'plasma', 'inferno', 'magma', 'cividis', 'tab10', etc.

# How many channels (excluding first column)
num_channels = len(df.columns[1:])

plt.figure(figsize=(4, 3))
# Loop with colors
for i, channel in enumerate(df.columns[1:]):
    intensities = df[channel]
    #Noemalized values
    norm_intensities = intensities /chan2.max()
    
    #Raw values
    #norm_intensities = intensities
    # Avoid ZeroDivisionError: if only one channel, just pick middle of colormap
    if num_channels > 1:
        color = cmap(i / (num_channels - 1))   # evenly spaced
    else:
        color = cmap(0.5)

    plt.plot(x, norm_intensities, label=channel, alpha=0.8, linewidth=3, color='#A1CE6B')
    #plt.scatter(x, norm_intensities, alpha=1, color=color)


fzise=16
# Labels and styling
plt.xlabel("Distance (µm)", fontsize=fzise)
plt.ylabel("Norm: Intensity (a.u.)", fontsize=fzise)
plt.tick_params(axis="both", which="major", labelsize=fzise)
plt.legend(frameon=False)
plt.tight_layout()
plt.ylim(bottom=0, top=1.2)
# Save figure
plt.savefig("tube.pdf", dpi=300)
plt.show()


