import pandas as pd
import matplotlib.pyplot as plt

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

# Load data
df = pd.read_csv("DataL0.csv")

# First column = distance
x = df.iloc[:, 0]

# Choose a colormap
#cmap = plt.get_cmap("Accent")   # 'plasma', 'inferno', 'magma', 'cividis', 'tab10', etc.
my_colors = ['#800080', '#A1CE6B']
# How many channels (excluding first column)
num_channels = len(df.columns[1:])

plt.figure(figsize=(5, 3))
# Loop with colors
for i, channel in enumerate(df.columns[1:]):
    
    intensities = df[channel]
    #Noemalized values
    if i==0:
        norm_intensities = intensities /966 #max of the cristae channel
    else:
        norm_intensities = intensities /1458 #max of the mt-optoidr channel
    print(intensities.max())


    
    #Raw values
    #norm_intensities = intensities
    # Avoid ZeroDivisionError: if only one channel, just pick middle of colormap
    
    
    current_color = my_colors[i % len(my_colors)]

    plt.plot(x, norm_intensities, label=channel, alpha=0.8, linewidth=2, color=current_color)
    #plt.scatter(x, norm_intensities, alpha=1, color=current_color)


# Labels and styling
plt.xlabel("Distance (µm)", fontsize=15)
plt.ylabel("Norm: Intensity (a.u.)", fontsize=15)
plt.tick_params(axis="both", which="major", labelsize=15)
plt.xlim(-0.04, 1.55)
plt.ylim(0.3,1.1)
plt.legend(frameon=False)
plt.tight_layout()
# Save figure
plt.savefig("intensity_profile.pdf", dpi=300)
plt.show()


