import pandas as pd
import matplotlib.pyplot as plt
import glob
import matplotlib.cm as cm
import matplotlib.colors as colors

#this makes text editable
plt.rcParams['pdf.fonttype'] = 42

data_path = r"C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper#1\Commun. Biol revisions\Revision1\FIGURES\Fig3\Shape_change\singleDroplet_line\*.csv"

files = glob.glob(data_path)

fig, ax = plt.subplots(figsize=(5,3))

norm = colors.Normalize(vmin=0, vmax=(len(files) - 1)*2) #normalize to seconds x2 as frames are 2 secounds apart 

cmap1 = cm.Greens_r
cmap2 = cm.RdPu_r

i=0
for f in files:
    
    df = pd.read_csv(f)
    
    # First column = distance
    x = df['Distance'] 
    y_cristae = df['Cristae'] /df['Cristae'].sum()
    y_mtopto = df['mt-optoIDR'] /df['mt-optoIDR'].sum()
    
    
    y_cristae=y_cristae/y_mtopto.max()
    y_mtopto=y_mtopto/y_mtopto.max()
    
    #mtOptoIDR peak 
    peak_idx = y_mtopto.idxmax()
    peak_x = x.loc[peak_idx]
    
    x_centered = x - peak_x
    
    
    #different color maps for cristea and mtoptoIDR
    
    color_mtopto = cmap1(norm(i))
    color_cristae = cmap2(norm(i))
    
    ax.plot(
        x_centered,
        y_mtopto,
        color=color_mtopto,
        alpha=0.8,
        linewidth=1.5
    )

    ax.plot(
        x_centered,
        y_cristae,
        color=color_cristae,
        alpha=0.8,
        linewidth=1.5,
        linestyle="--",
    )
    print(i)
    i=i+2
    

sm1 = cm.ScalarMappable(norm=norm, cmap=cmap1)
sm1.set_array([])

sm2 = cm.ScalarMappable(norm=norm, cmap=cmap2)
sm2.set_array([])


cbar = plt.colorbar(sm1, ax=ax)
cbar = plt.colorbar(sm2, ax=ax)
cbar.set_label("Time (s)")


#ax.axvline(0, color="k", linestyle="--", linewidth=1)  # peak
ax.set_xlabel("Distance (µm)")
ax.set_ylabel("Intensity (a.u.)")
plt.xlim(-0.45,0.45)
plt.tight_layout()
fig.savefig('int_sigle_drop.pdf', dpi=600)
plt.show()
    
    
    
    








