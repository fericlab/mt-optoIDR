#Box or violin plot
import matplotlib.pyplot as plt
import pandas as pd


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



data_WT = pd.read_csv('C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/Data_Analysis/Size_comparison_analysis/data_20241026_WT_original.csv')
data_Drp1k38a = pd.read_csv('C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/Data_Analysis/Size_comparison_analysis/data_Drp1k38a_sizes.csv')


WT_major_raw = data_WT['puncta major']
WT_minor_raw = data_WT['puncta minor']

WT_major = data_WT.groupby('file number')['puncta major'].mean()
WT_minor = data_WT.groupby('file number')['puncta minor'].mean()

Drp1_major = data_Drp1k38a.groupby('file number')['puncta major'].mean()
Drp1_minor = data_Drp1k38a.groupby('file number')['puncta minor'].mean()


categories = ['WTmajor', 'WTminor', 'Drp1K38Amajor', 'Drp1K38Aminor']


########### violin plots
data_list = [WT_major_raw, WT_minor_raw, Drp1_major, Drp1_minor] #without averaging for each cell in WT

#data_list = [WT_major, WT_minor, Drp1_major, Drp1_minor] #with averaging for each cell in WT


vp = plt.violinplot(dataset=data_list, showmeans=True, showmedians=False, showextrema=True)
plt.ylabel(" Diameter (µm)", **font_properties)
plt.xticks(ticks=[1, 2, 3, 4], labels=['WT_major', 'WT_minor', 'Drp1_major', 'Drp1_monor',])
plt.tick_params(axis='x', labelsize=tick_font_size)  # Font size for x-axis tick labels
plt.tick_params(axis='y', labelsize=tick_font_size)

colors = ['#564F8E', '#564F8E', 'gray', 'gray']

for i, body in enumerate(vp['bodies']):  
    body.set_facecolor(colors[i])  # Fill color
    body.set_edgecolor('black')    # Outline color
    body.set_alpha(0.7)            # Transparency

# Customize median and mean lines
vp['cmeans'].set_color('black')  # Mean line color
#vp['cmedians'].set_color('black') # Median line color
#vp['cmedians'].set_linewidth(2)   # Median line thickness

vp['cbars'].set_color('black')  # Change whisker line color
vp['cbars'].set_linewidth(1)  # Adjust thickness

# Change the top and bottom caps (if needed)
vp['cmins'].set_color('black')  # Bottom cap
vp['cmaxes'].set_color('black')  # Top cap



plt.xticks(rotation=60)
plt.tight_layout()
plt.savefig('size_comparison_violin_plot.pdf')
plt.show()


print('######################################')
print(f'WT_major mean:{WT_major.mean():.4f}, std:{WT_major.std():.4f}' )
print(f'WT_minor mean: {WT_minor.mean():.4f}, std:{WT_minor.std():.4f}')

print(f'Drp1_major mean: {Drp1_major.mean():.4f}, std:{Drp1_major.std():.4f}' )
print(f'Drp1_minor mean: {Drp1_minor.mean():.4f}, std: {Drp1_minor.std():.4f}')




