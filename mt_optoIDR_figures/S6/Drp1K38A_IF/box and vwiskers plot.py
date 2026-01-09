import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns

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

###################################


# Import data
wt = pd.read_csv('intensity_data_0.csv')
drp1 = pd.read_csv('intensity_data_1.csv')



# Combine data into one DataFrame
wt['group'] = 'WT'
drp1['group'] = 'Drp1K38A'


wt_m = wt.groupby(['file number', 'group']).mean().reset_index()
drp1_m = drp1.groupby(['file number', 'group']).mean().reset_index()


combined = pd.concat([wt_m, drp1_m], ignore_index=True)

# Set styles
box_colors = {'WT': '#5F2C96', 'Drp1K38A': '#131313'}     # Box fill colors
dot_colors = {'WT': '#564F8E', 'Drp1K38A': '#000000'}     # Dot colors

plt.figure(figsize=(5, 4))

# Draw boxplot with transparent fill
sns.boxplot(
    x='group', y='average network int', data=combined,
    palette=box_colors,
    width=0.3,
    linewidth=1.5,
    fliersize=0,
    boxprops=dict(edgecolor='black', linewidth=1.5, alpha=0.7),
    medianprops=dict(color='black', linewidth=1.2),
    whiskerprops=dict(color='black', linewidth=1.2),
    capprops=dict(color='black', linewidth=1.2)
)


# Overlay individual points, colored by group
for group in combined['group'].unique():
    subset = combined[combined['group'] == group]
    sns.stripplot(
        x='group', y='average network int', data=subset,
        jitter=0.1, size=5, alpha=0.8,
        color=dot_colors[group]
    )



wt_means = wt.groupby('file number')['average network int'].mean()
wt_mean = wt_means.mean()
print('Wild type network average intensity:', wt_mean)


drp1_means = drp1.groupby('file number')['average network int'].mean()
drp1_mean = drp1_means.mean()
print('Drp1K38A network average intensity:', drp1_mean)


plt.annotate(f'wt_mean={wt_mean:.2f},\ndrp1_mean={drp1_mean:.2f}', xy=(0.1,0.5), xycoords='axes fraction')
# Styling
plt.title('Box wiskers plot with Individual Data Points')
plt.ylabel('Average Intensity (a.u.)', **font_properties)
# KEEP box around the plot
sns.despine(trim=False, top=False, right=False)
plt.tight_layout()
plt.savefig('network_int_box_wiskers_plot.pdf')
plt.show()





###statistical test

data_wt = wt_m['average network int']
data_drp1 = drp1_m['average network int']

import scipy.stats as stats
t_stat, p_value = stats.ttest_ind(data_wt, data_drp1)

# Print the results
print(f"T-statistic: {t_stat}")
print(f"P-value: {p_value}")

plt.hist(data_wt, bins=5)



# Perform Mann-Whitney U test
statistic, p_value = stats.mannwhitneyu(data_wt, data_drp1)
print(f"Mann-Whitney U test: U = {statistic:.2f}, p = {p_value:.4f}")







