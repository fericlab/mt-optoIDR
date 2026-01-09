#plotting distributions as a box-dot plot

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42  

font_properties = {
    'family': 'sans-serif',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}




# path to your csv files
files = glob.glob(r"C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper#1\Commun. Biol revisions\Revision1\FIGURES\Fig2\DifferentIDRs\Plots/*.csv")

palette = {
    "puncta major": "#8C89CE",
    "puncta minor": "#A3A3A3"
}

palette_2 = {
    "puncta major": "#6A66AD",
    "puncta minor": "gray"
}


all_data = []


for f in files:
    df = pd.read_csv(f)

    # keep only the two columns
    df = df[["puncta major", "puncta minor"]]

    # long format
    df_long = df.melt(
        var_name="Type",        # major / minor
        value_name="Value"
    )

    # label by file
    df_long["Dataset"] = os.path.splitext(os.path.basename(f))[0]

    all_data.append(df_long)

data = pd.concat(all_data, ignore_index=True)

plt.figure(figsize=(7, 4))
sns.violinplot(
    data=data,
    x="Dataset",
    y="Value",
    hue="Type",
    linecolor="black",
    palette=palette,
    width=0.8,
    linewidth=0.5,
    cut=0,
    zorder=0
)


# Overlay all data points
sns.stripplot(
    data=data,
    x="Dataset",
    y="Value",
    hue="Type",
    palette=palette_2,
    dodge=True,
    jitter=0.1,
    alpha=0.2,
    size=4,
    linewidth=0,
    zorder=0
)








# Fix duplicate legend
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
n = data["Type"].nunique()
plt.legend(handles[:n], labels[:n], title="Type", fontsize=14, frameon=False)



plt.tick_params(axis='both', labelsize=14)
plt.ylabel('Size(µm)', **font_properties)
plt.xlabel('IDR', **font_properties)
plt.ylim(0.1,0.55)
plt.tight_layout()
plt.savefig('Violine_plot.pdf', dpi=600)
plt.show()