#bar and data plot
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob, os
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42  

font_properties = {
    'family': 'sans-serif',
    'size': 15,
    'weight': 'normal',
    'style': 'normal'
}



# -------- SETTINGS --------
data_path = r"C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper#1\Commun. Biol revisions\Revision1\FIGURES\Fig2\DifferentIDRs\Plots\*.csv"

palette = {
    "partition coeff": "#8C89CE"
}

palette_2 = {
    "partition coeff": "#564F8E"
}


# -------- LOAD & FORMAT DATA --------
files = glob.glob(data_path)
all_data = []

for f in files:
    df = pd.read_csv(f)
    df = df[["partition coeff"]]

    df_long = df.melt(var_name="Type", value_name="Value")
    df_long["Dataset"] = os.path.splitext(os.path.basename(f))[0]

    all_data.append(df_long)

data = pd.concat(all_data, ignore_index=True)

# -------- PLOT --------


plt.figure(figsize=(7, 4))
sns.violinplot(
    data=data,
    x="Dataset",
    y="Value",
    hue="Type",
    linecolor="black",
    palette=palette,
    width=0.4,
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
plt.ylabel('Partition coeff', **font_properties)
plt.xlabel('IDR', **font_properties)
plt.tight_layout()
plt.savefig('Partition_coeff_violine_plot.pdf', dpi=600)
plt.show()