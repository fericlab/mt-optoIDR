import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper\Figures\Fig2\Fusion\Line_profiles\Line_profile_data.csv')

frames = ['f3', 'f4', 'f5', 'f6', 'f7', 'f8']

max_int_f3 = df['f3'].max()

for i in frames:
    intensities = df[i]
    len(intensities)
    
    x=df.iloc[0:len(intensities), 0]
    
    norm_intensities = intensities/(max_int_f3)
    
    plt.plot(x, norm_intensities)
    plt.xlabel('Distance (um)')
    plt.ylabel('Normalized intensity (a.u)')
    plt.title(f'Frame {i}')
    
    plt.ylim(0, 1.7)
    
    plt.savefig(rf'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper\Figures\Fig2\Fusion\Line_profiles\Line_profile_{i}.svg')
    plt.show()
    
    


