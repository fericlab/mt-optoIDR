import pandas as pd
import matplotlib.pyplot as plt
import plotly.io as pio # todisplay plots in browser 
import plotly.graph_objects as go
pio.renderers.default='browser'
from plotly.offline import plot

df=pd.read_csv(r'C:/Users/sanja/OneDrive - The Pennsylvania State University/Research/paper#1/New_Figures/S3/Phase_daigram/data.csv')
df.head()


colors = df['Droplets']

x=df['I Max Airyscan mCherry (Preactivation)']
x_norm=x/x.max()
y=df['I Max Airyscan HaloSiR (Preactivation)']
y_norm=y/y.max()

plt.figure(figsize=(6, 6))  # 8 inches wide, 6 inches tall
plt.scatter(x=x_norm,y=y_norm , c=colors, alpha=0.8, cmap='winter')
plt.tick_params(axis='both', which='major', labelsize=14)


plt.xlabel('Normalized intensity of IDR-mCherry-CRY2olig (a.u)')
plt.ylabel('Normalized intensity of Halo-MTS (SiR HaloX) (a.u)')

plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('phase_daigram_norm.pdf')
plt.show()



#plt.scatter(x=x,y=y , c=colors, alpha=0.8, cmap='winter')
#plt.savefig('phase_daigram.svg')
#plt.show()



############Plotly html graphs###########
fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, opacity=0.8, mode='markers',
                         marker=dict(size=15, 
                                     color=df['Droplets'],  
                                     colorscale='Viridis',  
                                     showscale=False  
    )))
fig.update_layout(
    template='none',
    title='Phase diagram',
    xaxis_title='mCh intensity (a.u)',
    yaxis_title='HaloX Sir intensity (a.u)'
)


# Save the plot as an HTML file
#plot(fig, filename= r'C:\Users\sanja\OneDrive - The Pennsylvania State University\Research\paper\Figures\S2\Phase_daigram\Phase_digram.html', auto_open=False) 
fig.show()