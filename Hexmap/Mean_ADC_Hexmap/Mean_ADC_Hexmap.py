import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('/home/pravesh/Desktop/pad_ch_positions.txt', sep='\t') #reading pad position to data frame
mean= pd.read_csv('/home/pravesh/Desktop/Cross-Talk/New folder/Run1695495152/Mean_ADC.txt', sep='\t') #input txt file having pedestal mean

# Scatter plot with hexagonal shapes and colormap
plt.figure(figsize=(13, 10))
scatter = plt.scatter(
    df['x'], df['y'], 
    c=mean['L2_Mean_ADC'],  # Colors determined by Column 3
    marker='h',  # Hexagonal marker
    cmap='jet',  # Choose a colormap
    s=1000,  # Marker size
    alpha=1  # Transparency
)

# Add labels to each scatter point from Column 0
for i in range(len(df)):
    plt.text(
        df['x'][i],  # X position
        df['y'][i],  # Y position
        int(mean['L2_Mean_ADC'][i]),  # Text (value from Column 0)
        fontsize=10,  # Font size for the labels
        color='white',  # Color of the text
        ha='center',  # Horizontal alignment
        va='center'  # Vertical alignment
    )

# Remove ticks and axis lines
plt.xticks([])  # Remove x-axis ticks
plt.yticks([])  # Remove y-axis ticks
plt.gca().spines['top'].set_visible(False)  # Hide the top axis line
plt.gca().spines['right'].set_visible(False)  # Hide the right axis line
plt.gca().spines['left'].set_visible(False)  # Hide the left axis line
plt.gca().spines['bottom'].set_visible(False)  # Hide the bottom axis line

# Add labels and title
plt.title('Run1695495152 (Pedestal) \nLayer 2 Mean ADC', fontsize=20)
# Add color bar to show the mapping of colors
cbar = plt.colorbar(scatter, label='Mean ADC Values')
cbar.set_label('Mean ADC Values', fontsize=18)
cbar.ax.tick_params(labelsize=15) 
plt.savefig("Pedestal_mean_ADC.jpg")
plt.show()
