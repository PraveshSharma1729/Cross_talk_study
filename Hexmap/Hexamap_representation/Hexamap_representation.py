import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('/home/pravesh/Desktop/pad_ch_positions.txt', sep='\t') #reading pad position to data frame

# Scatter plot with hexagonal shapes and colormap
plt.figure(figsize=(10, 10))
scatter = plt.scatter(
    df['x'], df['y'], 
    c=df['Econde_Rx'],  # Colors determined by Column 3
    marker='h',  # Hexagonal marker
    cmap='viridis',  # Choose a colormap
    s=1000,  # Marker size
    alpha=0.9  # Transparency
)

# Add labels to each scatter point from Column 0
for i in range(len(df)):
    plt.text(
        df['x'][i],  # X position
        df['y'][i],  # Y position
        df['N_channel'][i],  # Text (value from Column 0)
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
plt.title('Scatter Plot with Hexagonal Markers and Colormap', fontsize=20)
plt.savefig("Hexmap.jpg")
plt.show()


