import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

# Load the CSV data into a pandas DataFrame
data = pd.read_csv('wave_data.csv', names=['x', 'y', 'z', 'value'])

# Get unique z levels
z_levels = data['z'].unique()

# Define the grid where the data will be interpolated
grid_x, grid_y = np.mgrid[min(data['x']):max(data['x']):100j, min(data['y']):max(data['y']):100j]

# Create a single figure
fig = plt.figure()

# Determine the number of rows and columns for the subplot grid
rows = int(np.ceil(np.sqrt(len(z_levels))))
cols = int(np.ceil(len(z_levels) / rows))

# Create a 3D surface plot for each z level
for i, z in enumerate(z_levels):
    z_data = data[data['z'] == z]

    # Interpolate the values for a smoother plot
    grid_z = griddata((z_data['x'], z_data['y']), z_data['value'], (grid_x, grid_y), method='cubic')

    ax = fig.add_subplot(rows, cols, i+1, projection='3d')
    surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap='coolwarm')

    # Label the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Value')
    ax.set_title(f'Z = {z}')

    # Add a colorbar to the subplot
    fig.colorbar(surf, ax=ax)

# Display the figure
plt.show()
