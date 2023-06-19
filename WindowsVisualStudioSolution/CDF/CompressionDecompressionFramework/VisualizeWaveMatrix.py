import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the CSV data into a pandas DataFrame
data = pd.read_csv('wave_data.csv', names=['x', 'y', 'z', 'value'])

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
scatter = ax.scatter(data['x'], data['y'], data['z'], c=data['value'], cmap='coolwarm')

# Label the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Add a colorbar to the scatter plot
fig.colorbar(scatter)

# Show the plot
plt.show()
