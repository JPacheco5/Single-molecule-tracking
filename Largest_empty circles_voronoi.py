#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 19:23:39 2024

@author: jep160
"""

from PIL import Image
import numpy as np
from scipy.spatial import Delaunay, distance, Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import pandas as pd

def circumcenter(vertex1, vertex2, vertex3):
    """
    Calculate the circumcenter of a triangle defined by three vertices.

    Parameters:
        vertex1, vertex2, vertex3: Tuple or array-like, representing the coordinates of three vertices.

    Returns:
        Tuple, coordinates of the circumcenter.
    """
    x1, y1 = vertex1
    x2, y2 = vertex2
    x3, y3 = vertex3

    D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

    Ux = ((x1**2 + y1**2) * (y2 - y3) + (x2**2 + y2**2) * (y3 - y1) + (x3**2 + y3**2) * (y1 - y2)) / D
    Uy = ((x1**2 + y1**2) * (x3 - x2) + (x2**2 + y2**2) * (x1 - x3) + (x3**2 + y3**2) * (x2 - x1)) / D

    return Ux, Uy

# Load binary image
img = Image.open('/your directory/mask.tif').convert('L')
z = np.array(img, copy=True)
z[z == 255] = 1

# Add a frame of 1-pixel width along all edges
z_padded = np.pad(z, 1, mode='constant', constant_values=1)

# Extract coordinates of points (where binary_image == 1)
points_y, points_x = np.where(z_padded == 1)
points = np.column_stack((points_x, points_y))

# Compute Delaunay triangulation
tri = Delaunay(points)

# Compute Voronoi diagram
vor = Voronoi(points)

# Calculate circumcenter distances and radii
circumcenter_distances = []
radii_data = {'Circumcenter': [], 'Radius': []}

for i, simplex in enumerate(tri.simplices):
    circumcenter_coords = circumcenter(*tri.points[simplex])
    circumcenter_coords = np.round(circumcenter_coords).astype(int)
    
    # Calculate distances from the circumcenter to all vertices
    vertex_distances = [distance.euclidean(circumcenter_coords, vertex) for vertex in tri.points[simplex]]

    circumcenter_distances.append(vertex_distances)

    # Find the minimum distance to any vertex as the radius
    min_distance = min(vertex_distances)
    radii_data['Circumcenter'].append(circumcenter_coords)
    radii_data['Radius'].append(min_distance)

# Convert the data to a pandas DataFrame
radii_df = pd.DataFrame(radii_data)

# Find the largest empty circle
largest_empty_circle = radii_df.loc[radii_df['Radius'].idxmax()]

# Plot the binary image
plt.imshow(z_padded, cmap='hot', origin='upper')

# Plot Delaunay triangulation
plt.triplot(points[:, 0], points[:, 1], tri.simplices, alpha=0.2, linewidth=0.1)

# Plot the Voronoi diagram
voronoi_plot_2d(vor, show_vertices=False, line_colors='r', line_width=0.5, alpha=0.1)

# Plot the largest empty circle
circle = plt.Circle(largest_empty_circle['Circumcenter'], largest_empty_circle['Radius'],
                    edgecolor='green', facecolor='none', linewidth=2)

plt.gca().add_patch(circle)

# Invert the y-axis after plotting the image
plt.gca().invert_yaxis()

# Set aspect ratio to 'equal'
plt.gca().set_aspect('equal', adjustable='box')

# Set limits on the x and y axes
plt.xlim(0, 100)  # Replace min_x_limit and max_x_limit with your desired values
plt.ylim(100, 0)  # Replace min_y_limit and max_y_limit with your desired values

# Save the figure with the desired DPI
plt.savefig('/your directory/output_plot.png', dpi=300, bbox_inches='tight')
plt.show()

# Save DataFrame to CSV
radii_df['Largest_empty_circle'] = largest_empty_circle[1]
radii_df.to_csv('/your directory/Radius.csv')
