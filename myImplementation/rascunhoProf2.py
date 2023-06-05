import matplotlib.pyplot as plt
import numpy as np
from shapely.wkt import loads
from shapely.geometry import Polygon
from shapely.affinity import translate, scale, affine_transform
import math
import geopandas as gpd

import ourSelectedFunctions as our

# ### INPUT data: method 1
# Give the coordinates of the source and target
pointsSrc = np.array([[2, 3], [6, 3], [5, 0.5], [3, 0.5]])
pointsTrg = np.array([[3, 1], [5, 2], [5, -3], [3, -3]])
polySrc = Polygon(pointsSrc)
polyTrg = Polygon(pointsTrg)

# ### Ploting source, target and interpolation
xmin = ymin = -6
xmax = ymax = 6

fig1, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)

# plots the source polygon
polyGPDSimp = gpd.GeoSeries(polySrc)
ax1.set_title('source')
ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
polyGPDSimp.plot(ax=ax1, color="blue", alpha=0.5)

# plots the target polygon
polyGPDSimp = gpd.GeoSeries(polyTrg)
ax2.set_title('target')
ax2.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
polyGPDSimp.plot(ax=ax2, color="red", alpha=0.5)

# plots the source and target polygons aligned. 
# This is the case you need to compute the distances to find the best correspondences.
transformed_polygon = our.interpolate_fixed(polySrc, polyTrg, t=1)
polyGPDSimp = gpd.GeoSeries(transformed_polygon)
ax3.set_title('Source aligned in target position')
ax3.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
#ax3.set(xlim=(0, 6), ylim=(-6, 0))
polyGPDSimp.plot(ax=ax3, color="green", alpha=0.5)

# plots the estimatedpolygon at time t
transformed_polygon = our.interpolate_deformable(polySrc, polyTrg, t=0.5)
polyGPDSimp = gpd.GeoSeries(transformed_polygon)
ax4.set_title('Estimated target')
ax4.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
#ax3.set(xlim=(0, 6), ylim=(-6, 0))
polyGPDSimp.plot(ax=ax4, color="green", alpha=0.5)

plt.show()





fig1, ax5 = plt.subplots()

for frame in range(25):
    ax5.clear()
    # Rotate the triangle by an angle at each frame
    transformed_polygon = our.interpolate_deformable(polySrc, polyTrg, t=frame / 24)
    polyGPDSimp = gpd.GeoSeries(transformed_polygon)
    ax5.set_title('Estimated target')
    ax5.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    #ax3.set(xlim=(0, 6), ylim=(-6, 0))
    polyGPDSimp.plot(ax=ax5, color="green", alpha=0.5)
    
    plt.pause(0.1)
plt.show()