import matplotlib.pyplot as plt
import numpy as np
from shapely.wkt import loads
from shapely.geometry import Polygon
from shapely.affinity import translate, scale, affine_transform
import math
import geopandas as gpd

def register_point_set(sourceOrigin, t, angle, dx, dy):

    #csource = source.centroid
    #cx_target, cy_target = target.centroid
    #translation_vector = (-csource.x,
    #                      -csource.y)

    polyTransformed = affine_transform(sourceOrigin, [math.cos(angle * t/32), -(math.sin(angle * t/32)), math.sin(angle * t/32), math.cos(angle * t/32), dx * t / 32, dy * t / 32])
    
    return polyTransformed

pointsSrc = [(2,3),(6,3),(4,0)]
angle = math.pi / 4
polySrc = Polygon([[p[0], p[1]] for p in pointsSrc])
dx = polySrc.centroid.x
dy = polySrc.centroid.y

polySrcOrigin = affine_transform(polySrc, [1, 0, 0, 1, dx, dy])
polyTrgOrigin = affine_transform(polySrcOrigin, [math.cos(angle), -(math.sin(angle)), math.sin(angle), math.cos(angle), dx, dy])
#print(pointsTrg)
polySrc = Polygon([[p[0], p[1]] for p in pointsSrc])
#print(polySrc)

fig1, ax1 = plt.subplots()
for t in range(33):
    polyTranformed = register_point_set(polySrc, t, angle, dx, dy)
#    print(polyTranformed)
    x, y = polyTranformed.exterior.xy
    
    # Plot the polygon
    polyGPDSimp = gpd.GeoSeries(polyTranformed) #transforma em gpd
    polyGPDSimp.plot(ax=ax1, color="blue", alpha=0.5)
    #plt.figure()
    plt.title('Rotated Polygon')
plt.show()
#x_target = [4, 8, 6]
#y_target = [2, 2, -1]

