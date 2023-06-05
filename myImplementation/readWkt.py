import matplotlib.pyplot as plt
import numpy as np
from shapely.wkt import loads
from shapely.geometry import Polygon
from shapely.affinity import translate, scale, affine_transform
import math

def show_polygon_from_wkt(wkt_file_path):
    # Read the WKT file and create a polygon object
    with open(wkt_file_path, 'r') as f:
        wkt = f.read()
    polygon = loads(wkt)

    # Extract the coordinates from the polygon
    x, y = polygon.exterior.xy

    # Create a matplotlib figure and axis
    fig, ax = plt.subplots()

    # Plot the polygon
    ax.plot(x, y)

    # Set labels and title
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Polygon from WKT')

def show_polygon_from_coords(points):    
    polygon = Polygon(points)
    
    # Extract the exterior coordinates from the Polygon
    x, y = polygon.exterior.xy
    
    # Plot the polygon
    plt.figure()
    plt.plot(x, y)
    plt.fill(x, y, alpha=0.5)
    plt.title('Rotated Polygon')

def register_point_set(source_wkt, target_wkt): #funcao para alinhar os poligonos caso haja rotacao
    # Parse the WKT strings to create Shapely polygon objects
    with open(source_wkt, 'r') as f:
        wkt = f.read()
    source_polygon = loads(wkt)

    with open(target_wkt, 'r') as f:
        wkt = f.read()
    target_polygon = loads(wkt)
    
    #centroide dos pols
    source_centroid = source_polygon.centroid
    target_centroid = target_polygon.centroid

    #vetor translacao para alinhar centroides (subtracao do prof)
    translation_vector = (target_centroid.x - source_centroid.x,
                          target_centroid.y - source_centroid.y)

    #translacao do target (funcao do shapely)
    translated_target_polygon = affine_transform(target_polygon, [1, 0, 0, 1, translation_vector[0], translation_vector[1]])

    #calculo do angulo de rotacao
    source_angle = math.atan2(source_polygon.exterior.coords[0][1] - source_centroid.y,
                              source_polygon.exterior.coords[0][0] - source_centroid.x)
    target_angle = math.atan2(translated_target_polygon.exterior.coords[0][1] - target_centroid.y,
                              translated_target_polygon.exterior.coords[0][0] - target_centroid.x)
    rotation_angle = source_angle - target_angle
    print(rotation_angle)
    if math.isclose(rotation_angle, 0.0, abs_tol=1): #se os angulos nao forem mais do que 20 graus de diferenca nao roda
        print("no rotation")
        return translated_target_polygon
    else: #caso contrario roda
        #faz a rotacao do target (apos translacao aplicada)
        print("rotation")
        aligned_target_polygon = affine_transform(translated_target_polygon, [math.cos(rotation_angle), -math.sin(rotation_angle), math.sin(rotation_angle), math.cos(rotation_angle), 0, 0])
        return aligned_target_polygon

show_polygon_from_wkt('g25.wkt')
show_polygon_from_wkt('g26.wkt')
points = register_point_set("g25.wkt","g26.wkt")
show_polygon_from_coords(points)
plt.show()