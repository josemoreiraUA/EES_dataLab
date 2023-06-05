import numpy as np
from scipy.spatial import procrustes
from shapely.geometry import Polygon, Point
from shapely.affinity import affine_transform
import geopandas as gpd

############# SVD
# Apply the same translation and rotation to all vertices of the source polygon

def svdAlignment(src, trg):
    # Compute the centroids of the polygons
    src_centroid = np.mean(src, axis=0)
    trg_centroid = np.mean(trg, axis=0)

    # Center the polygons by subtracting the centroids
    src_centered = src - src_centroid
    trg_centered = trg - trg_centroid

    # Compute the covariance matrix
    covariance_matrix = src_centered.T @ trg_centered

    # Perform Singular Value Decomposition
    U, W, Vt = np.linalg.svd(covariance_matrix)

    # Calculate the rotation matrix
    rotation_matrix = Vt.T @ U.T

    # Calculate the translation vector
    translation_vector = trg_centroid - src_centroid

    # Compute the distances between the vertices in the source and target aligned (skew / deformation)
    print("2.", src_centered)
    print("1.", trg_centered)
    print("3.", src_centered @ rotation_matrix)
    skew_matrix = trg_centered - (src_centered @ rotation_matrix)

    # The result
    return translation_vector, rotation_matrix, skew_matrix


def interpolate_fixed(src, trg, t):
    assert 0 <= t <= 1, 'parameter t \in [0, 1]'  
    #assert type(source) == 'shapely.geometry.polygon.Polygon'
    print("inside func: ", src.exterior.coords.xy)
    
    translation, rotation, _ = svdAlignment(np.array(src.exterior.coords.xy), np.array(trg.exterior.coords.xy))

    # Compute the centroids of the polygons
    src_centroid = np.mean(np.array(src.exterior.coords.xy), axis=0)
 
    # Center the source by subtracting the centroid
    src_centered = np.array(src.exterior.coords.xy) - src_centroid
    
    # Extract the cosine and sine values
    cos_theta = rotation[0, 0]
    sin_theta = rotation[1, 0]

    # Calculate the angle in radians
    theta = np.arctan2(sin_theta, cos_theta)

    # Convert the angle to degrees
    theta_deg_at_t = np.rad2deg(theta) * t

    # Convert angles from degrees to radians
    angle_rad = np.deg2rad(theta_deg_at_t) 
    
    rotation_matrice = np.array([ [np.cos(angle_rad), -np.sin(angle_rad)], [np.sin(angle_rad), np.cos(angle_rad)] ])

    #transformation = np.dot(np.array(source.exterior.coords),rotation_matrice.T) + translation * t
    transformation = np.dot(src_centered,rotation_matrice.T) + translation * t + src_centroid
    transformed_polygon = Polygon(transformation)    
    
    # alternative
    #transformed_polygon = affine_transform(source, [np.cos(angle_rad), -np.sin(angle_rad), np.sin(angle_rad), np.cos(angle_rad), translation[0] * t, translation[1] * t])

    return transformed_polygon

def interpolate_deformable(src, trg, t):
    assert 0 <= t <= 1, 'parameter t \in [0, 1]'  
    #assert type(source) == 'shapely.geometry.polygon.Polygon'

    translation, rotation, deformation = svdAlignment(np.array(src.exterior.coords), np.array(trg.exterior.coords))

    # Compute the centroids of the polygons
    src_centroid = np.mean(np.array(src.exterior.coords), axis=0)
 
    # Center the source by subtracting the centroid
    src_centered = np.array(src.exterior.coords) - src_centroid
    
    # Extract the cosine and sine values
    cos_theta = rotation[0, 0]
    sin_theta = rotation[1, 0]

    # Calculate the angle in radians
    theta = np.arctan2(sin_theta, cos_theta)

    # Convert the angle to degrees
    theta_deg_at_t = np.rad2deg(theta) * t
    print("> theta_deg", theta_deg_at_t) 

    # Convert angles from degrees to radians
    angle_rad = np.deg2rad(theta_deg_at_t) 
    
    rotation_matrice = np.array([ [np.cos(angle_rad), -np.sin(angle_rad)], [np.sin(angle_rad), np.cos(angle_rad)] ])

    #transformation = np.dot(np.array(source.exterior.coords),rotation_matrice.T) + translation * t
    transformation = np.dot(src_centered + deformation * t,rotation_matrice.T) + (translation * t + src_centroid) #- deformation * t
    transformed_polygon = Polygon(transformation)    
    
    # alternative
    #transformed_polygon = affine_transform(source, [np.cos(angle_rad), -np.sin(angle_rad), np.sin(angle_rad), np.cos(angle_rad), translation[0] * t, translation[1] * t])

    return transformed_polygon

