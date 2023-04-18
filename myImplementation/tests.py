def getFeaturePoints(points, dist_range):
    feature_points = []

    for point in points:
        if point[0]+point[1] > dist_range:
            feature_points.append(point)

    return feature_points

arrPts = [(0,1),(1,1),(3,3),(2,1)]

# Open file for writing
f= open('testsss.txt', 'w')
for i in range(10, 50):
    for j in range(150,300):
        arrE = getFeaturePoints(arrPts, i)
        totalPoints = len(arrE)
        f.write(f"dist_value: {i} and total_points: {totalPoints}\n")
