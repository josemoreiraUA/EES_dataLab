import numpy as np
from shapely.geometry import Polygon, LineString, Point
import shapely.geometry as sg
from shapely.ops import unary_union
import shapely.wkt
import geopandas as gpd
import matplotlib.pyplot as plt
import math
from numpy.linalg import eig

'''#arrSrcPoints = [(1,2),(2,4),(3,6),(4,8),(5,10),(6,12),(7,14),(8,16),(9,18),(10,20)]
arrSrcPoints = [(1,2),(4,8),(5,10),(6,12),(7,14),(8,16),(9,18),(10,20)]
arrTrgPoints = [(1,2),(4,8),(5,20),(16,32),(19,21),(36,72),(49,98),(64,128),(100,200)]
arrFPSrc = [(1,2),(4,8),(7,14),(8,16)]
arrFPTrg = [(1,2),(4,8),(16,32),(49,98)]
arrCorrFP = [[(1,2),(4,8)],[(4,8),(16,32)],[(7,14),(49,98)],[(8,16),(64,128)]]
#arrCorrP = arrCorrFP.copy()
arrCorrP = []
#if len(arrSrcPoints) > len(arrTrgPoints):
for corr in range(len(arrCorrFP)-1):

    if corr == 0: #na primeira correspondencia
        idxSource = []
        idxTarget = []
        arrSrc = []
        arrTrg = []
        idxS = arrSrcPoints.index(arrCorrFP[corr][0])
        idxT = arrTrgPoints.index(arrCorrFP[corr][1])
        if idxS != 0: #se o idx do primeiro feature point source, for diferente de 0 no source points, quer dizer que ha pontos antes
            for k in range(idxS): #para cada idx ate ao primeiro FP source
                idxSource.append(k) #guarda neste array
        if idxT != 0: #se o idx do primeiro feature point target for diferente de 0 no target points, quer dizer que ha pontos antes
            for k in range(idxT):
                idxTarget.append(k)
        arrSrc = [arrSrcPoints[j] for j in idxSource] #arrSrc sao os numero do arrSrcPoints
        arrTrg = [arrTrgPoints[j] for j in idxTarget]

        if len(idxSource) > 0 or len(idxTarget) > 0: #se houver pontos no source ou no target

            if len(idxSource) > len(idxTarget):
                fct = len(idxTarget) / len(idxSource)
                for i, s in enumerate(arrSrc):
                    if fct==0: #se o fator for 0, isto e', se o arrT=0 (nao houver pontos entre os 2 fp do source)
                        arrCorrP.append([s,arrCorrFP[corr][1]])
                    else:
                        corrIdx = int(i * fct)
                        arrCorrP.append([s, arrTrg[corrIdx]])
            else:
                fct = len(idxSource) / len(idxTarget)
                for i, s in enumerate(arrTrg):
                    if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                        arrCorrP.append([arrCorrFP[corr][0],s])
                    else:
                        corrIdx = int(i * fct)
                        #print([arrS[corrIdx],s])
                        arrCorrP.append([arrSrc[corrIdx],s])

    arrCorrP.append(arrCorrFP[corr])
    idxS = arrSrcPoints.index(arrCorrFP[corr][0])
    idxSnext = arrSrcPoints.index(arrCorrFP[corr+1][0])
    idxT = arrTrgPoints.index((arrCorrFP[corr][1]))
    idxTnext = arrTrgPoints.index((arrCorrFP[corr+1][1]))

    numbersSrc = []
    for i in range(idxS+1,idxSnext):
        numbersSrc.append(i)

    numbersTrg = []
    for i in range(idxT+1,idxTnext):
        numbersTrg.append(i)
    
    arrS = [arrSrcPoints[j] for j in numbersSrc]
    arrT = [arrTrgPoints[j] for j in numbersTrg]

    if len(arrS) < len(arrT):
        print("arrS<arrT")
        fator = len(arrS) / len(arrT)
        print("fator",fator)
        for i, s in enumerate(arrT):
            if fator==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                arrCorrP.append([arrCorrFP[corr][0],s])
            else:
                print(i)
                print(s)
                corrIdx = int(i * fator)
                #print([arrS[corrIdx],s])
                arrCorrP.append([arrS[corrIdx],s])
                    
    elif len(arrS) > len(arrT):
        print("arrS>arrT")
        fator = len(arrT) / len(arrS)
        for i, s in enumerate(arrS):
            if fator==0: #se o fator for 0, isto e', se o arrT=0 (nao houver pontos entre os 2 fp do source)
                arrCorrP.append([s,arrCorrFP[corr][1]])
            else:
                print(i)
                print(s)
                corrIdx = int(i * fator)
                #print(([s, arrT[corrIdx]]))
                arrCorrP.append([s, arrT[corrIdx]])

    else:
        print("arrS=arrT")
        fator = 1
        for i, s in enumerate(arrS):
            print(i)
            print(s)
            corrIdx = int(i * fator)
            #print([s, arrT[corrIdx]])
            arrCorrP.append([s, arrT[corrIdx]])

    if(corr == (len(arrCorrFP) - 2)):
        arrCorrP.append(arrCorrFP[corr+1])

        idxSource = []
        idxTarget = []
        arrSrc = []
        arrTrg = []
        idxS = arrSrcPoints.index(arrCorrFP[corr+1][0])
        idxT = arrTrgPoints.index(arrCorrFP[corr+1][1])
        
        if idxS != len(arrSrcPoints)-1: #se o idx do ultimo feature point source, for diferente do idx do ultimo point source, quer dizer que ha pontos depois
            for l in range(idxS+1, len(arrSrcPoints)): #para cada idx ate ao ultimo FP source
                idxSource.append(l) #guarda neste array
        if idxT != len(arrTrgPoints)-1: #se o idx do ultimo feature point target, for diferente do idx do ultimo point target, quer dizer que ha pontos depois
            for l in range(idxT+1, len(arrTrgPoints)):
                idxTarget.append(l)

        arrSrc = [arrSrcPoints[j] for j in idxSource] #arrSrc sao os numero do arrSrcPoints
        arrTrg = [arrTrgPoints[j] for j in idxTarget]

        if len(idxSource) > 0 or len(idxTarget) > 0: #se houver pontos no source ou no target

            if len(idxSource) > len(idxTarget):
                fct = len(idxTarget) / len(idxSource)
                for i, s in enumerate(arrSrc):
                    if fct==0: #se o fator for 0, isto e', se o arrT=0 (nao houver pontos entre os 2 fp do source)
                        arrCorrP.append([s,arrCorrFP[corr][1]])
                    else:
                        corrIdx = int(i * fct)
                        arrCorrP.append([s, arrTrg[corrIdx]])
            else:
                fct = len(idxSource) / len(idxTarget)
                for i, s in enumerate(arrTrg):
                    if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                        arrCorrP.append([arrCorrFP[corr][0],s])
                    else:
                        corrIdx = int(i * fct)
                        arrCorrP.append([arrSrc[corrIdx],s])

    print("::::::::")
print(arrCorrP)'''

def divide_line_segment(start_point, end_point, num_parts, arrPts, arrC):
    # Calculate the length of the line segment
    segment_length = ((end_point[0] - start_point[0]) ** 2 + (end_point[1] - start_point[1]) ** 2) ** 0.5
    
    # Calculate the distance between each point
    interval_length = segment_length / num_parts
    
    # Calculate the direction vector of the line segment
    dx = end_point[0] - start_point[0]
    dy = end_point[1] - start_point[1]
    
    # Normalize the direction vector
    direction_vector = (dx / segment_length, dy / segment_length)
    
    # Calculate the starting points for each of the equal parts of the line segment
    result = []
    for i in range(num_parts):
        # Calculate the distance from the start point to the current part
        distance = interval_length * i
        
        # Calculate the coordinates of the starting point of the current part
        x = start_point[0] + distance * direction_vector[0]
        y = start_point[1] + distance * direction_vector[1]
        
        result.append((x, y))
    print(result)
    corr = [] #correspondencia entre source e target    
    for i in range(len(arrPts)):
        corr=[arrPts[i],(result[i])]
        arrC.append(corr)

    return arrC

start_point = (1, 2)
end_point = (2, 3)
num_parts = 5
arrPointsTrg = [(1,1),(2,1),(3,1),(4,1),(5,1)]
arrCorr = []
starting_points = divide_line_segment(start_point, end_point, num_parts, arrPointsTrg, arrCorr)
print(starting_points)

