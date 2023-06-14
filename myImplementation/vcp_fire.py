from shapely.geometry import Polygon
import shapely.geometry as sg
from shapely.ops import unary_union
import shapely.wkt
import geopandas as gpd
import matplotlib.pyplot as plt
import math
from ros import Ros
from featurepoint import FeaturePoint
from polygonCorrespondence import PolygonCorrespondence
import re
import time

#variaveis globais
arrPolyWKT = [] #array com os poligonos wkt
arrPolyGPD = [] #array com os poligonos gpd (para efeitos de visualizacao)

def readWKT(wktFile): #le os ficheiros wkt e guarda os poligonos num array

    f = open(wktFile) 
    wktCoords = f.read() #wktCoords fica com as coordenadas do poligono

    wktPol = shapely.wkt.loads(wktCoords) #cria um poligono wkt
    arrPolyWKT.append(wktPol) #acrescenta o poligono do ficheiro ao array de poligonos
    gpdPol = gpd.GeoSeries([wktPol])
    arrPolyGPD.append(gpdPol) #igual ao wkt mas para fazer graficos
    
def getFeaturePoints(polygon, dmin, amax): #obtem os featurePoints do poligono

    arrFeaturePoints = []
    arestas = []
    dmax = 0 #distancia maxima e' considerada a maior aresta do poligono

    for j in range(len(polygon.exterior.coords)-1): #para cada ponto no poligono, em ordem a calcular tamanho medio das arestas (para o dmin e dmax)
        arestas.append(math.dist(polygon.exterior.coords[j], polygon.exterior.coords[j+1])) #distancia entre os 2 pontos
        if math.dist(polygon.exterior.coords[j], polygon.exterior.coords[j+1]) > dmax:
            dmax = math.dist(polygon.exterior.coords[j], polygon.exterior.coords[j+1]) #dmax e a maior aresta do poligono

    for k in range(1, len(polygon.exterior.coords)-1): #para cada ponto do poligono, em ordem a descobrir se e feature point
        if k == len(polygon.exterior.coords)-1: #se o k for a ultima posicao do array, usamos o primeiro ponto como Pi+
            distAB = math.dist(polygon.exterior.coords[k-1], polygon.exterior.coords[k]) #distancia Pi- e Pi
            distBC = math.dist(polygon.exterior.coords[k], polygon.exterior.coords[0]) #distancia Pi e Pi+
            distAC = math.dist(polygon.exterior.coords[k-1], polygon.exterior.coords[0]) #distancia Pi- e Pi+
        elif k == 0: #se o k for a primeira posicao do array, usamos o ultimo ponto como Pi-
            distAB = math.dist(polygon.exterior.coords[len(polygon.exterior.coords)-1], polygon.exterior.coords[k]) #distancia Pi- e Pi
            distBC = math.dist(polygon.exterior.coords[k], polygon.exterior.coords[k+1]) #distancia Pi e Pi+
            distAC = math.dist(polygon.exterior.coords[len(polygon.exterior.coords)-1], polygon.exterior.coords[0]) #distancia Pi- e Pi+
        else: #senao faz normalmente (Pi+ e' o seguinte e Pi- e' o anterior)
            distAB = math.dist(polygon.exterior.coords[k-1], polygon.exterior.coords[k]) #distancia Pi- e Pi
            distBC = math.dist(polygon.exterior.coords[k], polygon.exterior.coords[k+1]) #distancia Pi e Pi+
            distAC = math.dist(polygon.exterior.coords[k-1], polygon.exterior.coords[k+1]) #distancia Pi- e Pi+
        valueABC = (pow(distBC,2) + pow(distAB,2) - pow(distAC,2)) / (2*distBC*distAB)
        
        if valueABC < -1.0: #para os casos tipo 1.0000002
            valueABC = -1.0

        angleABCc = math.acos(valueABC) #arccos para calculo do angulo alpha
        angleABC = math.degrees(angleABCc)
        if (distBC >= dmin and distBC <= dmax) and (distAB >= dmin and distAB <= dmax) and (angleABC <= amax): #se o ponto preencher todos os requisitos
            arrFeaturePoints.append(polygon.exterior.coords[k]) #e adicionado ao array de feature point
        arrFeaturePointsWithoutDuplicates = []
        [arrFeaturePointsWithoutDuplicates.append(x) for x in arrFeaturePoints if x not in arrFeaturePointsWithoutDuplicates] #garante q nao ha pontos repetidos

    return arrFeaturePointsWithoutDuplicates

def getROS(arrOfFeaturePoints, idxPoint): #obtem as regions of support do poligono
    ros = [] #region of support
    rol = [] #region of left
    ror = [] #region of right

    if idxPoint==0: #se estiver na primeira posicao, o feature point anterior e o ultimo do array
        ros.append(arrOfFeaturePoints[len(arrOfFeaturePoints)-1])
        ros.append(arrOfFeaturePoints[idxPoint])
        ros.append(arrOfFeaturePoints[idxPoint+1])
        rol.append(arrOfFeaturePoints[len(arrOfFeaturePoints)-1])
        rol.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[idxPoint+1])
        
    elif idxPoint == len(arrOfFeaturePoints)-1: #se estiver na ultima posicao, o feature point seguinte e o primeiro do array
        ros.append(arrOfFeaturePoints[idxPoint-1])
        ros.append(arrOfFeaturePoints[idxPoint])
        ros.append(arrOfFeaturePoints[0])
        rol.append(arrOfFeaturePoints[idxPoint-1])
        rol.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[0])
        
    else:
        ros.append(arrOfFeaturePoints[idxPoint-1])
        ros.append(arrOfFeaturePoints[idxPoint])
        ros.append(arrOfFeaturePoints[idxPoint+1])
        rol.append(arrOfFeaturePoints[idxPoint-1])
        rol.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[idxPoint])
        ror.append(arrOfFeaturePoints[idxPoint+1])

    return ros, rol, ror

def divide_line_segment_src(start_point, end_point, num_parts, arrPts, arrC): #funcao para dividir o segmento em partes iguais quando nao ha pontos normais entre fp (para evitar repeticao nas correspondencias)
    # Tamanho do segmento de reta
    segment_length = ((end_point[0] - start_point[0]) ** 2 + (end_point[1] - start_point[1]) ** 2) ** 0.5
    
    # Distancia entre cada ponto
    interval_length = segment_length / num_parts
    
    # vetor direcao do segmento de reta
    dx = end_point[0] - start_point[0]
    dy = end_point[1] - start_point[1]
    
    # normalizacao do vetor direcao
    direction_vector = (dx / segment_length, dy / segment_length)
    
    # pontos iniciais para cada segmento
    result = []
    for i in range(num_parts):
        # distancia do ponto inicial ate ao atual
        distance = interval_length * i
        
        # coordenadas de inicio do segmento atual
        x = start_point[0] + distance * direction_vector[0]
        y = start_point[1] + distance * direction_vector[1]
        
        result.append((x, y))
        
    corr = [] #correspondencia entre source e target    
    for i in range(len(arrPts)):
        corr=[(result[i]),arrPts[i]]
        arrC.append(corr)
        
    return arrC

def divide_line_segment_trg(start_point, end_point, num_parts, arrPts, arrC): #funcao para dividir o segmento em partes iguais quando nao ha pontos normais entre fp (para evitar repeticao nas correspondencias)
    # igual ao anterior mas para o target
    segment_length = ((end_point[0] - start_point[0]) ** 2 + (end_point[1] - start_point[1]) ** 2) ** 0.5
    
    interval_length = segment_length / num_parts
    
    dx = end_point[0] - start_point[0]
    dy = end_point[1] - start_point[1]
    
    direction_vector = (dx / segment_length, dy / segment_length)
    
    result = []
    for i in range(num_parts):
        distance = interval_length * i
        
        x = start_point[0] + distance * direction_vector[0]
        y = start_point[1] + distance * direction_vector[1]
        
        result.append((x, y))
        
    corr = [] #correspondencia entre source e target    
    for i in range(len(arrPts)):
        corr=[arrPts[i],(result[i])]
        arrC.append(corr)
        
    return arrC

def getIntermediateCorrespondences(polSrc, polTrg, arrCorrFP): #pontos pol src, pontos pol trg, arr de correspondencias
    s = ""
    for corr in arrCorrFP:
        s = s + str(corr) #para cada correspondencia, concatena na mesma string

    pattern = r"\(\d+\.\d+,\s*\d+\.\d+\)"
    matches = re.findall(pattern, s) #fica so com os pontos
    
    arrCorrFP = [] #array onde vao ficar as correspondencias entre feature points
    
    for i in range(0, len(matches), 2): #para cada ponto
        x = tuple(map(float, matches[i][1:-1].split(",")))
        y = tuple(map(float, matches[i+1][1:-1].split(",")))
        arrCorrFP.append([x,y]) #array no formato [[(x1s,y1s),(x1t,y1t)],[(x2s,y2s),(x2t,y2t)]] -> arrCorrFP
    
    arrCorrFPtmp = arrCorrFP.copy() #usamos um tmp so pra eliminar fp com mais que uma correspondencia

    for i in range(1,len(arrCorrFP)): #para cada correspondencia entre FP
        if arrCorrFP[i][0] == arrCorrFP[i-1][0]: #se o ponto source for igual ao anterior
            arrCorrFPtmp.remove(arrCorrFP[i]) #remove essa correspondencia, e os pontos passam a ser tratados como pontos normais (e nao FP)
        if arrCorrFP[i][1] == arrCorrFP[i-1][1]: #se o ponto target for igual ao anterior
            arrCorrFPtmp.remove(arrCorrFP[i]) #remove essa correspondencia, e os pontos passam a ser tratados como pontos normais (e nao FP)

    arrCorrFP = arrCorrFPtmp.copy() #o array de correspondencias passa a ser o tmp
    arrCorrP = [] #array onde vao ficar as correspondencias entre pontos
    
    arrSrcPoints = polSrc.exterior.coords[:]
    arrTrgPoints = polTrg.exterior.coords[:]
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
                            arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([s, arrTrg[corrIdx]])
                else:
                    fct = len(idxSource) / len(idxTarget)
                    for i, s in enumerate(arrTrg):
                        if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                            arrCorrP = divide_line_segment_trg(arrCorrFP[corr][1],arrCorrFP[corr+1][1],len(arrS),arrS,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([arrSrc[corrIdx],s])

        arrCorrP.append(arrCorrFP[corr]) #adiciona a correspondencia (entre fp) ao novo array de correspondencias entre pontos
        
        idxS = arrSrcPoints.index(arrCorrFP[corr][0]) #obtem o idx do feature point ponto no arr Source
        idxSnext = arrSrcPoints.index(arrCorrFP[corr+1][0]) #obtem o idx do proximo feature point no arr Source
        idxT = arrTrgPoints.index((arrCorrFP[corr][1])) #igual mas para o target
        idxTnext = arrTrgPoints.index((arrCorrFP[corr+1][1]))
        
        numbersSrc = []
        for i in range(idxS+1,idxSnext):
            numbersSrc.append(i)  #idx dos pontos nao feature point no arrSrc

        numbersTrg = []
        for i in range(idxT+1,idxTnext):
            numbersTrg.append(i)

        arrS = [arrSrcPoints[j] for j in numbersSrc] #arr de pontos que nao sao fp 
        arrT = [arrTrgPoints[j] for j in numbersTrg]

        if (len(arrS)==0) and (len(arrT)) == 0: #se nao houver pontos entre os 2 fp src e trg, passa à proxima correspondencia
            continue
        
        if len(arrS) < len(arrT):
            fator = (len(arrS)) / (len(arrT))
            for i, s in enumerate(arrT): #para cada ponto no array de pontos nao fp
                if fator==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                    arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                    break
                else:
                    corrIdx = int(i * fator) #a correspondencia e' feita com
                    arrCorrP.append([arrS[corrIdx],s]) #e' adicionado ao array de correspondencias entre pontos a correspond 
                        
        else:
            fator = len(arrT) / len(arrS)
            for i, s in enumerate(arrS):
                if fator==0: #se o fator for 0, isto e', se o arrT=0 (nao houver pontos entre os 2 fp do target)
                    arrCorrP = divide_line_segment_trg(arrCorrFP[corr][1],arrCorrFP[corr+1][1],len(arrS),arrS,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                    break
                else:
                    corrIdx = int(i * fator)
                    arrCorrP.append([s, arrT[corrIdx]])
        
        if(corr == (len(arrCorrFP) - 2)): #se ja estiver na penultima correspondencia
            arrCorrP.append(arrCorrFP[corr+1]) #adiciona a ultima correspondencia entre fp

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
                            arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([s, arrTrg[corrIdx]])
                else:
                    fct = len(idxSource) / len(idxTarget)
                    for i, s in enumerate(arrTrg):
                        if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                            arrCorrP = divide_line_segment_trg(arrCorrFP[corr][1],arrCorrFP[corr+1][1],len(arrS),arrS,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([arrSrc[corrIdx],s])

    return arrCorrP    

def getJaccardIndex(source, target, correspondences):
    
    middlePoints = [] #pontos do poligono medio
    jaccardIndex=0
    if(len(correspondences)==0): #se nao houver correspondencias acaba (para quando nao faz o catch do zerodivisionerror)
       jaccardIndex = 0
       return jaccardIndex
    
    for j in range(len(correspondences)): #para cada correspondencia
        ptSrc = correspondences[j][0]
        ptTrg = correspondences[j][1]
        coordX = ptSrc[0] + ((ptTrg[0] - ptSrc[0]) / 2) #calcula o ponto medio
        coordY = ptSrc[1] + ((ptTrg[1] - ptSrc[1]) / 2)
        middlePoints.append((coordX, coordY)) #e coloca no array 
    midPoly = Polygon([p[0], p[1]] for p in middlePoints)

    #para evitar erros de topologia nos poligonos mal formados
    #se calhar devia fazer-se um try except a passar estes a' frente
    if(not source.is_valid):
        source = source.buffer(0)
    elif(not target.is_valid):
        target = target.buffer(0)
    elif(not midPoly.is_valid):
        midPoly = midPoly.buffer(0)

    #calculo do jaccard index do source para o middle
    intersectionSM = source.intersection(midPoly).area
    unionSM = unary_union([source, midPoly]).area
    jaccardIndexSM = intersectionSM/unionSM

    #calculo do jaccard index do middle para o target
    intersectionMT = midPoly.intersection(target).area
    unionMT = unary_union([midPoly, target]).area
    jaccardIndexMT = intersectionMT/unionMT
    
    jaccardIndex = (jaccardIndexSM + jaccardIndexMT) / 2
    print(jaccardIndex)

    return jaccardIndex

arrPontos = []
#wktFiles = ["dataset/fireSPTDatalab/f_spt_dl0.wkt","dataset/fireSPTDatalab/f_spt_dl1.wkt","dataset/fireSPTDatalab/f_spt_dl2.wkt","dataset/fireSPTDatalab/f_spt_dl3.wkt","dataset/fireSPTDatalab/f_spt_dl4.wkt","dataset/fireSPTDatalab/f_spt_dl5.wkt","dataset/fireSPTDatalab/f_spt_dl6.wkt","dataset/fireSPTDatalab/f_spt_dl7.wkt","dataset/fireSPTDatalab/f_spt_dl8.wkt"]
#wktFiles = ["simp0 (1).wkt","simp0 (2).wkt","simp0 (3).wkt","simp0 (4).wkt","simp0 (5).wkt","simp0 (6).wkt","simp0 (7).wkt","simp0 (8).wkt","simp0 (9).wkt"]
'''wktFiles = [
    "tigas226/polygon_1.wkt", "tigas226/polygon_2.wkt", "tigas226/polygon_3.wkt", "tigas226/polygon_4.wkt", "tigas226/polygon_5.wkt", 
    "tigas226/polygon_6.wkt", "tigas226/polygon_7.wkt", "tigas226/polygon_8.wkt", "tigas226/polygon_9.wkt", "tigas226/polygon_10.wkt", 
    "tigas226/polygon_11.wkt", "tigas226/polygon_12.wkt", "tigas226/polygon_13.wkt", "tigas226/polygon_14.wkt", "tigas226/polygon_15.wkt", 
    "tigas226/polygon_16.wkt", "tigas226/polygon_17.wkt", "tigas226/polygon_18.wkt", "tigas226/polygon_19.wkt", "tigas226/polygon_20.wkt", 
    "tigas226/polygon_21.wkt", "tigas226/polygon_22.wkt", "tigas226/polygon_23.wkt", "tigas226/polygon_24.wkt", "tigas226/polygon_25.wkt", 
    "tigas226/polygon_26.wkt", "tigas226/polygon_27.wkt", "tigas226/polygon_28.wkt", "tigas226/polygon_29.wkt", "tigas226/polygon_30.wkt", 
    "tigas226/polygon_31.wkt", "tigas226/polygon_32.wkt", "tigas226/polygon_33.wkt", "tigas226/polygon_34.wkt", "tigas226/polygon_35.wkt", 
    "tigas226/polygon_36.wkt", "tigas226/polygon_37.wkt", "tigas226/polygon_38.wkt", "tigas226/polygon_39.wkt", "tigas226/polygon_40.wkt", 
    "tigas226/polygon_41.wkt", "tigas226/polygon_42.wkt", "tigas226/polygon_43.wkt", "tigas226/polygon_44.wkt", "tigas226/polygon_45.wkt", 
    "tigas226/polygon_46.wkt", "tigas226/polygon_47.wkt", "tigas226/polygon_48.wkt", "tigas226/polygon_49.wkt", "tigas226/polygon_50.wkt", 
    "tigas226/polygon_51.wkt", "tigas226/polygon_52.wkt", "tigas226/polygon_53.wkt", "tigas226/polygon_54.wkt", "tigas226/polygon_55.wkt", 
    "tigas226/polygon_56.wkt", "tigas226/polygon_57.wkt", "tigas226/polygon_58.wkt", "tigas226/polygon_59.wkt", "tigas226/polygon_60.wkt", 
    "tigas226/polygon_61.wkt", "tigas226/polygon_62.wkt", "tigas226/polygon_63.wkt", "tigas226/polygon_64.wkt", "tigas226/polygon_65.wkt", 
    "tigas226/polygon_66.wkt", "tigas226/polygon_67.wkt", "tigas226/polygon_68.wkt", "tigas226/polygon_69.wkt", "tigas226/polygon_70.wkt", 
    "tigas226/polygon_71.wkt", "tigas226/polygon_72.wkt", "tigas226/polygon_73.wkt", "tigas226/polygon_74.wkt", "tigas226/polygon_75.wkt", 
    "tigas226/polygon_76.wkt", "tigas226/polygon_77.wkt", "tigas226/polygon_78.wkt", "tigas226/polygon_79.wkt", "tigas226/polygon_80.wkt", 
    "tigas226/polygon_81.wkt", "tigas226/polygon_82.wkt", "tigas226/polygon_83.wkt", "tigas226/polygon_84.wkt", "tigas226/polygon_85.wkt", 
    "tigas226/polygon_86.wkt", "tigas226/polygon_87.wkt", "tigas226/polygon_88.wkt", "tigas226/polygon_89.wkt", "tigas226/polygon_90.wkt", 
    "tigas226/polygon_91.wkt", "tigas226/polygon_92.wkt", "tigas226/polygon_93.wkt", "tigas226/polygon_94.wkt", "tigas226/polygon_95.wkt", 
    "tigas226/polygon_96.wkt", "tigas226/polygon_97.wkt", "tigas226/polygon_98.wkt", "tigas226/polygon_99.wkt", "tigas226/polygon_100.wkt", 
    "tigas226/polygon_101.wkt", "tigas226/polygon_102.wkt", "tigas226/polygon_103.wkt", "tigas226/polygon_104.wkt", "tigas226/polygon_105.wkt", 
    "tigas226/polygon_106.wkt", "tigas226/polygon_107.wkt", "tigas226/polygon_108.wkt", "tigas226/polygon_109.wkt", "tigas226/polygon_110.wkt", 
    "tigas226/polygon_111.wkt", "tigas226/polygon_112.wkt", "tigas226/polygon_113.wkt", "tigas226/polygon_114.wkt", "tigas226/polygon_115.wkt", 
    "tigas226/polygon_116.wkt", "tigas226/polygon_117.wkt", "tigas226/polygon_118.wkt", "tigas226/polygon_119.wkt", "tigas226/polygon_120.wkt", 
    "tigas226/polygon_121.wkt", "tigas226/polygon_122.wkt", "tigas226/polygon_123.wkt", "tigas226/polygon_124.wkt", "tigas226/polygon_125.wkt", 
    "tigas226/polygon_126.wkt", "tigas226/polygon_127.wkt", "tigas226/polygon_128.wkt", "tigas226/polygon_129.wkt", "tigas226/polygon_130.wkt", 
    "tigas226/polygon_131.wkt", "tigas226/polygon_132.wkt", "tigas226/polygon_133.wkt", "tigas226/polygon_134.wkt", "tigas226/polygon_135.wkt", 
    "tigas226/polygon_136.wkt", "tigas226/polygon_137.wkt", "tigas226/polygon_138.wkt", "tigas226/polygon_139.wkt", "tigas226/polygon_140.wkt", 
    "tigas226/polygon_141.wkt", "tigas226/polygon_142.wkt", "tigas226/polygon_143.wkt", "tigas226/polygon_144.wkt", "tigas226/polygon_145.wkt", 
    "tigas226/polygon_146.wkt", "tigas226/polygon_147.wkt", "tigas226/polygon_148.wkt", "tigas226/polygon_149.wkt", "tigas226/polygon_150.wkt", 
    "tigas226/polygon_151.wkt", "tigas226/polygon_152.wkt", "tigas226/polygon_153.wkt", "tigas226/polygon_154.wkt", "tigas226/polygon_155.wkt", 
    "tigas226/polygon_156.wkt", "tigas226/polygon_157.wkt", "tigas226/polygon_158.wkt", "tigas226/polygon_159.wkt", "tigas226/polygon_160.wkt", 
    "tigas226/polygon_161.wkt", "tigas226/polygon_162.wkt", "tigas226/polygon_163.wkt", "tigas226/polygon_164.wkt", "tigas226/polygon_165.wkt", 
    "tigas226/polygon_166.wkt", "tigas226/polygon_167.wkt", "tigas226/polygon_168.wkt", "tigas226/polygon_169.wkt", "tigas226/polygon_170.wkt", 
    "tigas226/polygon_171.wkt", "tigas226/polygon_172.wkt", "tigas226/polygon_173.wkt", "tigas226/polygon_174.wkt", "tigas226/polygon_175.wkt", 
    "tigas226/polygon_176.wkt", "tigas226/polygon_177.wkt", "tigas226/polygon_178.wkt", "tigas226/polygon_179.wkt", "tigas226/polygon_180.wkt", 
    "tigas226/polygon_181.wkt", "tigas226/polygon_182.wkt", "tigas226/polygon_183.wkt", "tigas226/polygon_184.wkt", "tigas226/polygon_185.wkt", 
    "tigas226/polygon_186.wkt", "tigas226/polygon_187.wkt", "tigas226/polygon_188.wkt", "tigas226/polygon_189.wkt", "tigas226/polygon_190.wkt", 
    "tigas226/polygon_191.wkt", "tigas226/polygon_192.wkt", "tigas226/polygon_193.wkt", "tigas226/polygon_194.wkt", "tigas226/polygon_195.wkt", 
    "tigas226/polygon_196.wkt", "tigas226/polygon_197.wkt", "tigas226/polygon_198.wkt", "tigas226/polygon_199.wkt", "tigas226/polygon_200.wkt", 
    "tigas226/polygon_201.wkt", "tigas226/polygon_202.wkt", "tigas226/polygon_203.wkt", "tigas226/polygon_204.wkt", "tigas226/polygon_205.wkt", 
    "tigas226/polygon_206.wkt", "tigas226/polygon_207.wkt", "tigas226/polygon_208.wkt", "tigas226/polygon_209.wkt", "tigas226/polygon_210.wkt", 
    "tigas226/polygon_211.wkt", "tigas226/polygon_212.wkt", "tigas226/polygon_213.wkt", "tigas226/polygon_214.wkt", "tigas226/polygon_215.wkt", 
    "tigas226/polygon_216.wkt", "tigas226/polygon_217.wkt", "tigas226/polygon_218.wkt", "tigas226/polygon_219.wkt", "tigas226/polygon_220.wkt", 
    "tigas226/polygon_221.wkt", "tigas226/polygon_222.wkt", "tigas226/polygon_223.wkt", "tigas226/polygon_224.wkt", "tigas226/polygon_225.wkt", 
    "tigas226/polygon_226.wkt"
]'''
wktFiles = ["dataset/wkt/test1.wkt", "dataset/wkt/test2.wkt", "dataset/wkt/test3.wkt"]
valuesTiago=[] #array onde sao armazenados os valores da solucao (distancia e angulo)
correspondences = []
correspondencesBetweenAllPoints = [] #correspondencias entre todos os pontos, nao so feature points
jaccard = 0 #jaccard index medio entre S->M e M->T

for file in wktFiles:
    readWKT(file) #popula o array de poligonos (arrPolyWKT)

f=open('testResults.txt', 'w') #nome do ficheiro onde vai ser guardada a solucao
for pol in range(len(arrPolyWKT)-1): #para cada par de poligonos
    minJaccard = 0 
    for angulo in range(120,181,10): #testa para valores de ang  
        for distancia in range(1,20): #e de dist

            f.write(f"pol: {pol}, ang: {angulo}, dist: {distancia} \n")

            try:
                arrFPObjSource = [] #a cada poligono os feature points sao resetados
                arrFPObjTarget = []
                correspondences = []

                if((angulo==170 or angulo==180) and distancia<3): #para dist inferior a 3 so faz angulos menores que 170
                    f.write(f"Too many values\n")
                    continue

                start_time = time.time()
                featurePointsSource = getFeaturePoints(arrPolyWKT[pol], distancia, angulo)
                featurePointsTarget = getFeaturePoints(arrPolyWKT[pol+1], distancia, angulo)

                if time.time() - start_time > 120: #se demorar mais de 2 minutos a calcular correspondencias pode passar ao proximo par de valores (deve passar ao prox ang)
                    f.write(f"to many feature points after feature points. skipping to next cycle \n")
                    continue

                maxi = 0
                mini = 123445612456

                for i in range(len(featurePointsSource)):
                    ros, rol, ror = getROS(featurePointsSource, i) #calcula as regions do ponto
                    regOfSupport = Ros(ros, featurePointsSource[i],arrPolyWKT[pol]) #cria um objeto ROS
                    featPoint = FeaturePoint(featurePointsSource[i],arrPolyWKT[pol], ros, rol, ror) #cria um objeto FP

                    if i == len(featurePointsSource)-1:
                        valueR = math.dist(featurePointsSource[i], featurePointsSource[0])
                    else:
                        valueR = math.dist(featurePointsSource[i], featurePointsSource[i+1])
                        
                    if valueR > maxi:
                        maxi = valueR
                    if valueR < mini:
                        mini = valueR

                    featPoint.rorSize = valueR

                    if i == 0:
                        valueL = math.dist(featurePointsSource[i], featurePointsSource[len(featurePointsSource)-1])
                    else:
                        valueL = math.dist(featurePointsSource[i], featurePointsSource[i-1])

                    if valueL > maxi:
                        maxi = valueL
                    if valueL < mini:
                        mini = valueL

                    featPoint.rolSize = valueL
                    arrFPObjSource.append(featPoint) #adiciona o FP ao array de objetos FP

                maxiTarget = 0
                miniTarget = 123445612456

                for i in range(len(featurePointsTarget)):
                    ros, rol, ror = getROS(featurePointsTarget, i)
                    regOfSupport = Ros(ros, featurePointsTarget[i],arrPolyWKT[pol+1])
                    featPoint = FeaturePoint(featurePointsTarget[i],arrPolyWKT[pol+1], ros, rol, ror)

                    if i == len(featurePointsTarget)-1:
                        valueR = math.dist(featurePointsTarget[i], featurePointsTarget[0])
                    else:
                        valueR = math.dist(featurePointsTarget[i], featurePointsTarget[i+1])
                        
                    if valueR > maxiTarget:
                        maxiTarget = valueR
                    if valueR < miniTarget:
                        miniTarget = valueR

                    featPoint.rorSize = valueR

                    if i == 0:
                        valueL = math.dist(featurePointsTarget[i], featurePointsTarget[len(featurePointsTarget)-1])
                    else:
                        valueL = math.dist(featurePointsTarget[i], featurePointsTarget[i-1])

                    if valueL > maxiTarget:
                        maxiTarget = valueL
                    if valueL < miniTarget:
                        miniTarget = valueL

                    featPoint.rolSize = valueL
                    arrFPObjTarget.append(featPoint) #adiciona o FP ao array de objetos FP
                
                for fp in arrFPObjSource: #calcula as features de cada FP do source
                    fp.rolSizeNorm = (fp.rolSize-mini) / (maxi - mini)
                    fp.rorSizeNorm = (fp.rorSize-mini) / (maxi - mini)
                    fp.getFeatures(fp.ros)

                for fp in arrFPObjTarget: #calcula as features de cada FP do target
                    fp.rolSizeNorm = (fp.rolSize-miniTarget) / (maxiTarget - miniTarget)
                    fp.rorSizeNorm = (fp.rorSize-miniTarget) / (maxiTarget - miniTarget)
                    fp.getFeatures(fp.ros)
                
                if time.time() - start_time > 120: #se demorar mais de 2 minutos a calcular correspondencias pode passar ao proximo par de valores (deve passar ao prox ang)
                    f.write(f"to many feature points after correspondences. skipping to next cycle \n")
                    continue

                #calcula as correspondencias entre cada FP do poligono atual e do seguinte
                # compute correspondences
                pc1 = PolygonCorrespondence(arrFPObjSource, arrFPObjTarget, .3, .3, .4, 1)
                correspondencesBetweenAllPoints = getIntermediateCorrespondences(arrPolyWKT[pol], arrPolyWKT[pol+1], pc1.getFPCorrespondences(3))
                correspondences.append(correspondencesBetweenAllPoints)

            except ZeroDivisionError:
                f.write(f"Not enough feature points. Skipping to next cycle.\n")
                continue  # Skip to next cycle
            except IndexError:
                f.write(f"Not enough feature points. Skipping to next cycle.\n")
                continue  # Skip to next cycle

            totalDistEntrePols = 0 #dist entre correspondencias
            numCorr=0 #numero de correspondencias
            minNumCorr= 1234567890 #minimo numero de correspondencias
           
            try:
                numFP = (len(arrFPObjSource) + len(arrFPObjTarget))/2 #num medio de fp's
                jaccard = getJaccardIndex(arrPolyWKT[pol], arrPolyWKT[pol+1], correspondences[0])
                for j in range(len(correspondences[0])): #soma das correspondencias entre os 2 poligonos atuais
                    numCorr=numCorr+1
                    totalDistEntrePols = totalDistEntrePols+math.dist(correspondences[0][j][0], correspondences[0][j][1])

                if jaccard > minJaccard: #se o jaccard for menor que o minimo jaccard anterior
                    if (jaccard-minJaccard < 0.1): #se a diferença for ao nível das centesimas, entao verifica o numCorr
                        if numCorr/numFP < minNumCorr: #se o valor por correspondencia deste jaccard for menor que o anterior
                            minJaccard=jaccard #o jaccard passa a ser o novo minimo
                            minNumCorr = numCorr/numFP
                            distFinal=distancia
                            angFinal=angulo
                        else: 
                            minJaccard=minJaccard #o minimo jaccard anterior continua a ser o minimo
                            
                    else:
                        minJaccard = jaccard #o jaccard passa a ser o novo minimo
                        minNumCorr = numCorr/numFP #o valor medio por correspondencia passa a ser o novo minimo
                        distFinal=distancia
                        angFinal=angulo

                    f.write(f"min dist: {distancia} max angle: {angulo} valor medio por corr: {minNumCorr} jaccard: {minJaccard} \n")#so escreve isto se for um minimo
                
                f.flush()
            except ZeroDivisionError:
                f.write(f"Not enough correspondences. Skipping to next cycle. \n")

    valuesTiago.append([distFinal,angFinal])
    f.write(f"values: {valuesTiago}\n")
    f.write(f"--------------- NEW POLYGON ---------------\n")