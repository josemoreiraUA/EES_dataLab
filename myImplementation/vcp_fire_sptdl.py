import numpy as np
from shapely.geometry import Polygon, LineString, Point
import shapely.geometry as sg
from shapely.ops import unary_union
import shapely.wkt
import geopandas as gpd
import matplotlib.pyplot as plt
import math
from numpy.linalg import eig
from ros import Ros
from featurepoint import FeaturePoint
from polygonCorrespondence import PolygonCorrespondence
import re
import time

#variaveis globais
arrPolyWKT = [] #array com os poligonos wkt
arrPolyGPD = [] #array com os poligonos gpd (para efeitos de visualizacao)

def get_angles(vec_1,vec_2): #obter o angulo entre 2 vetores (graus)
    
    dot = np.dot(vec_1, vec_2) 
    det = np.cross(vec_1,vec_2) 
    angle_in_rad = np.arctan2(det,dot) 
    return np.degrees(angle_in_rad)


def simplify_by_angle(poly_in, deg_tol): #usando trigonometria remover os pontos que formam angulos inferiores a um valor

    ext_poly_coords = poly_in.exterior.coords[:] #obtem as coordenadas exteriores do poligono (vertices que definem a forma)
    vector_rep = np.diff(ext_poly_coords, axis=0) #calcula a diferenca entre as coordenadas => (436, 33) (436, 34) = (0, 1)
    num_vectors = len(vector_rep) #num vetores existentes = numero de diferencas 
    angles_list = []
    for i in range(0, num_vectors): #para cada vetor
        angles_list.append(np.abs(get_angles(vector_rep[i], vector_rep[(i + 1) % num_vectors]))) #calcula o angulo entre os vetores

    #   get mask satisfying tolerance

    thresh_vals_by_deg = np.where(np.array(angles_list) > deg_tol) #guarda o indice dos vetores que formam angulos superiores ao angulo dado
    new_idx = list(thresh_vals_by_deg[0] + 1)
    new_vertices = [ext_poly_coords[idx] for idx in new_idx] #vai buscar ao array de vertices apenas os indices que correspondem a angulos maiores que o dado

    return Polygon(new_vertices)

def readWKT(wktFile): #le os ficheiros wkt e guarda os poligonos num array

    f = open(wktFile) 
    wktCoords = f.read() #wktCoords fica com as coordenadas do poligono

    wktPol = shapely.wkt.loads(wktCoords) #cria um poligono wkt
    wktPolSimp = simplify_by_angle(wktPol, 45) #simplifica o poligono (usar quando for testado com poligonos grandes)
    arrPolyWKT.append(wktPol) #acrescenta o poligono do ficheiro ao array de poligonos (deve acrescentar o wktPolSimp quando forem poligonos grandes)
    gpdPol = gpd.GeoSeries([wktPol]) #trocar para wktPolSimp quando usar o meu algoritmo de simplificacao 
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

        angleABCc = math.acos(valueABC) #arccos para calculo do angulo alpha´
        #fazer print dos angulos para ver se ha negativos
        angleABC = math.degrees(angleABCc) #150 < angleABC < 360-angleABC
        #|| angleABC >= 360-angleABC
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
    #print("CORRESPONDENCE BETWEEN ONLY FEATURE POINTS")
    #print("TARGET POINTS")
    #for point in polSrc.exterior.coords:
        #print(point)
    #print("SOURCE POINTS")
    #for point in polTrg.exterior.coords:
        #print(point)
    print("PTS SRC", len(polSrc.exterior.coords))
    print("PTS TRG", len(polTrg.exterior.coords))

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
        #print(arrCorrFP[i])
        #print(arrCorrFP[i][0])
        #print(arrCorrFP[i][1])

        if arrCorrFP[i][0] == arrCorrFP[i-1][0]: #se o ponto source for igual ao anterior
            arrCorrFPtmp.remove(arrCorrFP[i]) #remove essa correspondencia, e os pontos passam a ser tratados como pontos normais (e nao FP)
        if arrCorrFP[i][1] == arrCorrFP[i-1][1]: #se o ponto target for igual ao anterior
            arrCorrFPtmp.remove(arrCorrFP[i]) #remove essa correspondencia, e os pontos passam a ser tratados como pontos normais (e nao FP)

    arrCorrFP = arrCorrFPtmp.copy() #o array de correspondencias passa a ser o tmp
    #for corr in arrCorrFP:
        #print(corr)
    #print(len(polSrc.exterior.coords))
    #print(len(polTrg.exterior.coords))

    arrCorrP = [] #array onde vao ficar as correspondencias entre pontos
    
    arrSrcPoints = polSrc.exterior.coords[:]
    arrTrgPoints = polTrg.exterior.coords[:] #aqui se calhar temos que adicionar os pontos que sao repetidos
    for corr in range(len(arrCorrFP)-1): #para cada correspondencia

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
                            #arrCorrP.append([s,arrCorrFP[corr][1]])
                            arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([s, arrTrg[corrIdx]])
                else:
                    fct = len(idxSource) / len(idxTarget)
                    for i, s in enumerate(arrTrg):
                        if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                            #arrCorrP.append([arrCorrFP[corr][0],s])
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
            fator = (len(arrS)) / (len(arrT)) #garantir q nao da zero
            for i, s in enumerate(arrT): #para cada ponto no array de pontos nao fp
                if fator==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                    #arrCorrP.append([arrCorrFP[corr][0],s]) #o fp ponto atual do source fica com o ponto do target
                    arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                    break
                else:
                    corrIdx = int(i * fator) #a correspondencia e' feita com
                    arrCorrP.append([arrS[corrIdx],s]) #e' adicionado ao array de correspondencias entre pontos a correspond 
                        
        else:
            fator = len(arrT) / len(arrS)
            for i, s in enumerate(arrS):
                if fator==0: #se o fator for 0, isto e', se o arrT=0 (nao houver pontos entre os 2 fp do target)
                    #arrCorrP.append([s,arrCorrFP[corr][1]]) #o fp ponto atual do target fica com o ponto do source
                    arrCorrP = divide_line_segment_trg(arrCorrFP[corr][1],arrCorrFP[corr+1][1],len(arrS),arrS,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                    break
                else:
                    corrIdx = int(i * fator)
                    arrCorrP.append([s, arrT[corrIdx]])
        
        if(corr == (len(arrCorrFP) - 2)): #se ja estiver na penultima correspondencia
            arrCorrP.append(arrCorrFP[corr+1]) #adiciona a ultima correspondencia entre fp (isto tem que se ver melhor porque se houver mais pontos depois desse feature point nao sao incluidos nas novas correspondencias)

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
                            #arrCorrP.append([s,arrCorrFP[corr][1]])
                            arrCorrP = divide_line_segment_src(arrCorrFP[corr][0],arrCorrFP[corr+1][0],len(arrT),arrT,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([s, arrTrg[corrIdx]])
                else:
                    fct = len(idxSource) / len(idxTarget)
                    for i, s in enumerate(arrTrg):
                        if fct==0: #se o fator for 0, isto e', se o arrS=0 (nao houver pontos entre os 2 fp do source)
                            #arrCorrP.append([arrCorrFP[corr][0],s])
                            arrCorrP = divide_line_segment_trg(arrCorrFP[corr][1],arrCorrFP[corr+1][1],len(arrS),arrS,arrCorrP) #ponto inicial, ponto final, numero de segmentos, pontos para correspondencia, array de correspondencias
                            break
                        else:
                            corrIdx = int(i * fct)
                            arrCorrP.append([arrSrc[corrIdx],s])

    #print("CORRESPONDENCE BETWEEN ALL POINTS")
    #for corr in arrCorrP:
        #print(corr)

    return arrCorrP    

def renderCorrespondences(arrOfCorrespondences):
    #print("ARR CORR: ", arrOfCorrespondences)
    fig1, ax1 = plt.subplots()
    arrOfPoints = [[] for _ in range(len(arrOfCorrespondences))] #arrOfPoints e' um array tridimensional (transicao entre poligonos -> correspondencias -> pontos )
    #print(arrOfCorrespondences)
    for i in range(len(arrOfCorrespondences)): #para cada par de poligonos (correspondencias entre poligonos)
        for j in range(len(arrOfCorrespondences[i])): #para cada correspondencia
            arr=[] #a cada novo conjunto de correspondencias recomeca
            ptSrc = arrOfCorrespondences[i][j][0]
            ptTrg = arrOfCorrespondences[i][j][1]
            #print(ptSrc, ptTrg)
            pointDist =  math.dist(ptSrc, ptTrg) #calcula a diferenca entre o ponto atual e o proximo
            arr.append(ptSrc) #adiciona o ponto source
            for k in range(1, 31): #divide a diferenca em 30 partes, para fazer a interpolacao
                #if pointDist >= 0.0:
                coordX = ptSrc[0] + (k * (ptTrg[0] - ptSrc[0]) / 31)
                coordY = ptSrc[1] + (k * (ptTrg[1] - ptSrc[1]) / 31)
                arr.append((coordX, coordY)) #adiciona os pontos intermedios
            arr.append(ptTrg) #adiciona o ponto target
            arrOfPoints[i].append(arr) #cada conjunto de pontos e' uma posicao no arrOfPoints
    
    arrOfPointsForPolygon = []
    arrOfPolygons = []
    
    for k in range(len(arrOfPoints)): #para cada conjunto de pontos entre poligonos
        for j in range(len(arrOfPoints[0][0])): #para cada conjunto de pontos intermedios (sempre 32) 
            arrOfPointsForPolygon = [] #reset do array a cada novo poligono
            for i in range(len(arrOfPoints[k])): #acede a mesma posicao para todos os conjuntos, e os pontos que ai estao vao formar um poligono
                #print("k,i,j", k,i,j)
                #print(arrOfPoints[k][i][j])
                arrOfPointsForPolygon.append(arrOfPoints[k][i][j]) #acede ao ponto j, da correspondencia i, na relacao entre poligonos k
                #esta a quebrar aqui, na transicao do t7 para o t8
            poly = Polygon([p[0], p[1]] for p in arrOfPointsForPolygon) #para cada conjunto de pontos cria um poligono
            arrOfPolygons.append(poly) #acrescenta a um array de poligonos'''
    #print(i)
    
    source = arrOfPolygons[0]
    target = arrOfPolygons[31]
    xxOuterSrc,yyOuterSrc = source.exterior.xy
    xxOuterTrg,yyOuterTrg = target.exterior.xy

    cnt=0    
    for polygon in arrOfPolygons: #para cada poligono
        plt.plot(xxOuterSrc, yyOuterSrc, color="green") #e mostra no grafico
        plt.plot(xxOuterTrg, yyOuterTrg, color="red") #e mostra no grafico
        polyGPDSimp = gpd.GeoSeries(polygon) #transforma em gpd
        polyGPDSimp.plot(ax=ax1, color="blue") #e mostra no grafico (nao precisa do alpha)
        plt.pause(0.2) #durante meio segundo
        cnt+=1
        if cnt%32==0:
            xxOuterSrc=xxOuterTrg
            yyOuterSrc=yyOuterTrg
            if  cnt+32 < len(arrOfPolygons):
                target = arrOfPolygons[cnt+32]
            else:
                target = arrOfPolygons[-1]
            xxOuterTrg,yyOuterTrg = target.exterior.xy
            plt.pause(1.2) #durante meio segundo
        plt.clf() #apaga
        ax1 = plt.gca() #para evitar o arrasto da imagem

    
    source = arrOfPolygons[0]
    target = arrOfPolygons[31]

    cnt=0    
    for x, y in source.exterior.coords:
        plt.scatter(x, y, color='green')

    for x, y in target.exterior.coords:
        plt.scatter(x, y, color='red')

    source=target
    for corr in arrOfCorrespondences[0]:
        line = sg.LineString(corr)
        plt.plot(*line.xy, color="blue", alpha=0.5)
    plt.show()

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
    #print(jaccardIndexSM)

    #calculo do jaccard index do middle para o target
    intersectionMT = midPoly.intersection(target).area
    unionMT = unary_union([midPoly, target]).area
    jaccardIndexMT = intersectionMT/unionMT
    #print(jaccardIndexMT)
    
    jaccardIndex = (jaccardIndexSM + jaccardIndexMT) / 2
    print(jaccardIndex)

    #sempre a dar 0...
    return jaccardIndex

arrPontos = []
wktFiles = ["f_spt_dl0.wkt","f_spt_dl1.wkt","f_spt_dl2.wkt","f_spt_dl3.wkt","f_spt_dl4.wkt","f_spt_dl5.wkt","f_spt_dl6.wkt","f_spt_dl7.wkt","f_spt_dl8.wkt"]
valuesTiago=[]
#values = [[1, 130], [6, 120], [7, 150], [10, 180], [10, 180], [1, 160], [8, 160], [9, 170]]
correspondences = []
correspondencesBetweenAllPoints = [] #correspondencias entre todos os pontos, nao so feature points
jaccard = 0 #jaccard index medio entre S->M e M->T

for file in wktFiles:
    readWKT(file) #popula o array de poligonos (arrPolyWKT)

f=open('f_spt_dl.txt', 'w')
for pol in range(len(arrPolyWKT)-1): #para cada par de poligonos
    minJaccard = 0 
    for angulo in range(120,181,10): #testa para valores de ang  
        for distancia in range(1,20): #e de dist

            f.write(f"pol: {pol}, ang: {angulo}, dist: {distancia} \n")

            try:
                arrFPObjSource = [] #a cada poligono os feature points sao resetados
                arrFPObjTarget = []
                correspondences = []

                if((angulo==170 or angulo==180) and distancia<3):
                    f.write(f"Too many values\n")
                    continue

                start_time = time.time()
                #featurePointsSource = getFeaturePoints(arrPolyWKT[pol], valuesTiago[pol][0], valuesTiago[pol][1]) #array de feature points do poligono origem (5,180) para poligonos grandes 
                featurePointsSource = getFeaturePoints(arrPolyWKT[pol], distancia, angulo)
                #featurePointsTarget = getFeaturePoints(arrPolyWKT[pol+1], valuesTiago[pol][0], valuesTiago[pol][1]) #array de feature points do poligono destino
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
                
                #o fp.featSize por vezes da negativo porque o mini > valueRouL
                for fp in arrFPObjSource: #calcula as features de cada FP do source
                    fp.rolSizeNorm = (fp.rolSize-mini) / (maxi - mini)
                    fp.rorSizeNorm = (fp.rorSize-mini) / (maxi - mini)
                    fp.getFeatures(fp.ros, arrPolyWKT[pol].length)

                for fp in arrFPObjTarget: #calcula as features de cada FP do target
                    fp.rolSizeNorm = (fp.rolSize-miniTarget) / (maxiTarget - miniTarget)
                    fp.rorSizeNorm = (fp.rorSize-miniTarget) / (maxiTarget - miniTarget)
                    fp.getFeatures(fp.ros, arrPolyWKT[pol+1].length)
                
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

                    f.write(f"min dist: {distancia} max angle: {angulo} valor medio por corr: {minNumCorr} jaccard: {minJaccard} \n")
                
                f.flush()
            except ZeroDivisionError:
                f.write(f"Not enough correspondences. Skipping to next cycle. \n")

    valuesTiago.append([distFinal,angFinal])
    f.write(f"values: {valuesTiago}\n")
    f.write(f"--------------- NEW POLYGON ---------------\n")
#renderCorrespondences(correspondences) #chama a funcao para mostrar a evolucao dos poligonos 