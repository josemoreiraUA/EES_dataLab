from shapely.geometry import Polygon
import geopandas as gpd
import matplotlib.pyplot as plt
import shapely.geometry as sg

def renderCorrespondences(arrOfCorrespondences): #funcao para fazer e mostrar a interpolacao dos poligonos

    fig1, ax1 = plt.subplots()
    arrOfPoints = [[] for _ in range(len(arrOfCorrespondences))] #arrOfPoints e' um array tridimensional (transicao entre poligonos -> correspondencias -> pontos )
    
    for i in range(len(arrOfCorrespondences)): #para cada par de poligonos (correspondencias entre poligonos)
        for j in range(len(arrOfCorrespondences[i])): #para cada correspondencia
            arr=[] #a cada novo conjunto de correspondencias recomeca
            ptSrc = arrOfCorrespondences[i][j][0]
            ptTrg = arrOfCorrespondences[i][j][1]
            arr.append(ptSrc) #adiciona o ponto source
            for k in range(1, 31): #divide a diferenca em 30 partes, para fazer a interpolacao
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
                arrOfPointsForPolygon.append(arrOfPoints[k][i][j]) #acede ao ponto j, da correspondencia i, na relacao entre poligonos k
            poly = Polygon([p[0], p[1]] for p in arrOfPointsForPolygon) #para cada conjunto de pontos cria um poligono
            arrOfPolygons.append(poly) #acrescenta a um array de poligonos'''
    
    #while True:
    source = arrOfPolygons[0]
    target = arrOfPolygons[31]
    xxOuterSrc,yyOuterSrc = source.exterior.xy
    xxOuterTrg,yyOuterTrg = target.exterior.xy

    cnt=0

    for polygon in arrOfPolygons:
        ax1.plot(xxOuterSrc, yyOuterSrc, color="green")  # Plot the green line
        ax1.plot(xxOuterTrg, yyOuterTrg, color="red")  # Plot the red line

        polyGPDSimp = gpd.GeoSeries(polygon)
        polyGPDSimp.plot(ax=ax1, color="blue")  # Plot the polygon

        plt.pause(0.2)

        cnt += 1
        if cnt % 32 == 0:
            xxOuterSrc = xxOuterTrg
            yyOuterSrc = yyOuterTrg

            if cnt + 32 < len(arrOfPolygons):
                target = arrOfPolygons[cnt + 32]
            else:
                target = arrOfPolygons[-1]
            
            xxOuterTrg, yyOuterTrg = target.exterior.xy

            plt.pause(0.5)

        ax1.clear()  # Clear the plot

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

#array de correspondencias que vem do programa vcp_fire (copiado do ficheiro txt gerado)
#correspondencesJaccDist = [] #insert correspondences in order to test
#renderCorrespondences(correspondencesJaccDist) #chama a funcao para mostrar a evolucao dos poligonos 