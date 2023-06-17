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
correspondences = [[[(181.952, 9.0), (224.0, 19.0)], [(135.954, 12.001), (189.948, 9.0)], [(131.953, 15.0), (136.948, 12.001)], [(116.954, 13.001), (76.947, 18.0)], [(64.0, 22.001), (67.947, 22.001)], [(0.952, 16.001), (0.947, 17.001)], [(0.952, 197.001), (0.947, 203.0)], [(56.72, 208.422), (75.627, 208.422)], [(117.952, 218.0), (75.627, 208.422)], [(148.954, 227.0), (151.249, 227.37)], [(241.954, 186.0), (189.065, 227.37)], [(283.602, 202.736), (243.937, 215.005)], [(316.687, 215.057), (301.947, 224.007)], [(349.772, 254.844), (331.947, 213.0)], [(402.954, 250.012), (354.947, 251.005)], [(402.954, 250.012), (453.948, 256.005)], [(477.358, 240.637), (521.786, 251.06)], [(554.882, 240.637), (567.164, 265.284)], [(553.935, 267.192), (567.166, 284.253)], [(600.261, 274.752), (610.654, 293.688)], [(718.431, 284.243), (726.934, 293.688)], [(790.283, 271.908), (812.962, 284.247)], [(838.5, 252.947), (846.08, 273.824)], [(875.366, 224.531), (873.468, 238.738)], [(964.244, 208.428), (978.935, 216.0)], [(1014.953, 156.001), (1065.387, 160.107)], [(1024.0, 130.001), (1105.088, 162.0)], [(1110.953, 149.0), (1109.816, 120.317)], [(1141.966, 124.107), (1159.914, 113.684)], [(1173.16, 99.477), (1246.882, 100.424)], [(1139.126, 72.947), (1155.185, 83.37)], [(1085.247, 64.424), (1155.185, 83.37)], [(993.953, 51.0), (1103.948, 76.001)], [(850.954, 51.0), (1016.947, 65.0)], [(810.952, 48.001), (941.947, 59.0)], [(787.954, 53.001), (921.947, 59.0)], [(788.954, 56.0), (915.948, 63.0)], [(774.953, 75.001), (889.947, 82.0)], [(769.953, 75.001), (854.947, 70.001)], [(764.954, 71.001), (796.948, 63.0)], [(637.954, 58.001), (755.948, 85.001)], [(633.953, 63.0), (708.947, 64.0)], [(555.952, 65.0), (708.947, 64.0)], [(549.953, 62.001), (708.947, 58.001)], [(502.953, 51.0), (689.948, 59.0)], [(498.954, 50.0), (559.948, 61.001)], [(346.953, 32.0), (460.947, 43.001)], [(320.952, 34.001), (424.948, 48.001)], [(250.953, 28.0), (258.947, 23.0)]], [[(224.0, 19.0), (120.987, 19.0)], [(189.948, 9.0), (120.987, 19.0)], [(136.948, 12.001), (120.987, 19.0)], [(76.947, 18.0), (0.986, 12.001)], [(67.947, 22.001), (0.986, 12.001)], [(0.947, 17.001), (0.986, 12.001)], [(0.947, 203.0), (0.986, 194.0)], [(75.627, 208.422), (0.986, 194.0)], [(151.249, 227.37), (85.986, 214.012)], [(151.249, 227.37), (128.986, 215.005)], [(189.065, 227.37), (188.124, 234.947)], [(243.937, 215.005), (338.988, 255.005)], [(301.947, 224.007), (352.0, 246.005)], [(301.947, 224.007), (407.987, 242.005)], [(331.947, 213.0), (550.988, 238.005)], [(331.947, 213.0), (548.987, 244.0)], [(354.947, 251.005), (605.02, 303.16)], [(453.948, 256.005), (750.599, 300.322)], [(521.786, 251.06), (853.644, 289.899)], [(567.164, 265.284), (939.987, 228.005)], [(567.164, 265.284), (951.987, 204.0)], [(567.166, 284.253), (964.987, 190.0)], [(610.654, 293.688), (1048.987, 178.001)], [(726.934, 293.688), (1048.987, 178.001)], [(812.962, 284.247), (1142.987, 144.0)], [(846.08, 273.824), (1200.987, 102.001)], [(846.08, 273.824), (1199.987, 69.0)], [(859.774, 256.281), (1166.987, 65.0)], [(873.468, 238.738), (1160.988, 52.001)], [(978.935, 216.0), (912.987, 53.001)], [(1065.387, 160.107), (900.987, 57.001)]], [[(120.987, 19.0), (133.993, 28.001)], [(120.987, 19.0), (72.992, 29.0)], [(120.987, 19.0), (32.993, 19.0)], [(120.987, 19.0), (24.992, 21.0)], [(0.986, 12.001), (0.993, 22.001)], [(0.986, 194.0), (0.993, 203.001)], [(85.986, 214.012), (72.992, 216.0)], [(85.986, 214.012), (149.993, 222.011)], [(128.986, 215.005), (209.992, 240.0)], [(128.986, 215.005), (283.993, 224.0)], [(188.124, 234.947), (288.0, 220.009)], [(338.988, 255.005), (288.0, 218.005)], [(338.988, 255.005), (321.992, 208.005)], [(352.0, 246.005), (333.994, 214.0)], [(352.0, 246.005), (338.993, 244.005)], [(407.987, 242.005), (392.992, 263.007)], [(550.988, 238.005), (453.993, 255.0)], [(550.988, 238.005), (492.991, 238.007)], [(548.987, 244.0), (611.343, 338.038)], [(605.02, 303.16), (893.132, 311.302)], [(750.599, 300.322), (942.993, 235.0)], [(853.644, 289.899), (942.993, 235.0)], [(939.987, 228.005), (942.993, 235.0)], [(951.987, 204.0), (942.993, 235.0)], [(964.987, 190.0), (1102.993, 184.001)], [(1048.987, 178.001), (1102.993, 184.001)], [(1142.987, 144.0), (1102.993, 184.001)], [(1200.987, 102.001), (1193.993, 118.001)]]]
renderCorrespondences(correspondences) #chama a funcao para mostrar a evolucao dos poligonos 

#1 calculo dos params
#2 calculo das correspondencias com base nos params
#3 render das correspondencias