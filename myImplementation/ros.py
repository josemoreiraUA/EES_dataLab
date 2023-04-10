import numpy as np
from shapely.geometry import Point
from numpy.linalg import eig
import math

class Ros:
    def __init__(self, points, featurePoint, polygon): #criar um objeto: ros1 = ros(arrFP, fp1)
        #inicializacao de variaveis
        self.points = points
        self.featurePoint = featurePoint
        self.polygon = polygon
        self.center = Point(0,0)
        self.bisector = [0 for i in range(2)] 
        self.covMatrix = [[0 for _ in range(2)] for _ in range(2)]
        self.eigenValue = 0
        self.eigenVector = None
        self.normalValue = 0
        self.tangentValue = 0
        self.normalVector = [[0 for _ in range(2)] for _ in range(2)]
        self.tangentVector = [[0 for _ in range(2)] for _ in range(2)]
        self.featureVariation = 0
        self.featureSideVariation = 0
        self.featureSize = 0
        #chamada das funcoes
        if len(self.points) < 3: #se for rol ou ror
            self.getThirdPointForRolr() #acrescenta um terceiro ponto (medio) ao array de pontos da region
        self.getCenter()
        self.getCovMatrix()
        self.getEigenValuesAndVectors()
        self.getNormalAndTangent()
        self.getBisector()

    '''def getCenter(self): #centro da ros
        rosCenterX = 0
        rosCenterY = 0
        cntPointInRos = 0

        for point in self.points: #para cada ponto da ros
            rosCenterX += point[0]  #faz a soma das coordenadas dos pontos
            rosCenterY += point[1] 
            cntPointInRos += 1
            
        rosCenterX = rosCenterX / (cntPointInRos+1) #e divide pelo total de pontos + 1
        rosCenterY = rosCenterY / (cntPointInRos+1)

        self.center = Point(rosCenterX, rosCenterY)'''
    
    def getCenter(self): ##nova versao (len(points), em vez de cntPointInRos+1)
        rosCenterX = 0
        rosCenterY = 0

        for point in self.points: #para cada ponto da ros
            rosCenterX += point[0]  #faz a soma das coordenadas dos pontos
            rosCenterY += point[1] 
            
        rosCenterX = rosCenterX / (len(self.points)) #e divide pelo total de pontos + 1
        rosCenterY = rosCenterY / (len(self.points))

        self.center = Point(rosCenterX, rosCenterY)

    '''def getCovMatrix(self):
        cntPointInRos = 0
        x = 0
        y = 0 

        for point in self.points: #calcula a covMatriz da ros
            
            x = point[0] - self.center.x
            y = point[1] - self.center.y
            cntPointInRos += 1

            self.covMatrix[0][0] += x * x
            self.covMatrix[0][1] += x * y
            self.covMatrix[1][0] += y * x
            self.covMatrix[1][1] += y * y

        self.covMatrix[0][0] /= cntPointInRos + 1
        self.covMatrix[0][1] /= cntPointInRos + 1
        self.covMatrix[1][0] /= cntPointInRos + 1
        self.covMatrix[1][1] /= cntPointInRos + 1'''
    
    def getCovMatrix(self): #nova versao
        x = 0
        y = 0 

        for point in self.points: #calcula a covMatriz da ros
            
            x = point[0] - self.center.x
            y = point[1] - self.center.y

            self.covMatrix[0][0] += x * x
            self.covMatrix[0][1] += x * y
            self.covMatrix[1][0] += y * x
            self.covMatrix[1][1] += y * y

        self.covMatrix[0][0] /= len(self.points)
        self.covMatrix[0][1] /= len(self.points)
        self.covMatrix[1][0] /= len(self.points)
        self.covMatrix[1][1] /= len(self.points)
    
    def getEigenValuesAndVectors(self):
        self.covMatrix = np.array(self.covMatrix)
        self.eigenValue, self.eigenVector = eig(self.covMatrix)

    def getBisector(self): 
        
        self.bisector[0] = self.points[0][0] + abs(self.points[0][0] - self.points[2][0]) /2
        self.bisector[1] = self.points[0][1] + abs(self.points[0][1] - self.points[2][1]) /2
        self.bisector[0] -= self.points[1][0]
        self.bisector[1] -= self.points[1][1]
        return self.bisector
    
    '''def getBisector(self): #nova funcao, com base nas formulas. nao esta a fazer diferenca no resultado final
        p1, p2, p3 = self.points
        v1 = np.array([p1[0] - p2[0], p1[1] - p2[1]]) 
        v2 = np.array([p3[0] - p2[0], p3[1] - p2[1]])
        v1 /= np.linalg.norm(v1) #normalizacao dos valores
        v2 /= np.linalg.norm(v2)
        dot_product = np.clip(np.dot(v1, v2), -1, 1) #garante que o valor esta entre 1 e -1, para o calculo do cosseno
        angle = math.acos(dot_product)
        bisector = v1 + v2
        bisector /= np.linalg.norm(bisector)
        bisector *= math.sqrt(2 * (1 + math.cos(angle)))
        bisector += np.array([p2[0], p2[1]])
        self.bisector = bisector.tolist()
        return self.bisector'''

    def getNormalAndTangent(self): #calcula o vetor normal e tangente
    
        #bisector = self.getBisector()
        #rosCenter = getRosCenter(ros)
        
        nX = self.eigenVector[0][0] + (self.points[1][0] - self.center.x) #ros[1] = feature point
        nY = self.eigenVector[0][1] + (self.points[1][1] - self.center.y)
        dot_ev1 = nX * self.bisector[0] + nY * self.bisector[1] #produto escalar entre o eigenvector1 e o bisector

        nX = self.eigenVector[1][0] + (self.points[1][0] - self.center.x)
        nY = self.eigenVector[1][1] + (self.points[1][1] - self.center.y)
        dot_ev2 = nX * self.bisector[0] + nY * self.bisector[1] #produto escalar entre o eigenvector2 e o bisector

        #definicao de qual dos eigen vector e' normal e tangente
        self.normalVector[0][0] = self.center.x
        self.normalVector[0][1] = self.center.y
        self.tangentVector[0][0] = self.center.x
        self.tangentVector[0][1] = self.center.y

        #o eigenvector com maior produto escalar e' o vetor normal
        if (dot_ev1 > dot_ev2):
            #eigen vector 1 => normal vector
            self.normalVector[1][0] = self.center.x + self.eigenVector[0][0]
            self.normalVector[1][1] = self.center.y + self.eigenVector[0][1]
            self.normalValue = self.eigenValue[0]

            #eigen vector 2 => tangent vector
            self.tangentVector[1][0] = self.center.x + self.eigenVector[1][0]
            self.tangentVector[1][1] = self.center.y + self.eigenVector[1][1]
            self.tangentValue = self.eigenValue[1]
        else:
            #eigen vector 2 => normal vector
            self.normalVector[1][0] = self.center.x + self.eigenVector[1][0]
            self.normalVector[1][1] = self.center.y + self.eigenVector[1][1]
            self.normalValue = self.eigenValue[1]

            #eigen vector 1 => tangent vector
            self.tangentVector[1][0] = self.center.x + self.eigenVector[0][0]
            self.tangentVector[1][1] = self.center.y + self.eigenVector[0][1]
            self.tangentValue = self.eigenValue[0]
        
    def getThirdPointForRolr(self): #acrescenta um terceiro ponto a region of left ou right

        pointInPolIdx = [] #index do ponto da region of left, no poligono
        middlePoint = () #ponto medio a acrescentar a rolr 
        middlePointX = 0
        middlePointY = 0

        for i in range(len(self.polygon.exterior.coords)): #para descobrir o indice dos feature points no poligono
            for point in self.points:
                if point == self.polygon.exterior.coords[i]:
                    pointInPolIdx.append(i)
        
        if(abs(pointInPolIdx[0]-pointInPolIdx[1])) > 1: #se houver pontos entre os feature points
            middlePoint = self.polygon.exterior.coords[int(((pointInPolIdx[0]+pointInPolIdx[1])/2))]  #o ponto medio e' o ponto do poligono com o indice medio
        else: #senao o ponto medio e' o ponto medio do segmento de reta, visto que como conecta dois pontos seguidos (sem pontos pelo meio) nao vai ser uma curva
            middlePointX = (self.points[0][0] + self.points[1][0])/2
            middlePointY = (self.points[0][1] + self.points[1][1])/2
            middlePoint = (middlePointX, middlePointY) #pode dar 0 (em ultimo caso usamos o indice seguinte ou anterior)
        if self.points[0] == self.featurePoint: #quer dizer que e' ror
            self.points.insert(0, middlePoint) #o ponto medio fica na 1a posicao para o fp continuar no meio
        elif self.points[1] == self.featurePoint: #quer dizer que e' rol
            self.points.append(middlePoint) #o ponto medio fica na ultima posicao para o fp continuar no meio
    