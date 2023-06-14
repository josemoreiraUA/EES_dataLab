import numpy as np
from correspondence import Correspondence

class PolygonCorrespondence:
    def __init__(self, arrFeaturePointsSource, arrFeaturePointsTarget, featVariation, featSideVariation, featSizeVariation, weightSimilarity):
        self.arrFeaturePointsSource = arrFeaturePointsSource
        self.arrFeaturePointsTarget = arrFeaturePointsTarget
        self.weightSimilarity = weightSimilarity

        #similarityCosts
        self.simCostsGraph = np.zeros((len(arrFeaturePointsSource), len(arrFeaturePointsTarget)))
        for s in range(len(arrFeaturePointsSource)):
            for t in range(len(arrFeaturePointsTarget)):
                self.simCostsGraph[s][t] = arrFeaturePointsSource[s].similarityCost(arrFeaturePointsTarget[t], featVariation, featSideVariation, featSizeVariation)

        #discardCosts
        self.discardCostSource = []
        self.discardCostTarget = []
        for i in range(len(arrFeaturePointsSource)):
            self.discardCostSource.append(arrFeaturePointsSource[i].discardCost(featVariation, featSideVariation, featSizeVariation))
        for j in range(len(arrFeaturePointsTarget)):
            self.discardCostTarget.append(arrFeaturePointsTarget[j].discardCost(featVariation, featSideVariation, featSizeVariation))

    def getFPCorrespondences(self, skips):
        min = 123456789
        correspondences = []

        for s_init in range(len(self.arrFeaturePointsSource)): #para cada ponto source
            
            for t_init in range(len(self.arrFeaturePointsTarget)): #para cada ponto target
                tmp = self.getFeaturePointCorrespondencesByPath(s_init, t_init, skips) #o tmp e' uma solucao
                total = 0 
                for c in tmp: #para cada correspondencia no tmp
                    total += c.cost #adiciona o seu custo

                if min > total: #se o custo dessa solucao (tmp) for menor que o min
                    min = total
                    correspondences = tmp #o tmp passa a ser a solucao

        #ORDER BY SOURCE INDEX
        correspondences.sort(key=lambda x: x.s_i)

        for i in range(len(correspondences)-1):
            if abs(correspondences[i].s_i - correspondences[i+1].s_i) > 1: #se a diferenca entre indices for > 1 quer dizer q nao sao consecutivos
                interval = abs(correspondences[i].s_i - correspondences[i+1].s_i)
                for j in range(1,interval):
                    c = Correspondence()
                    c.point_s = self.arrFeaturePointsSource[i+j]
                    c.point_s_coordinates = (self.arrFeaturePointsSource[i+j].point)
                    c.s_i = i+j
                    c.point_t = correspondences[i].point_t
                    c.t_i = correspondences[i].t_i
                    c.point_t_coordinates = (correspondences[i].point_t_coordinates)
                    correspondences.insert(i+j,c)

        return correspondences
    
    def getPath(self, s_init, t_init, skips):
        correspondences = []
        s_i, t_i = s_init, t_init
        delta_s, delta_t = 0, 0
        while True:
            # direct correspondence
            s_next, t_next = s_i, t_i #indices dos pontos do poligono
            min = self.delta(s_i, s_next, t_i, t_next)
            # skips
            nr_s_skips, nr_t_skips = 0, 0  # number of skips
            for s_skips in range(skips + 1): #para os skips de source
                for t_skips in range(skips + 1): #e de target
                    #calcula os custos entre o source+skip e o target+skip. % e' para dar a volta quando chega ao ultimo
                    tmp = self.delta(s_i, (s_i + s_skips) % len(self.arrFeaturePointsSource), t_i, (t_i + t_skips) % len(self.arrFeaturePointsTarget)) #calcula os custos entre [src, trg+2]
                    if tmp < min: #se o custo entre target e source for minimo, esse par fica na solucao final
                        min = tmp
                        nr_s_skips = s_skips
                        nr_t_skips = t_skips
                        s_next, t_next = (s_i + s_skips) % len(self.arrFeaturePointsSource), (t_i + t_skips) % len(self.arrFeaturePointsTarget)
            # increment already processed vertexes
            delta_s += nr_s_skips
            delta_t += nr_t_skips
            if delta_s > len(self.arrFeaturePointsSource) or delta_t > len(self.arrFeaturePointsTarget):
                break

            c = Correspondence()
            c.point_s = self.arrFeaturePointsSource[s_next]
            c.point_s_coordinates = (self.arrFeaturePointsSource[s_next].point)
            c.s_i = s_next
            c.point_t = self.arrFeaturePointsTarget[t_next]
            c.t_i = t_next
            c.point_t_coordinates = (self.arrFeaturePointsTarget[t_next].point) 
            c.cost = min
            correspondences.append(c)
            # move to next values
            s_i, t_i = (s_next + 1) % len(self.arrFeaturePointsSource), (t_next + 1) % len(self.arrFeaturePointsTarget)
            delta_s += 1
            delta_t += 1
            
            if delta_s >= len(self.arrFeaturePointsSource) or delta_t >= len(self.arrFeaturePointsTarget): #condicao de paragem
                break

        return correspondences
    
    def delta(self, s_begin, s_end, t_begin, t_end):
        s_begin, s_end = s_begin % len(self.arrFeaturePointsSource), s_end % len(self.arrFeaturePointsSource)
        t_begin, t_end = t_begin % len(self.arrFeaturePointsTarget), t_end % len(self.arrFeaturePointsTarget)
        s_dist, t_dist = abs(s_end - s_begin), abs(t_end - t_begin)
        costs = 0
        for i in range(s_dist):
            costs += self.discardCostSource[(s_begin + i) % len(self.arrFeaturePointsSource)]
        for j in range(t_dist):
            costs += self.discardCostTarget[(t_begin + j) % len(self.arrFeaturePointsTarget)]
        costs += self.weightSimilarity * self.simCostsGraph[s_end][t_end]
        return costs

    def getFeaturePointCorrespondencesByPath(self, s_init, t_init, skips):
        return self.getPath(s_init, t_init, skips)