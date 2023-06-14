from ros import Ros

class FeaturePoint:
    def __init__(self, point, polygon, ros, rol, ror):
        self.point = point
        self.ros = ros
        self.rol = rol
        self.ror = ror
        self.polygon = polygon
        self.prev = None
        self.next = None
        self.featVariation = 0
        self.featVariationL = 0
        self.featVariationR = 0
        self.featSideVariation = 0
        self.featSize = 0
        self.normVal = 0
        self.normVect = None
        self.tangVal = 0
        self.tangVect = None
        self.rolSize = 0
        self.rorSize = 0
        self.rolSizeNorm = 0
        self.rorSizeNorm = 0
        self.maxSize = 0
        self.minSize = 0
    
    def getFeatures(self, arrPointsInROS):
        rosObj = Ros(arrPointsInROS, self.point, self.polygon)
        self.normVal = rosObj.normalValue
        self.tangVal = rosObj.tangentValue
        ex = 0

        bx = abs(arrPointsInROS[1][0]-arrPointsInROS[0][0])
        by = abs(arrPointsInROS[1][1]-arrPointsInROS[0][1])
        cx = abs(arrPointsInROS[0][0]-arrPointsInROS[2][0])
        cy = abs(arrPointsInROS[0][1]-arrPointsInROS[2][1])

        if (bx*cy - by*cx) >= 0: #feature point e' convexo
            ex = -1 
        else:
            ex = 1

        #featVariation
        self.featVariation = (ex*self.normVal)/(self.normVal+self.tangVal)
        #featSideVariation
        rol = Ros(arrPointsInROS[:2], self.point, self.polygon)
        self.featVariationL = abs(rol.normalValue / (rol.normalValue+rol.tangentValue))
        ror = Ros(arrPointsInROS[1:], self.point, self.polygon)
        self.featVariationR = abs(ror.normalValue / (ror.normalValue+ror.tangentValue))
        self.featSideVariation = (self.featVariationR + self.featVariationL) / 2
        #featSize
        self.featSize = ((self.rolSizeNorm + self.rorSizeNorm) / 2) #normalizacao
        
    def similarityCost(self, featPoint, wVariation, wSideVariation, wSizeVariation):

        deltaVariation = wVariation*abs(self.featVariation-featPoint.featVariation)
        deltaSideVariation = wSideVariation*((abs(self.featVariationL-featPoint.featVariationL))+(abs(self.featVariationR-featPoint.featVariationR)))/2
        deltaSize = wSizeVariation*((abs(self.featSize-featPoint.featSize)))/2

        if self.featSize > featPoint.featSize:
            max = self.featSize
        else:
            max = featPoint.featSize
    
        cost = max * (deltaVariation + deltaSideVariation + deltaSize)

        return cost

    def discardCost(self, wVariation, wSideVariation, wSizeVariation):
        omegaSize = self.featSize
        self.featVariation = wVariation*abs(self.featVariation)
        self.featSideVariation = wSideVariation*abs(self.featSideVariation)
        self.featSize = wSizeVariation*abs(self.featSize)

        cost = omegaSize * (self.featVariation + self.featSideVariation + self.featSize)

        return cost
    
    
    