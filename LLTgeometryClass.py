# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:44:27 2020

@author: n.tatomir
"""
import numpy as np
import math
import sys

class wingGeom():
    
    def __init__(self, b, Lam, beta, cRoot, cTip, Nspan, alpha, wakeLenFactor, HSbreak, typeSpacing, typeEvalPt):
        self.b = b
        self.Lam = Lam
        self.beta = beta
        self.cRoot = cRoot
        self.cTip = cTip
        self.Nspan = Nspan
        self.alpha = alpha
        self.wakeLenFactor = wakeLenFactor
        self.HSbreak = HSbreak
        self.typeSpacing = typeSpacing
        self.typeEvalPt = typeEvalPt
       
        self.theta = np.linspace(0,180,int(self.Nspan/2+1))
        
        #------PRIVATE-ATTRIBUTES---------------------------------------------#
        self.__wingGeomPts = self.__getWingGeomPts()
        self.__wingChordwisePts = self.__getWingChordwisePts()
        self.__wingDistributionPts = self.__getWingDistributionPts()
        self.__groupedHorseShoePts = self.__getGroupedHorseShoePts()
        
        #------PUBLIC-ATTRIBUTES----------------------------------------------#
        self.panelPts = self.getGroupedPanelPts()
        self.panelNormVec = self.getPanelNormalVec()
        self.panelArea = self.getPanelArea()
        self.boundVortexCenterPts = self.getBoundVortexCenterPts()
        self.evaluationPts = self.getEvaluationPts()
        self.horseShoeVortexPts = self.getHorseShoeVortexPts()

#---------------PRIVATE-FUNCTIONS---------------------------------------------#
        
    def __getWingGeomPts(self):
        wingGeomPts = np.zeros((7,3))
        wingGeomPts[0,:] =  0.5*self.b*math.tan( math.radians(self.Lam) ),           -0.5*self.b, 0.5*self.b*math.tan( math.radians(self.beta) )
        wingGeomPts[1,:] =  0. , 0. , 0.
        wingGeomPts[2,:] =  wingGeomPts[0,:] + [0, self.b, 0]
        wingGeomPts[3,:] =  wingGeomPts[2,:] + [self.cTip, 0, 0] 
        wingGeomPts[4,:] =  self.cRoot, 0., 0.
        wingGeomPts[5,:] =  wingGeomPts[3,:] + [0, -self.b, 0]
        wingGeomPts[6,:] =  wingGeomPts[0,:]
        return wingGeomPts # returns an ARRAY
    
    
    def __getWingChordwisePts(self):
        wingChordwisePts = np.zeros((2,6,3))
        
        wingChordwisePts[:,0,:] =  self.__wingGeomPts[0,:], self.__wingGeomPts[1,:]
        wingChordwisePts[:,1,:] =  self.__wingGeomPts[0,:] + [0.25*self.cTip,0,0], self.__wingGeomPts[1,:] + [0.25*self.cRoot,0,0]
        wingChordwisePts[:,2,:] =  self.__wingGeomPts[0,:] + [0.75*self.cTip,0,0], self.__wingGeomPts[1,:] + [0.75*self.cRoot,0,0]
        wingChordwisePts[:,3,:] =  self.__wingGeomPts[0,:] + [self.cTip,0,0], self.__wingGeomPts[1,:] + [self.cRoot,0,0]
        wingChordwisePts[:,4,:] =  wingChordwisePts[0,3,:] + [self.HSbreak*self.cTip,0,0] , wingChordwisePts[1,3,:] + [self.HSbreak*self.cRoot,0,0]
        wingChordwisePts[:,5,:] =  wingChordwisePts[0,4,:] + [self.wakeLenFactor*self.b , 0 , self.wakeLenFactor*self.b*math.tan(math.radians(self.alpha))] , wingChordwisePts[1,4,:] + [ self.wakeLenFactor*self.b , 0 , self.wakeLenFactor*self.b*math.tan(math.radians(self.alpha)) ]
        
        return wingChordwisePts
    
    
    def __getWingDistributionPts(self):
        
        wingDistributionPts = np.zeros((6, self.Nspan+1, 3))
        
        for i in range(6):
            pt0 = self.__wingChordwisePts[0,i,:]
            pt1 = self.__wingChordwisePts[1,i,:]
            vec = pt1 - pt0
                            
            if self.typeSpacing == "uniform":
                for j in range(self.Nspan+1):
                    if j <= int(self.Nspan/2):
                        f = 1/(self.Nspan/2)
                        wingDistributionPts[i,j,:] = pt0 + j*f*vec
                    else:
                        wingDistributionPts[i,j,:] = wingDistributionPts[i,self.Nspan-j,0], -wingDistributionPts[i,self.Nspan-j,1], wingDistributionPts[i,self.Nspan-j,2]
                    
            elif self.typeSpacing == "cosine":
                for j in range(self.Nspan+1):
                    if j <= int(self.Nspan/2):
                        wingDistributionPts[i,j,:] = pt0 + 0.5*vec*(1-math.cos(math.radians(self.theta[j])))
                    else:
                        wingDistributionPts[i,j,:] = wingDistributionPts[i,self.Nspan-j,0], -wingDistributionPts[i,self.Nspan-j,1], wingDistributionPts[i,self.Nspan-j,2]
            else:
                sys.exit('WARNING: typeSpacing must be either "uniform" or "cosine"!!!')
        return wingDistributionPts
    
    
    def __getGroupedHorseShoePts(self):
        groupedHorseShoePts = np.zeros((3, self.Nspan +1 , 3))
        
        groupedHorseShoePts[0,:,:] = self.__wingDistributionPts[1,:,:]
        groupedHorseShoePts[1,:,:] = self.__wingDistributionPts[4,:,:]
        groupedHorseShoePts[2,:,:] = self.__wingDistributionPts[5,:,:]
        
        return groupedHorseShoePts


#-----------------------PUBLIC-FUNCTIONS-----------------------------------------#    
    
    def getGroupedPanelPts(self):
        panelPts = np.zeros((2,self.Nspan+1,3))
        
        panelPts[0,:,:] = self.__wingDistributionPts[0,:,:]
        panelPts[1,:,:] = self.__wingDistributionPts[3,:,:]

        return panelPts   
    
    
    def getPanelNormalVec(self):
        panelNormVec = np.zeros((self.Nspan, 3))
        
        for i in range(self.Nspan):
            Ak = self.panelPts[1,i+1,:] - self.panelPts[0,i,:] 
            Bk = self.panelPts[0,i+1,:] - self.panelPts[1,i,:]
            panelNormVec[i,:] =  np.cross(Ak,Bk)/np.linalg.norm(np.cross(Ak,Bk))
        
        return panelNormVec
    
    
    def getPanelArea(self):
        panelArea = np.zeros(self.Nspan)
        
        for i in range(self.Nspan):
            Ak = self.panelPts[1,i+1,:] - self.panelPts[0,i,:] 
            Bk = self.panelPts[0,i+1,:] - self.panelPts[1,i,:]
            panelArea[i] =  0.5*np.linalg.norm(np.cross(Ak,Bk))
        
        return panelArea
    
    
    def getBoundVortexCenterPts(self):
        boundVortexCenterPts = np.zeros((self.Nspan,3))
        
        for i in range(self.Nspan):
            boundVortexCenterPts[i,:] = self.__groupedHorseShoePts[0,i,:] + 0.5*(self.__groupedHorseShoePts[0,i+1,:]-self.__groupedHorseShoePts[0,i,:]) 
      
        return boundVortexCenterPts             
            
    
    def getEvaluationPts(self):
        evaluationPts = np.zeros((self.Nspan,3))
        
        if self.typeSpacing == "uniform":
            if self.typeEvalPt == "center":
                for i in range(self.Nspan):
                    evaluationPts[i,:] = self.__wingDistributionPts[2,i,:] + 0.5*(self.__wingDistributionPts[2,i+1,:]-self.__wingDistributionPts[2,i,:]) 
            else:
                sys.exit('WARNING: If typeSpacing = "uniform", then typeEvalPt must be equal to "center"!!!')
        elif self.typeSpacing == "cosine":
            if self.typeEvalPt == "center":
                for i in range(self.Nspan):
                    evaluationPts[i,:] = self.__wingDistributionPts[2,i,:] + 0.5*(self.__wingDistributionPts[2,i+1,:]-self.__wingDistributionPts[2,i,:]) 
            elif self.typeEvalPt == "glauert":
                for i in range(self.Nspan):
                    pt0 = self.__wingDistributionPts[2,0,:]
                    pt1 = self.__wingDistributionPts[2,int(self.Nspan/2),:]
                    vec = pt1 - pt0
                    if i <= int(self.Nspan/2 - 1):
                        evaluationPts[i,:] = pt0 + 0.5*vec*(1-math.cos(math.radians((self.theta[i+1]+self.theta[i])/2)))
                    else:
                        evaluationPts[i,:] = evaluationPts[(self.Nspan-1)-i,0], -evaluationPts[(self.Nspan-1)-i,1], evaluationPts[(self.Nspan-1)-i,2]
            else:
                sys.exit('WARNING: typeEvalPt must be "center" or "glauert"!!!')
        
        return evaluationPts
    
    
    def getHorseShoeVortexPts(self):
        horseShoeVortexPts = np.zeros((self.Nspan, 6, 3))
        
        for i in range(self.Nspan):
            horseShoeVortexPts[i,0,:] = self.__groupedHorseShoePts[2,i,:]
            horseShoeVortexPts[i,1,:] = self.__groupedHorseShoePts[1,i,:]
            horseShoeVortexPts[i,2,:] = self.__groupedHorseShoePts[0,i,:]
            horseShoeVortexPts[i,3,:] = self.__groupedHorseShoePts[0,i+1,:]
            horseShoeVortexPts[i,4,:] = self.__groupedHorseShoePts[1,i+1,:]
            horseShoeVortexPts[i,5,:] = self.__groupedHorseShoePts[2,i+1,:]
        
        return horseShoeVortexPts
                
#----------------------------------------------------------------------------------------#

        
      