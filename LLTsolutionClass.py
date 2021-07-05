# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:01:45 2020

@author: n.tatomir
"""
import numpy as np
import math

class wingCalc():
    
    def __init__(self,geom,freestreamVel, rho, cRef):
        self.geom = geom
        
        #----------added-attributes-------------------------------------------#
        self.freestreamVel = freestreamVel
        self.rho = rho
        self.cRef = cRef
        
        #------------calculated-attributes------------------------------------#
        self.qInf = 0.5*self.rho*(self.freestreamVel**2)
        self.QinfVec = np.array([self.freestreamVel*math.cos(math.radians(self.geom.alpha)), 0., self.freestreamVel*math.sin(math.radians(self.geom.alpha)) ])
        self.S_actual = sum(self.geom.panelArea)
        self.S_projected = sum(self.geom.panelArea)*math.cos(math.radians(self.geom.beta))

#--------------ATTRIBUTES-FROM-FUNCTIONS-OF-THIS-CLASS----------------------------------#

        #-----------calculation-attributes------------------------------------#
        self.inducedVel = self.getInducedVel()
        self.influenceCoeff = self.getInfluenceCoeff()
        self.RHS = self.getRHS()
        self.gamma = self.getCirculation()
        
        #----------post-processing-attributes---------------------------------#
        self.horseShoeVortexInducedVel = self.getHorseShoeVortexInducedVel()
        self.effectiveInducedVel = self.getEffectiveInducedVel()
        
        self.panelForceVec = self.getPanelForceVec() # panel force vectors in WING COOR SYSTEM
        self.panelMomentVec = self.getPanelMomentVec() # panel moment vectors in WING COOR SYSTEM
        self.transformedPanelForceVec = self.getTransformedPanelForceVec() # panel force vectors in AERODYNAMIC COOR SYSTEM
        self.transformedPanelMomentVec = self.getTransformedPanelMomentVec() # panel moment vectors in AERODYNAMIC COOR SYSTEM
        
        self.resultantPanelForceVec = self.getResultantPanelForceVec() # Resultant force vector in WING COOR SYSTEM
        self.resultantPanelMomentVec = self.getResultantPanelMomentVec() # Resultant Moment vector in WING COOR SYSTEM
        self.transformedResultantPanelForceVec = self.getTransformedResultantPanelForceVec() # Resultant force vector in AERODYNAMIC COOR SYSTEM
        self.transformedResultantPanelMomentVec = self.getTransformedResultantPanelMomentVec() # Resultant Moment vector in AERODYNAMIC COOR SYSTEM
        
 
    def getInducedVel(self): 
        inducedVel = np.zeros((self.geom.Nspan, self.geom.Nspan, 3))
                
        gamma = 1.
        for i in range(self.geom.Nspan): # i defines collocation PT
            for j in range(self.geom.Nspan): # j defines HorseShoe Vortex
                for k in range(5): # k defines number of segments of 1 HorseShoe Vortex
                    
                    r0 = self.geom.horseShoeVortexPts[j,k+1,:] - self.geom.horseShoeVortexPts[j,k,:] 
                    r1 = self.geom.evaluationPts[i,:] - self.geom.horseShoeVortexPts[j,k,:]
                    r2 = self.geom.evaluationPts[i,:] - self.geom.horseShoeVortexPts[j,k+1,:]

                    r0r1 = np.vdot(r0,r1)
                    r0r2 = np.vdot(r0,r2)
        
                    norm_r1 = np.linalg.norm(r1)
                    norm_r2 = np.linalg.norm(r2)
        
                    r1_cross_r2 = np.cross(r1,r2)
                    norm_r1_cross_r2 = np.linalg.norm(r1_cross_r2)
                                        
                    if ( abs(norm_r1) < 1e-13 ) or ( abs(norm_r2) < 1e-13 ) or ( abs(norm_r1_cross_r2) < 1e-13 ):
                        #print ("Yes")
                        Kfactor = 0
                    else:
                        Kfactor = ( gamma / (4.*np.pi*(norm_r1_cross_r2**2) ) )*( (r0r1/norm_r1) - (r0r2/norm_r2) )  # added for when divisoins are really big
                    
                    inducedVel[i,j,:] = inducedVel[i,j,:] + r1_cross_r2*Kfactor
                    
        return inducedVel
    
  
    def getInfluenceCoeff(self):
        influenceCoeff = np.zeros((self.geom.Nspan, self.geom.Nspan))
        
        for i in range(self.geom.Nspan): # i defines collocation PT
            for j in range(self.geom.Nspan): # j defines HorseShoe Vortex
                influenceCoeff[i,j] = np.vdot(self.inducedVel[i,j,:],self.geom.panelNormVec[i,:])
            
        return influenceCoeff

  
    def getRHS(self): 
        RHS = np.zeros(self.geom.Nspan)
       
        for i in range(self.geom.Nspan):
            RHS[i] = -np.vdot(self.QinfVec[:],self.geom.panelNormVec[i,:])
        
        return RHS

   
    def getCirculation(self):
        gamma = np.matmul(np.linalg.inv(self.influenceCoeff),self.RHS)
        
        return gamma
    
   
    def getHorseShoeVortexInducedVel(self): # this method is used ONCE in getEffIndVel()
        horseShoeVortexInducedVel = np.zeros((self.geom.Nspan, self.geom.Nspan, 3))
        
        for i in range(self.geom.Nspan): # i defines bound vortex center PT
            for j in range(self.geom.Nspan): # j defines HorseShoe Vortex
                for k in range(5): # k defines number of segments of 1 HorseShoe Vortex
        
                    r0 = self.geom.horseShoeVortexPts[j,k+1,:] - self.geom.horseShoeVortexPts[j,k,:] 
                    r1 = self.geom.boundVortexCenterPts[i,:] - self.geom.horseShoeVortexPts[j,k,:]
                    r2 = self.geom.boundVortexCenterPts[i,:] - self.geom.horseShoeVortexPts[j,k+1,:]

                    r0r1 = np.vdot(r0,r1)
                    r0r2 = np.vdot(r0,r2)
        
                    norm_r1 = np.linalg.norm(r1)
                    norm_r2 = np.linalg.norm(r2)
        
                    r1_cross_r2 = np.cross(r1,r2)
                    norm_r1_cross_r2 = np.linalg.norm(r1_cross_r2)
                                        
                    if ( abs(norm_r1) < 1e-13 ) or ( abs(norm_r2) < 1e-13 ) or ( abs(norm_r1_cross_r2) < 1e-13 ):
                        #print ("Yes")
                        Kfactor = 0
                    else:
                        Kfactor = ( self.gamma[j] / (4.*np.pi*(norm_r1_cross_r2**2) ) )*( (r0r1/norm_r1) - (r0r2/norm_r2) )  # added for when divisoins are really big
                    
                    horseShoeVortexInducedVel[i,j,:] = horseShoeVortexInducedVel[i,j,:] + r1_cross_r2*Kfactor
                    
        return horseShoeVortexInducedVel
    
 
    def getEffectiveInducedVel(self): 
        effectiveInducedVel = np.zeros((self.geom.Nspan,3))
        
        for i in range(self.geom.Nspan):
            effectiveInducedVel[i,:] = sum(self.horseShoeVortexInducedVel[i]) +  self.QinfVec[:] 
            
        return effectiveInducedVel
    
  
# #----------CALCULATION-OF-FORCES-AND-MOMENTS----------------------------------#    
    
    def getPanelForceVec(self): 
        panelForceVec = np.zeros((self.geom.Nspan,3))
        
        for i in range(self.geom.Nspan):
            rBound = self.geom.horseShoeVortexPts[i,3,:] - self.geom.horseShoeVortexPts[i,2,:]
            panelForceVec[i,:] = self.rho*self.gamma[i]*np.cross(self.effectiveInducedVel[i,:],rBound)
        
        return panelForceVec
    
    
    def getPanelMomentVec(self): 
        panelMomentVec = np.zeros((self.geom.Nspan,3))
        
        for i in range(self.geom.Nspan):
            panelMomentVec[i,:] = np.cross(self.geom.boundVortexCenterPts[i,:],self.panelForceVec[i,:])
        
        return panelMomentVec
    
    
    def getTransformedPanelForceVec(self):
        T = np.array([[math.cos(math.radians(self.geom.alpha)),0.,math.sin(math.radians(self.geom.alpha))],
                      [0.,1.,0.],
                      [-math.sin(math.radians(self.geom.alpha)), 0., math.cos(math.radians(self.geom.alpha))]])
        
        transformedPanelForceVec = np.zeros((self.geom.Nspan,3))
        for i in range(self.geom.Nspan):
            transformedPanelForceVec[i,:] = np.matmul(T,self.panelForceVec[i,:])
        
        return transformedPanelForceVec
    

    def getTransformedPanelMomentVec(self):
        T = np.array([[math.cos(math.radians(self.geom.alpha)),0.,math.sin(math.radians(self.geom.alpha))],
                      [0.,1.,0.],
                      [-math.sin(math.radians(self.geom.alpha)), 0., math.cos(math.radians(self.geom.alpha))]])
        
        transformedPanelMomentVec = np.zeros((self.geom.Nspan,3))
        for i in range(self.geom.Nspan):
            transformedPanelMomentVec[i,:] = np.matmul(T,self.panelMomentVec[i,:])
        
        return transformedPanelMomentVec


    def getResultantPanelForceVec(self):
        resultantPanelForceVec = sum(self.panelForceVec)
        
        return resultantPanelForceVec
    
    
    def getResultantPanelMomentVec(self):
        resultantPanelMomentVec = sum(self.panelMomentVec)
        
        return resultantPanelMomentVec
    
    
    def getTransformedResultantPanelForceVec(self):
        transformedResultantPanelForceVec = sum(self.transformedPanelForceVec)
        
        return transformedResultantPanelForceVec


    def getTransformedResultantPanelMomentVec(self):
        transformedResultantPanelMomentVec = sum(self.transformedPanelMomentVec)
        
        return transformedResultantPanelMomentVec

#----------CALCULATION-OF-AERODYNAMIC-COEFFICIENTS----------------------------#

    def getInducedDragCoeff(self):
        CDind = self.transformedResultantPanelForceVec[0]/(self.qInf*self.S_projected)
        return CDind
    
    def getSideForceCoeff(self):
        CY = self.transformedResultantPanelForceVec[1]/(self.qInf*self.S_projected)
        return CY
    
    def getLiftCoeff(self):
        CL = self.transformedResultantPanelForceVec[2]/(self.qInf*self.S_projected)
        return CL
    
    def getMomentCoeff_X(self):
        CM_X = self.transformedResultantPanelMomentVec[0]/(self.qInf*self.S_projected*self.cRef)
        return CM_X
    
    def getMomentCoeff_Y(self):
        CM_Y = self.transformedResultantPanelMomentVec[1]/(self.qInf*self.S_projected*self.cRef)
        return CM_Y
    
    def getMomentCoeff_Z(self):
        CM_Z = self.transformedResultantPanelMomentVec[2]/(self.qInf*self.S_projected*self.cRef)
        return CM_Z

#----------CALCULATION-OF-DISTRIBUTED-VALUES----------------------------#

    def getLocalLiftCoeffDistribution(self):
        cl_dist = np.zeros(self.geom.Nspan)
        
        for i in range(self.geom.Nspan):
            l = self.geom.horseShoeVortexPts[i,3,:] - self.geom.horseShoeVortexPts[i,2,:] 
            cl_dist[i] = ( self.transformedPanelForceVec[i,2]/np.linalg.norm(l) ) / (self.qInf*self.cRef)          
        
        return cl_dist
    
    def getLocalInducedDragCoeffDistribution(self):
        cdInd_dist = np.zeros(self.geom.Nspan)
        
        for i in range(self.geom.Nspan):
            l = self.geom.horseShoeVortexPts[i,3,:] - self.geom.horseShoeVortexPts[i,2,:] 
            cdInd_dist[i] = ( self.transformedPanelForceVec[i,0]/np.linalg.norm(l) ) / (self.qInf*self.cRef)    
        
        return cdInd_dist
    
    def getLocalMomentCoeff_YDistribution(self):
        cmY_dist = np.zeros(self.geom.Nspan)
        
        for i in range(self.geom.Nspan):
            l = self.geom.horseShoeVortexPts[i,3,:] - self.geom.horseShoeVortexPts[i,2,:] 
            cmY_dist[i] = ( self.transformedPanelMomentVec[i,1]/np.linalg.norm(l) ) / (self.qInf*self.cRef*self.cRef)    
        
        return cmY_dist
    
    def getLocalInducedAOAdistribution(self):
        indAOAdist = np.zeros(self.geom.Nspan)
        
        for i in range(self.geom.Nspan):
            indAOAdist[i] = -self.geom.alpha + math.degrees(math.atan(self.effectiveInducedVel[i,2]/self.effectiveInducedVel[i,0]))
    
        return indAOAdist