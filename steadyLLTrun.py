# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:23:25 2020

@author: n.tatomir
"""
import numpy as np
import math 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from src.wingGeometry import wingGeometry
from src.solution import solution
from src.postProcess import postProcess

#------------------------------------------------------------------------------------------#

b= 10. # span (in meters)
wakeLenFactor = 40. # denotes how many factors of span b does the wake extend  
Lam = 20. # sweep angle (in degrees) -- range -90 to 90
beta = 0. # dihedral angle (in degrees) -- range -90 to 90
cRoot = 1. # wing root chord length
cTip = 1. # wing tip chord
Nspan = int(20) # TOTAL numbers of panels, must be EVEN
rho = 1. # density in kg/m3
freestreamVel = 1. # freestream angle of attack
alphaRange = np.linspace(5,10,5) # range of geometric AOA considered
cRef = 1. # reference chord length
HSbreak = 0.25 # location where HorseShoe vorticies break, for example 0.25 denotes the the HorseShoe vorticies break at 0.25 the local chord length behind the trailing edge
typeSpacing = "uniform" # "uniform" or "cosine"
typeEvalPt = "center" # "center" or "glauert" -----> NOTE: if UNIFORM spacing us used, then typeEvalPt must be "center"

#-------------------------------------------------------------------------------------------#
AOA, liftCoeff, inducedDragCoeff, momentCoeff_Y, gammaDistribution, liftCoeffDistribution, inducedAOAdistribution, inducedDragCoeffDistribution, momentCoeff_Ydistribution = [[] for j in range(9)]

for i in range(len(alphaRange)): # angle of attack in degrees
    alpha = alphaRange[i]
    
    wing = wingGeometry(b, Lam, beta, cRoot, cTip, Nspan, alpha, wakeLenFactor, HSbreak, typeSpacing, typeEvalPt)
    
    sol = solution(wing, freestreamVel, rho, cRef)

    AR = b**2/sol.S_projected
    print ("AR =", AR),    print ("AOA =", alpha)
    print ("CL =", sol.getLiftCoeff())
    print ("CDind =", sol.getInducedDragCoeff())
    print ("CMy =", sol.getMomentCoeff_Y())
    
    #-------------------Lists-of-Integral-Values------------------------------#    
    AOA.append( alpha )
    liftCoeff.append( sol.getLiftCoeff() )
    inducedDragCoeff.append( sol.getInducedDragCoeff() )
    momentCoeff_Y.append( sol.getMomentCoeff_Y() )
    #-------------------List-of-Distributed-Values----------------------------#
    gammaDistribution.append( sol.gamma )
    liftCoeffDistribution.append( sol.getLocalLiftCoeffDistribution() )
    inducedAOAdistribution.append( sol.getLocalInducedAOAdistribution() )
    inducedDragCoeffDistribution.append( sol.getLocalInducedDragCoeffDistribution() )
    momentCoeff_Ydistribution.append( sol.getLocalMomentCoeff_YDistribution() )

#-------------------------Post-Processing--------------------------------------#

PP = postProcess(wing, sol, AOA, liftCoeff, inducedDragCoeff, momentCoeff_Y, gammaDistribution, liftCoeffDistribution, inducedAOAdistribution, inducedDragCoeffDistribution,momentCoeff_Ydistribution)
PP.getIntegralAeroCoeff_textfile() # creates textfile of cl, cd, cm polars
PP.getDistributedAeroCoeff_textfile() # creates textfile of distributed cl, cd, induced aoa, cm and gamma along wingspan for a specific geometric AOA
PP.getPlotGeometry() # plots wing geometry for a specific geometric angle of attack
PP.getPlotAeroCoeffPolars() # plots cl, cd, cm polars 
PP.getPlotAeroCoeffDistribution() # plots distributed cl, cd, induced aoa, cm and gamma along wingspan for a specific geometric AOA