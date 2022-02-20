import numpy as np

from src.wingGeometry import wingGeometry
from src.solution import solution
from src.postProcess import postProcess

#---------------------INPUT-PARAMETERS------------------------------------------------------#

b= 20. # wing span (in meters)
wakeLenFactor = 40. # denotes how many factors of span b does the wake extend  
Lam = 15. # sweep angle (in degrees) -- range -90 to 90
beta = 0. # dihedral angle (in degrees) -- range -90 to 90
cRoot = 2. # wing root chord length
cTip = 1. # wing tip chord length
Nspan = int(80) # TOTAL numbers of panels, must be EVEN
rho = 1. # density in kg/m3
freestreamVel = 1. # freestream velocity in m/s
alphaRange = np.linspace(0,15,16) # range of geometric AOA considered
cRef = 2. # reference chord length (should generally be equal to cRoot)
HSbreak = 0.25 # location where HorseShoe vorticies break, for example 0.25 denotes the the HorseShoe vorticies break at 0.25 the local chord length behind the trailing edge
typeSpacing = "uniform" # "uniform" or "cosine"
typeEvalPt = "center" # "center" or "glauert" -----> NOTE: if UNIFORM spacing us used, then typeEvalPt must be "center"

#-------------------------------------------------------------------------------------------#
AOA, liftCoeff, inducedDragCoeff, momentCoeff_Y, gammaDistribution, liftCoeffDistribution, inducedAOAdistribution, inducedDragCoeffDistribution, momentCoeff_Ydistribution = [[] for j in range(9)]

print('#--------------------GEOMETRICAL-PARAMETERS--------------------#')
print('Wing Span =', b, '(m)')
print('Wing Root Chord =', cRoot, '(m)')
print('Wing Tip Chord =', cTip, '(m)')
print('Wing Sweep Angle =', Lam, '(in degrees)')
print('Wing Dihedral Angle =', beta, '(in degrees)')

for i in range(len(alphaRange)): # angle of attack in degrees

    #----------------------------SOLUTION------------------------------------# 
    alpha = alphaRange[i]

    #--------------------WING-GEOMETRY-CALCULATION---------------------------# 
    wing = wingGeometry(b, Lam, beta, cRoot, cTip, Nspan, alpha, wakeLenFactor, HSbreak, typeSpacing, typeEvalPt)

    #----------------------AERO-VALUES-CALCULATION---------------------------# 
    sol = solution(wing, freestreamVel, rho, cRef)

    AR = b**2/sol.S_projected
    if i==0:
        print('Wing Aspect Ratio =', np.round(AR, 4))
        print()
        print('#---------------------------SOLUTION---------------------------#')
        print ("AOA =", np.round(alpha, 4), "--> CL =", np.round(sol.getLiftCoeff(), 5 ),
               ", CDind =", np.round(sol.getInducedDragCoeff(),5), ", CMy =", np.round(sol.getMomentCoeff_Y(),5) )  
    else:
        print ("AOA =", np.round(alpha, 4), "--> CL =", np.round(sol.getLiftCoeff(), 5 ),
               ", CDind =", np.round(sol.getInducedDragCoeff(),5), ", CMy =", np.round(sol.getMomentCoeff_Y(),5) )  

    #----------------LISTS-OF-INTEGRAL-AERO-VALUES------------------------------#    
    AOA.append( alpha )
    liftCoeff.append( sol.getLiftCoeff() )
    inducedDragCoeff.append( sol.getInducedDragCoeff() )
    momentCoeff_Y.append( sol.getMomentCoeff_Y() )

    #----------------LIST-OF-DISTRIBUTED-AERO-VALUES----------------------------#
    gammaDistribution.append( sol.gamma )
    liftCoeffDistribution.append( sol.getLocalLiftCoeffDistribution() )
    inducedAOAdistribution.append( sol.getLocalInducedAOAdistribution() )
    inducedDragCoeffDistribution.append( sol.getLocalInducedDragCoeffDistribution() )
    momentCoeff_Ydistribution.append( sol.getLocalMomentCoeff_YDistribution() )

#-------------------------POST-PROCESSING--------------------------------------#

PP = postProcess(wing, sol, AOA, liftCoeff, inducedDragCoeff, momentCoeff_Y, gammaDistribution, liftCoeffDistribution, inducedAOAdistribution, inducedDragCoeffDistribution,momentCoeff_Ydistribution)
PP.getIntegralAeroCoeff_textfile() # creates textfile of cl, cd, cm polars
PP.getDistributedAeroCoeff_textfile() # creates textfile of distributed cl, cd, induced aoa, cm and gamma along wingspan for a specific geometric AOA
PP.getPlotGeometry() # plots wing geometry for a specific geometric angle of attack
PP.getPlotAeroCoeffPolars() # plots cl, cd, cm polars 
PP.getPlotAeroCoeffDistribution() # plots distributed cl, cd, induced aoa, cm and gamma along wingspan for a specific geometric AOA