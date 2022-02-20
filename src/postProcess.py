import numpy as np
import math 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

class postProcess():
    
    def __init__(self, geom, calc,  AOA, liftCoeff, inducedDragCoeff, momentCoeff_Y, gammaDistribution, liftCoeffDistribution, inducedAOAdistribution, inducedDragCoeffDistribution, momentCoeff_Ydistribution):
        self.geom = geom
        self.calc = calc
        self.AOA = AOA
        self.liftCoeff = liftCoeff
        self.inducedDragCoeff = inducedDragCoeff
        self.momentCoeff_Y = momentCoeff_Y
        self.gammaDistribution = gammaDistribution
        self.liftCoeffDistribution = liftCoeffDistribution
        self.inducedAOAdistribution = inducedAOAdistribution
        self.inducedDragCoeffDistribution = inducedDragCoeffDistribution
        self.momentCoeff_Ydistribution = momentCoeff_Ydistribution
        
    def getIntegralAeroCoeff_textfile(self):
        
        if len(self.AOA) > 1:
            with open("aeroIntegralValuesPolar.txt", "w") as file:
                file.write("\nWING CHARACTERISTICS:\n\n")
                
                file.write("span --> b = {}\n".format(self.geom.b))
                file.write("sweep angle --> Lambda = {}\n".format(self.geom.Lam))
                file.write("dihedral angle --> beta = {}\n\n".format(self.geom.beta))
                
                file.write("root chord --> cRoot = {}\n".format(self.geom.cRoot))
                file.write("tip chord --> cTip = {}\n".format(self.geom.cTip))
                file.write("taper ratio --> cRoot/cTip = {}\n\n".format(format(self.geom.cRoot/self.geom.cTip,'.4f')))
                
                file.write("Actual Wing Surface Area --> S_actual = {}\nProjected Wing Surface Area (onto xy plane of Wing Axis) --> S_proj = {}\n".format(format(self.calc.S_actual,'.4f'), format(self.calc.S_projected,'.4f')))
                file.write("Aspect Ratio --> AR (b^2/S_proj) = {}\n\n".format(format(self.geom.b**2/self.calc.S_projected,'.4f')))
                
                file.write("reference chord (= cRoot) --> c_ref = {}\n".format(self.geom.cRoot))
                file.write("Spacing Type --> {}\n".format(self.geom.typeSpacing))
                file.write("Evaluation Point Type --> {}\n".format(self.geom.typeEvalPt))
                file.write("HorseShoe Break Point Location --> {}\n".format(self.geom.HSbreak))
                file.write("number of Panels --> Nspan = {}\n\n\n\n".format(self.geom.Nspan))
                
                file.write("     Alpha (in degrees)                CL (-)             CDind (-)          CMy (-)\n")
                file.write("------------------------------------------------------------------------------------------\n")
                
                for i in range(len(self.AOA)):
                    file.write("          {}                      {}           {}           {}\n".format(format(self.AOA[i],'.3f'), format(self.liftCoeff[i],'.8f'), format(self.inducedDragCoeff[i],'.8f'), format(self.momentCoeff_Y[i],'.8f') ))    
        else:
            pass
    
     
    def getDistributedAeroCoeff_textfile(self):
        
        if len(self.AOA) == 1:
            with open("aeroDistributedValues.txt", "w") as file:
                file.write("\n\nThis file gives the distributed values along the wingspan for a single geometric angle of attack!")
                file.write("\nWING CHARACTERISTICS:\n\n")
                
                file.write("Angle of Attack (in degrees): {}\n".format(self.AOA[0]))
                file.write("span --> b = {}\n".format(self.geom.b))
                file.write("sweep angle --> Lambda = {}\n".format(self.geom.Lam))
                file.write("dihedral angle --> beta = {}\n\n".format(self.geom.beta))
            
                file.write("root chord --> cRoot = {}\n".format(self.geom.cRoot))
                file.write("tip chord --> cTip = {}\n".format(self.geom.cTip))
                file.write("taper ratio --> cRoot/cTip = {}\n\n".format(format(self.geom.cRoot/self.geom.cTip,'.4f')))
                
                file.write("Actual Wing Surface Area --> S_actual = {}\nProjected Wing Surface Area (onto xy plane of Wing Axis) --> S_proj = {}\n".format(format(self.calc.S_actual,'.4f'), format(self.calc.S_projected,'.4f')))
                file.write("Aspect Ratio --> AR (b^2/S_proj) = {}\n\n".format(format(self.geom.b**2/self.calc.S_projected,'.4f')))
                
                file.write("reference chord (= cRoot) --> c_ref = {}\n".format(self.geom.cRoot))
                file.write("Spacing Type --> {}\n".format(self.geom.typeSpacing))
                file.write("Evaluation Point Type --> {}\n".format(self.geom.typeEvalPt))
                file.write("HorseShoe Break Point Location --> {}\n".format(self.geom.HSbreak))
                file.write("number of Panels --> Nspan = {}\n\n\n\n".format(self.geom.Nspan))
                    
                file.write("-------------------INTEGRAL VALUES OF AERODYNAMIC COEFFICIENTS----------------------------\n")
                file.write("     Alpha (in degrees)                CL (-)             CDind (-)            CMy (-)\n")
                file.write("------------------------------------------------------------------------------------------\n")
                file.write("          {}                      {}           {}          {}\n\n\n\n".format(format(self.AOA[0],'.3f'), format(self.liftCoeff[0],'.8f'), format(self.inducedDragCoeff[0],'.8f'), format(self.momentCoeff_Y[0],'.8f') ))    
                
                file.write("---------------------------------------DISTRIBUTED VALUES OF AERODYNAMIC COEFFICIENTS---------------------------------------------------\n")           
                file.write("     span y (m)                Gamma (m^2/s)                indAOA (in degrees)             cl (-)          cd_ind (-)          cm_y (-)\n")
                file.write("----------------------------------------------------------------------------------------------------------------------------------------\n")
                
                for i in range(self.geom.Nspan):
                    file.write("    {}                   {}                   {}               {}         {}        {}\n".format(format(self.geom.boundVortexCenterPts[i,1],'.8f'), format(self.gammaDistribution[0][i],'.8f'), format(self.inducedAOAdistribution[0][i],'.8f'), format(self.liftCoeffDistribution[0][i],'.8f'), format(self.inducedDragCoeffDistribution[0][i],'.8f'), format(self.momentCoeff_Ydistribution[0][i],'.8f') ))    
        else:
            pass

    def getPlotGeometry(self):
        
        if len(self.AOA) == 1:
            subplotNum = 221
            subtitle = ['Front View', 'Top View', 'Side View #1', 'Side View #2']
            orientation_1 = [0,90,30,0]
            orientation_2 = [180,180,225,225+45]
            fig = plt.figure(figsize=(32.0,16.0))
            for i in range(4):
                ax = fig.add_subplot(subplotNum, projection = '3d')
            
                #----------PLOT-WING------------------------------------------------------#
                ax.plot(self.geom.panelPts[0,:,0], self.geom.panelPts[0,:,1], self.geom.panelPts[0,:,2], 'b--', linewidth = 0.5)
                ax.plot(self.geom.panelPts[1,:,0], self.geom.panelPts[1,:,1], self.geom.panelPts[1,:,2], 'b--', linewidth = 0.5)
                ax.plot(self.geom.panelPts[:,0,0], self.geom.panelPts[:,0,1], self.geom.panelPts[:,0,2], 'b--', linewidth = 0.5)
                ax.plot(self.geom.panelPts[:,-1,0], self.geom.panelPts[:,-1,1], self.geom.panelPts[:,-1,2], 'b--', linewidth = 0.5)
                
                #---------PLOT-PANEL------------------------------------------------------#
                for j in range(self.geom.Nspan+1):
                    ax.plot([self.geom.panelPts[0,j,0],self.geom.panelPts[1,j,0]], [self.geom.panelPts[0,j,1],self.geom.panelPts[1,j,1]] , [self.geom.panelPts[0,j,2],self.geom.panelPts[1,j,2]] , 'm-', linewidth = 0.25)
                
                #---------PLOT-BOUND-VORTEX-CENTER-POINTS---------------------------------#
                ax.plot(self.geom.boundVortexCenterPts[:,0], self.geom.boundVortexCenterPts[:,1], self.geom.boundVortexCenterPts[:,2], 'k.', markersize = 2)
                
                #---------PLOT-EVALUATION-POINTS------------------------------------------#
                ax.plot(self.geom.evaluationPts[:,0], self.geom.evaluationPts[:,1], self.geom.evaluationPts[:,2], 'r.', markersize = 2)
                
                #---------PLOT-HORSESHOE-VORTEX-------------------------------------------#
                for j in range(self.geom.Nspan):
                    if j == 0:
                        ax.plot(self.geom.horseShoeVortexPts[j,:,0],self.geom.horseShoeVortexPts[j,:,1],self.geom.horseShoeVortexPts[j,:,2],'g',linewidth = 0.5, label = 'HorseShoe Vortices')
                    else:
                        ax.plot(self.geom.horseShoeVortexPts[j,:,0],self.geom.horseShoeVortexPts[j,:,1],self.geom.horseShoeVortexPts[j,:,2],'g',linewidth = 0.5)
                    ax.plot(self.geom.horseShoeVortexPts[j,:,0],self.geom.horseShoeVortexPts[j,:,1],self.geom.horseShoeVortexPts[j,:,2],'g.',markersize = 2)
                
            #    for j in range(3):
            #        ax.plot(wing.grouped_HSpts[j,:,0],wing.grouped_HSpts[j,:,1],wing.grouped_HSpts[j,:,2], 'g.', markersize = 1)
                #-------PLOT-ORIGIN-------------------------------------------------------#
                ax.quiver(0,0,0,0.75,0,0, color = 'orange')
                ax.quiver(0,0,0,0,0.75,0, color = 'orange')
                ax.quiver(0,0,0,0,0,0.75, color = 'orange', label = "Origin of Wing XYZ Coordinate System")
                ax.quiver(0,0,0,0.75/math.cos(math.radians(self.AOA[0])),0,0.75*math.sin(math.radians(self.AOA[0]))/math.cos(math.radians(self.AOA[0])), color = 'c')
                ax.quiver(0,0,0,0,0.75,0, color = 'c')
                ax.quiver(0,0,0,-0.75*math.sin(math.radians(self.AOA[0]))/math.cos(math.radians(self.AOA[0])),0,0.75/math.cos(math.radians(self.AOA[0])), color = 'c', label = "Origin of Aerodynamic XYZ Coordinate System")
                #------SET-AXIS-LIMITS----------------------------------------------------#
                ax.set_xlim(-np.ceil(np.max(abs(self.geom.panelPts))),np.ceil(np.max(abs(self.geom.panelPts))))   
                ax.set_ylim(-np.ceil(np.max(abs(self.geom.panelPts)))-1,np.ceil(np.max(abs(self.geom.panelPts)))+1)
                ax.set_zlim(-np.ceil(np.max(abs(self.geom.panelPts))),np.ceil(np.max(abs(self.geom.panelPts))))
                ax.set_xlabel("X (m)", fontsize = 14)
                ax.set_ylabel("Y (m)", fontsize = 14)  
                ax.set_zlabel("Z (m)", fontsize = 14)
                ax.legend(loc = 'upper right', bbox_to_anchor=(0.3, 1.1))
                #--------MISCELENOUS------------------------------------------------------#
                ax.view_init(orientation_1[i], orientation_2[i])
                plt.title('{}'.format(subtitle[i]), fontsize = 18)
                subplotNum += 1
            plt.suptitle("Wing Geometry From Different Points of View",fontsize = 24)
            plt.savefig("wingGeometry.png", bbox_inches='tight', dpi = 250)
        else:
            pass
    
    
    def gridlines(self):
        plt.minorticks_on()
        plt.grid(zorder = 0, which='major', axis='x', linewidth=0.75, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='major', axis='y', linewidth=0.75, linestyle='-', color='0.75')
        plt.grid(zorder = 0, which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.75')
    
    def getPlotAeroCoeffPolars(self):
        if len(self.AOA) > 1:
            plt.figure(2,figsize=(24.0,8.0))        
            
            plt.subplot(131)
            plt.plot(self.AOA, self.liftCoeff, 'b-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Angle of Attack $\\alpha$ $(in$ $degrees)$', fontsize = 12)
            plt.ylabel('Lift Coefficient $C_L$ $(-)$', fontsize = 12)
            plt.title('$C_L$ vs. $\\alpha$', fontsize = 16)
    
            plt.subplot(132)
            plt.plot(self.AOA, self.inducedDragCoeff, 'r-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Angle of Attack $\\alpha$ $(in$ $degrees)$', fontsize = 12)
            plt.ylabel('Induced Drag Coefficient $C_{D,ind}$ $(-)$', fontsize = 12)
            plt.title('$C_{D,ind}$ vs. $\\alpha$', fontsize = 16)
    
            plt.subplot(133)
            plt.plot(self.AOA, self.momentCoeff_Y, 'g-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Angle of Attack $\\alpha$ $(in$ $degrees)$', fontsize = 12)
            plt.ylabel('Moment Coefficient about Aerodynamic Y axis $C_{M,y}$ $(-)$', fontsize = 12)
            plt.title('$C_{M,y}$ vs. $\\alpha$', fontsize = 16)
            
            plt.savefig("aeroCoeffPolars.png", bbox_inches='tight', dpi = 250)
        else:
            pass
        
    def getPlotAeroCoeffDistribution(self):
        if len(self.AOA) == 1:
            plt.figure(3,figsize=(32.0,12.0))        
            
            plt.subplot(231)
            plt.plot(self.geom.boundVortexCenterPts[:,1], self.gammaDistribution[0][:], 'k-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Wing Span $y$ $(m)$', fontsize = 12)
            plt.ylabel('Local Circulation $\\Gamma$ $(m^2/s)$', fontsize = 12)
            plt.title('Local $\\Gamma$ Distrubtion Along Wing Span', fontsize = 16)
            
            plt.subplot(232)
            plt.plot(self.geom.boundVortexCenterPts[:,1], self.inducedAOAdistribution[0][:], 'r--', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Wing Span $y$ $(m)$', fontsize = 12)
            plt.ylabel('Local Induced AOA $\\alpha_{ind}$ $(in$ $degrees)$', fontsize = 12)
            plt.title('Local $\\alpha_{ind}$ Distribtuion Along Wing Span', fontsize = 16)
            
            plt.subplot(233)
            plt.plot(self.geom.boundVortexCenterPts[:,1], self.inducedDragCoeffDistribution[0][:], 'r-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Wing Span $y$ $(m)$', fontsize = 12)
            plt.ylabel('Local Induced Drag Coefficient $c_{D,ind}$ $(-)$', fontsize = 12)
            plt.title('Local $c_{D,ind}$ Distribution Along Wing Span', fontsize = 16)
            
            plt.subplot(234)
            plt.plot(self.geom.boundVortexCenterPts[:,1], self.liftCoeffDistribution[0][:], 'b-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Wing Span $y$ $(m)$', fontsize = 12)
            plt.ylabel('Local Lift Coefficient $c_l$ $(-)$', fontsize = 12)
            plt.title('Local $c_l$ Distribution Along Wing Span', fontsize = 16)
            
            plt.subplot(235)
            plt.plot(self.geom.boundVortexCenterPts[:,1], self.momentCoeff_Ydistribution[0][:], 'g-', linewidth = 0.75)
            self.gridlines()
            plt.xlabel('Wing Span $y$ $(m)$', fontsize = 12)
            plt.ylabel('Local Moment Coefficient about Aerodynamic Y axis $c_{m,Y}$ $(-)$', fontsize = 12)
            plt.title('Local $c_{m,y}$ Distribution Along Wing Span', fontsize = 16)
            
           
            plt.savefig("aeroCoeffDistribution_AOA{}.png".format(self.AOA[0]), bbox_inches='tight', dpi = 250)
        else:
            pass
        
        
        