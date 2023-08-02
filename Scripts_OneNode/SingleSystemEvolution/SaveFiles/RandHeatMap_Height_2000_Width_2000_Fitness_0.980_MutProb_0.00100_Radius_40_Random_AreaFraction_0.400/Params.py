import os
from os import path

import sys

import numpy as np

import shutil,os
##############################################################################
##############################################################################
##############################################################################

#Probability of a mutation occuring
mutprob = 1e-3#1e-3#0.05

#Relative fitness of the mutant
fitness = 0.98#0.99#0.95

#Default system height
DomainHeight =2000#1000#1000

"""
#Default system width 
(can be changed a bit if the obstacles are periodicallty spaced)
"""
DomainWidth = 2000#1000#1000

#Obstacle Radius
Radius = 40#5

"""
#If the obstacles are placed randomly or on a lattice,
 this is the area fraction they would have
"""
AreaFraction = 0.4

"""
#When height limit reached, this is the amount by which we generate new space
 above
"""
ExtensionHeight = 20 

#The maximum possible height that the pests can explore
UpperDomainLimit = DomainHeight#150 

"""
#There are a variety of end conditions found in 
CoreFunctions/CoreSim.py RoughFront()
"""
EndCondition = 0#500



##############################################################################
##############################################################################
##############################################################################


#Type of Obstacle distribution, with corresponding statistics

#Lattice Stuff
Lattice = False
LatticeType = 0
BottomLeftxdisp = 0
BottomLeftydisp = 10
LatticeSurfSep = 12
"""
if Lattice:
    AreaFraction = (2*np.pi*Radius**2)/(np.sqrt(3)*(2*Radius+LatticeSurfSep)**2)
"""
#Comment out if want to focus on surface sep
LatticeSurfSep = np.ceil(
    Radius * np.sqrt(2*np.pi/(np.sqrt(3) * AreaFraction)) * 
    (1. - 2 * np.sqrt(np.sqrt(3) * AreaFraction/(2.* np.pi))))
print("LatticeSurfSep",LatticeSurfSep)


#Random
Random = True

#Single
Single = False

SaveFileDir = "SaveFiles/Height_%d_Width_%d_Fitness_%0.3f_MutProb_%0.5f"%(
    DomainHeight,DomainWidth,fitness,mutprob)

if (not Lattice) and (not Single) and Random:
    SaveFileDir += "_Radius_%d_Random_AreaFraction_%0.3f"%(
        Radius,AreaFraction)

elif (Lattice) and (not Single ) and (not Random): 
    SaveFileDir += "_Radius_%d_Lattice_LatticeSurfSep_%0.3f_AreaFraction_%0.3f"%(
        Radius,LatticeSurfSep,AreaFraction)

elif (not Lattice) and Single and (not Random):
    SaveFileDir += "_Radius_%d_SingleObstacle"%(Radius)

elif (not Lattice) and (not Single ) and (not Random):
    SaveFileDir += "_NoObstacles"

else:
    print("Misunderstood Obstacle Distribution")
    sys.exit()
os.mkdir(SaveFileDir)

shutil.copy("Params.py",SaveFileDir)
shutil.copy("Script.py",SaveFileDir)
