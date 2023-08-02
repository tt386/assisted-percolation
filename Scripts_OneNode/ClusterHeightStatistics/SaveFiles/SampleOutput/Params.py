import os
from os import path

import sys
import shutil
import numpy as np




##############################################################################
##############################################################################
##############################################################################


#Probability of a mutation occuring
mutprob = 0.

#Relative fitness of the mutant
fitness = 0.9

#Default system height
DomainHeight =500

"""
#Default system width
(can be changed a bit if the obstacles are periodicallty spaced)
"""
DomainWidth = 100

#Obstacle Radius
Radius = 60

"""
#If the obstacles are placed randomly or on a lattice,
 this is the area fraction they would have
"""
AreaFraction = 0.4

"""
#When height limit reached, this is the amount by which we generate new space
 above
"""
ExtensionHeight = 0

#The maximum possible height that the pests can explore
UpperDomainLimit = DomainHeight

"""
#There are a variety of end conditions found in
CoreFunctions/CoreSim.py RoughFront()
"""
EndCondition = 10000

SeededPos = (50,50)

Repeats = 100000


#Type of Obstacle distribution, with corresponding statistics

#Lattice Stuff
Lattice = False
LatticeType = 0
BottomLeftxdisp = 0
BottomLeftydisp = 5
LatticeSurfSep = 20

#Comment out if want to focus on surface sep
LatticeSurfSep = np.ceil(
    Radius * np.sqrt(2*np.pi/(np.sqrt(3) * AreaFraction)) *
    (1. - 2 * np.sqrt(np.sqrt(3) * AreaFraction/(2.* np.pi))))

#Random
Random = False

#Single
Single = False

SaveFileDir = ("SaveFiles/"
    "Height_%d_Width_%d_Fitness_%0.3f_Repeats_%d_SeededHeight_%d"%(
    DomainHeight,DomainWidth,fitness,Repeats,SeededPos[1]))



if Lattice or Random or Single:
    print("INCORRECT:Y APPLIED Obstacle Distribution")
    sys.exit()


if not os.path.isdir(SaveFileDir):
    os.mkdir(SaveFileDir)

shutil.copyfile("Params.py", SaveFileDir+"/Params.py")
#shutil.copyfile("Script.py", SaveFileDir+"/Script.py")
