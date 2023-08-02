#import Params
from Params import *

import sys
sys.path.insert(0,'../../CoreFunctions')

import CoreSim

#import ObstaclesFuns

import time

import os
from os import path

import copy

import numpy as np

import subprocess

from resource import getrusage, RUSAGE_SELF

import random
#############################################################################
###Find data size############################################################
#############################################################################
def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj,
             (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size













##############################################################################
##Initialise##################################################################
##############################################################################
initstarttime = time.time()

ObsCenterList,ObsCenterDict = [],{}

if Random:
    ObsCenterList,ObsCenterDict = CoreSim.Obstacles(
        Radius,AreaFraction,DomainHeight,DomainWidth)

elif Lattice:
    ObsCenterList,ObsCenterDict,DomainWidth = CoreSim.LatticeObstacles(
        Radius,DomainHeight,DomainWidth,LatticeType,LatticeSurfSep,
        BottomLeftxdisp,BottomLeftydisp)

elif Single:
    ObsCenterList,ObsCenterDict = CoreSim.SingleObstacle(
        Radius,DomainHeight,DomainWidth,int(DomainWidth/2),int(DomainHeight/2))

else:
    ObsCenterList,ObsCenterDict = CoreSim.Obstacles(
        Radius,0,DomainHeight,DomainWidth) 

PopulationDict,UnevolvedList,EvolvedList = CoreSim.Initialise(
    0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)

print("Obstacle Center List: " + str(ObsCenterList))

initendtime = time.time()

print(str(initendtime-initstarttime))



#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Define the (optional) random seed')
parser.add_argument('-r','--randomseed',type=int,required=False,help='The random seed')
args = parser.parse_args()

if args.randomseed is None:
    print("No seed provided")
    randomseed = random.randrange(sys.maxsize)

else:
    randomseed = int(args.randomseed)
    print("seed provided:",randomseed)

random.seed(randomseed)

##############################################################################
###Sim########################################################################
##############################################################################

DataSizeList = []

OutputData = []
SimStep = 0
ExpTime = 0.

#Data for plotting: this will be for the first frame
for key in PopulationDict:
    P = PopulationDict[key]
    newoutputrow = [0,0,0,0]
    newoutputrow[0] = 0                     #time
    newoutputrow[1] = P.GetPos()[0]         #x
    newoutputrow[2] = P.GetPos()[1]         #y
    newoutputrow[3] = P.GetMut()            #Species/colour
    OutputData.append(newoutputrow)




simstarttime = time.time()


#Below is the main loop executed amongst al the scripts. The system will loop
#util a predefined end condition is met, but by including it in every script
#there is the option to define specialist end conditions, such as stopping
#when there are no more mutants etc. We could also define things such as 
#turning the mutation rate on and off given certain conditions.

EndConditionMet = False
while EndConditionMet == False:

    #The 'time' elapsed in a single time step
    ExpTime += 1.0/(len(UnevolvedList) + fitness*len(EvolvedList))

    #For understanding the data usage
    DataSizeList.append((SimStep,getrusage(RUSAGE_SELF).ru_maxrss)) 

    #Some sample outputs to get a feel for what's going on
    if (SimStep % 10000) == 0:
        print("SimStep = " + str(SimStep))
        print("WT:Mutant ",len(UnevolvedList),len(EvolvedList))
        print("PopDictSize: " + str(DataSizeList[-1][1]))
    (PopulationDict, 
    EvolvedList, 
    UnevolvedList, 
    EndConditionMet, 
    ChildPest) = CoreSim.RoughFront(
        DomainHeight,DomainWidth,mutprob,fitness,PopulationDict,UnevolvedList,
        EvolvedList,Radius,ObsCenterList,ObsCenterDict,EndCondition)

    SimStep += 1

    #Add output to the list
    newoutputrow = [0,0,0,0]
    newoutputrow[0] = ExpTime                     #time
    newoutputrow[1] = ChildPest.GetPos()[0]      #x
    newoutputrow[2] = ChildPest.GetPos()[1]      #y
    newoutputrow[3] = ChildPest.GetMut()                     #color
    OutputData.append(newoutputrow)


    
    #Extend System if Necessary
    if (DomainHeight - ChildPest.GetPos()[1]) < (np.sqrt(3)/2)*Radius:
        if Random:
            ObsCenterList,ObsCenterDict,DomainHeight = CoreSim.ExtendSystem(
                Radius,AreaFraction,DomainHeight,DomainWidth,ObsCenterList,
                ObsCenterDict,ExtensionHeight)
        else:
            ObsCenterList,ObsCenterDict,DomainHeight = CoreSim.ExtendSystem(
                Radius,0,DomainHeight,DomainWidth,ObsCenterList,
                ObsCenterDict,ExtensionHeight)
        print("SYSTEM EXTENDED TO HEIGHT: " + str(DomainHeight))
    
    #Break if there's a problem
    """
    if len(UnevolvedList) == 0:
        EndConditionMet = True
    """
    if ChildPest.GetPos()[1] >= UpperDomainLimit-1:
        EndConditionMet = True

simendtime = time.time()

print("Time taken for sim: " + str(simendtime-simstarttime))





##############################################################################
##EndState Analysis###########################################################
##############################################################################
#Measure the area ratio of obstacles, so I can test how accurate my procedure
#of generating the obstacles is. Useful for random obstacle distribution as it
# will enable me to create an average.


MeasuredAreaRatio = 0
if Random: 
    MeasuredAreaRatio = CoreSim.MeasureAreaFraction(
        DomainWidth,DomainHeight,Radius,ObsCenterList,ObsCenterDict)

elif Lattice:
    MeasuredAreaRatio = CoreSim.MeasureAreaFractionLattice(
        DomainWidth,DomainHeight,Radius,ObsCenterList,ObsCenterDict,
        BottomLeftydisp)

print("Target Area Ratio: %0.3f, Measured Area Ratio: %0.3f"%(AreaFraction,
        MeasuredAreaRatio))

#####################
##Save Data##########
#####################
OutputDatafilename = SaveFileDir + '/timecourse.npz'
np.savez(OutputDatafilename,
    OutputData=OutputData,
    DomainHeight=DomainHeight,
    DomainWidth=DomainWidth,
    fitness=fitness,
    mutprob=mutprob,
    Radius=Radius,
    ObsCenterList=ObsCenterList,
    AreaFraction=AreaFraction,
    MeasuredAreaRatio=MeasuredAreaRatio,
    LatticeSurfSep=LatticeSurfSep,
    DataSizeList=DataSizeList,
    randomseed=randomseed)

with open(SaveFileDir+'/randomseed.txt', 'w') as f:
    f.write('seed:'+str(randomseed))


print("Final number of Evolved, Unevolved:",len(EvolvedList),
    len(UnevolvedList))

