
from Params import *

import sys
sys.path.insert(0,'../../CoreFunctions')

import CoreSim

import time

import random

import os
from os import path

import copy

import numpy as np

from argparse import ArgumentParser


#############################################################################
#############################################################################
#############################################################################



def InfectionProbability(
        DomainHeight,DomainWidth,Radius,mutprob,fitness,EndCondition,
        PopulationDict,UnevolvedList,EvolvedList,ObsCenterList,ObsCenterDict,
        R):
    """
    Run the rough front eden model in the presence of a single patch, and 
    record if the patch becomes invaded by M.
    
    ARGS:
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    Radius: int
        The radius of the patches.
    mutprob: float
        The probability that a given occupied site is a M.
    fitness: float
        The fitness of M: the ratio to which a M is chosen relative to a WT.
    EndCondition: int
        This is the option corresponding to the condition for which we will
        cease executing the rough front Eden Model.
    PopulationDict: dict
        A dictionary of pests: the keys are their locations within the system,
        and the values are the Pest class objects.
    UnevolvedList: list of tuples
        The coordinates of WT pests on the front.
    EvolvedList: list of tuples
        The coordinated of M pests on the front.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    R: int
        The enumerated repeat value of this repeat.

    RETURNS
    list of tuples:
        List of coordinates of the patches in the system.
    int:
        A signal as to why the simulation finished, so we can check if this
        should be included as a repeat.
    """


    ENDREASON = 0   #The reason why the system ends

    #Main Sim
    EndConditionMet = False
    SimStep = 0
    while EndConditionMet == False:
        ##########################################
        ###CORE###################################
        ##########################################
        if (SimStep % 10000) == 0:
            print("SimStep = " + str(SimStep))
            print("#Wt, #Mut: %d, %d"%(len(UnevolvedList),len(EvolvedList)))
        (PopulationDict,
         EvolvedList,
         UnevolvedList,
         EndConditionMet,
         ChildPest) = CoreSim.RoughFront(
            DomainHeight,DomainWidth,mutprob,fitness,PopulationDict,
            UnevolvedList,EvolvedList,Radius,ObsCenterList,ObsCenterDict,
            EndCondition)
        SimStep += 1


        ##########################################
        ###EndConditions##########################
        ##########################################
        ObstacleX = ObsCenterList[0][0]
        ObstacleY = ObsCenterList[0][1]

        if ((ChildPest.GetPos()[0]-ObstacleX)**2 + 
                (ChildPest.GetPos()[1]-ObstacleY)**2 < (Radius/2)**2):
            EndConditionMet = True
            print("Obstacle Got Infected.")
            ENDREASON = 1

        elif (len(EvolvedList) + len(UnevolvedList)) == 0:
            EndConditionMet = True
            print("Obstacle Not Infected.")
            ENDREASON = 2
        

    #Now, record the highest point
    return [ObsCenterList,ENDREASON]






#############################################################################
#############################################################################
#############################################################################
###Args for repeats
ParamPairsFileName = str(sys.argv[1])
ParamPairID = int(sys.argv[2])

with np.load(os.path.join(ParamPairsFileName)) as data:
    DataMatrix = data['DataMatrix']
    SaveFile = data['SaveFileDirName']
    randomseed = int(data['randomseed'])

random.seed(randomseed + ParamPairID)

fitness = DataMatrix[ParamPairID][0]
mutprob = np.exp(DataMatrix[ParamPairID][1])

print("Fitness",fitness)
print("MutProb",mutprob)

SpecificSaveFile = (str(SaveFile) + '/Fitness_%0.3f_MutProb_'%(fitness) + 
    '{:.1E}'.format(mutprob))

if os.path.isdir(SpecificSaveFile):
    os.rmdir(SpecificSaveFile)

os.mkdir(SpecificSaveFile)


#############################################################################
#############################################################################
#############################################################################

simstarttime = time.time()
####################
####################
####################

DataList = []

#Create Basic Background
#Initialisation

Centerposy = int(BottomBuffer + Radius)
Centerposx = int(DomainWidth/2)

ObsCenterList,ObsCenterDict = CoreSim.SingleObstacle(Radius,DomainHeight,DomainWidth,Centerposx,Centerposy)

InitialPopulationDict,InitialUnevolvedList,InitialEvolvedList = CoreSim.Initialise(0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)



for R in range(Repeats):
    print("Repeat",R)
    #Create Copies
    PopulationDict = copy.deepcopy(InitialPopulationDict)
    UnevolvedList = copy.deepcopy(InitialUnevolvedList)
    EvolvedList = copy.deepcopy(InitialEvolvedList)

    DataList.append(InfectionProbability(DomainHeight,DomainWidth,Radius,mutprob,fitness,EndCondition,PopulationDict,UnevolvedList,EvolvedList,ObsCenterList,ObsCenterDict,R))


####################
####################
####################

####################
####################
####################

simendtime = time.time()
print("Simulation Time: " + str(simendtime-simstarttime))

##################################################################################################
##################################################################################################
##################################################################################################

#SAVING
OutputDatafilename = SpecificSaveFile + '/DataFile.npz'
np.savez(OutputDatafilename,
    Repeats=Repeats,
    fitness=fitness,
    mutprob=mutprob,
    DomainHeight=DomainHeight,
    DomainWidth=DomainWidth,
    Radius=Radius,
    DataList=DataList,
    rseed=randomseed + ParamPairID)

