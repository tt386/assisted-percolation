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



def CapHeight(
        DomainHeight,DomainWidth,Radius,mutprob,fitness,EndCondition,
        PopulationDict,UnevolvedList,EvolvedList,ObsCenterList,ObsCenterDict,
        SeededPos,R):
    """
    Run the rough front eden model for a single patch, in which mutation is
    seeded to occur at a single point on the circumfrence of the patch.
    Measure the highest point of the M population before the system ends.

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


    RETURNS:
    list of tuples:
        List of coordinates of the patches in the system.
    int:
        The y-position of the maximum extent of the M population
    int:
        A signal as to why the simulation finished, so we can check if this
        should be included as a repeat.
    """

    ENDREASON = 0   #The reason why the system ends


    HighestMutant = 0

    MUTANTSTARTED = False

    #Main Sim
    EndConditionMet = False
    SimStep = 0
    while EndConditionMet == False:
        #####################################################################
        ###CORE##############################################################
        #####################################################################
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

        #Extend System if Necessary
        if (DomainHeight - ChildPest.GetPos()[1]) < (np.sqrt(3)/2)*Radius:
            ObsCenterList,ObsCenterDict,DomainHeight = CoreSim.ExtendSystem(
                Radius,AreaFraction,DomainHeight,DomainWidth,ObsCenterList,
                ObsCenterDict,ExtensionHeight)
            print("SYSTEM EXTENDED TO HEIGHT: " + str(DomainHeight))
        #####################################################################
        #####################################################################
        #####################################################################


        #Seed mutation
        if ChildPest.GetPos() == SeededPos:
            NewPest = CoreSim.Pest(
                ChildPest.GetPos(),ChildPest.GetParentPos(),0,
                ChildPest.GetGenotype(),1,Radius,ObsCenterList,ObsCenterDict,
                DomainHeight,DomainWidth,PopulationDict)
            NewPest.UpdateEmptyNeighbours(PopulationDict)
    
            PopulationDict[NewPest.GetPos()] = NewPest
            if ChildPest.GetPos() in UnevolvedList:
                UnevolvedList.remove(ChildPest.GetPos())
            if NewPest.GetActive():
                EvolvedList.append(NewPest.GetPos())

            MUTANTSTARTED = True


        #Updating max cap height
        if ChildPest.GetMut():
            if ChildPest.GetPos()[1] > HighestMutant:
               HighestMutant = ChildPest.GetPos()[1] 
            

        #####################################################################
        ###EndConditions#####################################################
        #####################################################################
        #if both lists become empty, then there's been a premature end due 
        #to WT being headed off by obsstacles.
        if (len(UnevolvedList) + len(EvolvedList)) == 0:
            EndConditionMet = True
            print("Ended Without mutation: WT cut off by pesticide.")
            ENDREASON = 1

        elif len(UnevolvedList) == 0:
            EndConditionMet = True
            print("Mutants Dominated the System")
            ENDREASON = 2


        if ChildPest.GetPos()[1] > 10000:
            EndConditionMet = True
            print("Ended because system got way too big.")
            ENDREASON = 3


        if MUTANTSTARTED:
            if len(EvolvedList) == 0:
                EndConditionMet = True
                print("Ended because mutants died")
                ENDREASON = 4
        

    #Now, record the highest point
    return [ObsCenterList,HighestMutant,ENDREASON]






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

SeededPos = DataMatrix[ParamPairID]

print("SeededPos",SeededPos)

SeededPos = tuple(SeededPos)

SpecificSaveFile = str(SaveFile) + '/xpos_%d_ypos_%d'%(SeededPos[0],
    SeededPos[1])

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

Centerposy = BottomBuffer + Radius
Centerposx = DomainWidth/2

ObsCenterList,ObsCenterDict = CoreSim.SingleObstacle(
    Radius,DomainHeight,DomainWidth,Centerposx,Centerposy)

(InitialPopulationDict,
InitialUnevolvedList,
InitialEvolvedList) = CoreSim.Initialise(
    0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)



for R in range(Repeats):
    print("Repeat",R)
    #Create Copies
    PopulationDict = copy.deepcopy(InitialPopulationDict)
    UnevolvedList = copy.deepcopy(InitialUnevolvedList)
    EvolvedList = copy.deepcopy(InitialEvolvedList)

    DataList.append(
        CapHeight(DomainHeight,DomainWidth,Radius,mutprob,fitness,
        EndCondition,PopulationDict,UnevolvedList,EvolvedList,
        ObsCenterList,ObsCenterDict,SeededPos,R))


####################
####################
####################

####################
####################
####################

simendtime = time.time()
print("Simulation Time: " + str(simendtime-simstarttime))

#############################################################################
#############################################################################
#############################################################################

#SAVING
OutputDatafilename = SpecificSaveFile + '/DataFile.npz'
np.savez(OutputDatafilename,
    Repeats=Repeats,
    fitness=fitness,
    mutprob=mutprob,
    DomainHeight=DomainHeight,
    DomainWidth=DomainWidth,
    Radius=Radius,
    SeededPos=list(SeededPos),
    DataList=DataList,
    rseed=randomseed + ParamPairID)

