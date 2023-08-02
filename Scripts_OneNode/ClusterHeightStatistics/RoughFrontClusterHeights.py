
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


##############################################################################
##############################################################################
##############################################################################



def ClusterHeight(
    DomainHeight,DomainWidth,Radius,mutprob,fitness,EndCondition,
    PopulationDict,UnevolvedList,EvolvedList,ObsCenterList,ObsCenterDict,
    SeededPos,R):
    """
    Keep track of the maximum and minimum extents of a cluster.

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
    UnevolvedList: list of tuples
        The coordinates of WT pests on the front.
    EvolvedList: list of tuples
        The coordinated of M pests on the front.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    SeededPos: tuple of ints
        The coordinates at which a Pest occupying this point shall become a M.
    R: int
        The enumerated repeat of the system.

    RETURNS:
    list of ints:
        [The y-pos of the highest M, the y-pos of the lowest M,
        the REASON the system ends (to provide an indication of
        if the repeat out to be counted)]
    """
    

    ENDREASON = 0   #The reason why the system ends


    HighestMutant = SeededPos[1]
    LowestMutant = SeededPos[1]

    MUTANTSTARTED = False

    #Main Sim
    EndConditionMet = False
    SimStep = 0
    while EndConditionMet == False:
        ##########################################
        ###CORE###################################
        ##########################################
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
        ##########################################
        ##########################################


        #Seed mutation
        if ChildPest.GetPos() == SeededPos:
            NewPest = CoreSim.Pest(
                ChildPest.GetPos(),ChildPest.GetParentPos(),0,0,1,Radius,
                ObsCenterList,ObsCenterDict,DomainHeight,DomainWidth,
                PopulationDict)

            NewPest.UpdateEmptyNeighbours(PopulationDict)
    
            PopulationDict[NewPest.GetPos()] = NewPest
            if ChildPest.GetPos() in UnevolvedList:
                UnevolvedList.remove(ChildPest.GetPos())
            if NewPest.GetActive():
                EvolvedList.append(NewPest.GetPos())

            MUTANTSTARTED = True


        #Updating max cap height
        if ChildPest.GetMut():
            childypos = ChildPest.GetPos()[1]

            if childypos > HighestMutant:
               HighestMutant = childypos
            
            if childypos < LowestMutant:
                LowestMutant = childypos

        ##########################################
        ###EndConditions##########################
        ##########################################
        #If both lists become empty, then there's been a premature end due to
        #WT being headed off by obstacles.
        
        if (len(UnevolvedList) + len(EvolvedList)) == 0:
            EndConditionMet = True
            print("Ended Without mutation: WT getting cut off by pesticide.")
            ENDREASON = 1

        elif len(UnevolvedList) == 0:
            EndConditionMet = True
            print("Mutants Dominated the System")
            ENDREASON = 2


        if ChildPest.GetPos()[1] > 1000:
            EndConditionMet = True
            print("Ended because system got way too big.")
            ENDREASON = 3


        if MUTANTSTARTED:
            if len(EvolvedList) == 0:
                EndConditionMet = True
                print("Ended because mutants died")
                ENDREASON = 4
        
    print("ClusterHeight: ",HighestMutant-LowestMutant + 1)
    #Now, record the highest point
    return [HighestMutant,LowestMutant,ENDREASON]






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

ObsCenterList,ObsCenterDict = CoreSim.Obstacles(Radius,0,DomainHeight,DomainWidth)

InitialPopulationDict,InitialUnevolvedList,InitialEvolvedList = CoreSim.Initialise(
    0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)



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



for R in range(Repeats):
    print("Repeat",R)
    #Create Copies
    PopulationDict = copy.deepcopy(InitialPopulationDict)
    UnevolvedList = copy.deepcopy(InitialUnevolvedList)
    EvolvedList = copy.deepcopy(InitialEvolvedList)

    DataList.append(ClusterHeight(
        DomainHeight,DomainWidth,Radius,mutprob,fitness,EndCondition,
        PopulationDict,UnevolvedList,EvolvedList,ObsCenterList,ObsCenterDict,
        SeededPos,R))


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
OutputDatafilename = SaveFileDir + '/RoughFrontDataFile.npz'
np.savez(OutputDatafilename,
    Repeats=Repeats,
    fitness=fitness,
    mutprob=mutprob,
    DomainHeight=DomainHeight,
    DomainWidth=DomainWidth,
    Radius=Radius,
    SeededPos=list(SeededPos),
    DataList=DataList,
    randomseed=randomseed)

with open(SaveFileDir+'/ROUGH_randomseed.txt', 'w') as f:
    f.write('seed:'+str(randomseed))

