
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



def DominationProb(
        ObsCenterList,ObsCenterDict,DomainHeight,DomainWidth,Radius,
        LatticeSurfSep,LatticeType,BottomLeftxdisp,BottomLeftydisp,mutprob,
        fitness,EndCondition,Repeat):
    """
    Exectute the rough front eden model for a Lattice patch distribuion,
    and only permit a single mutation to occur. End the simulation when either
    the M population fall to 0 after the mutation, or when the WT pop falls to
    0.

    ARGS:
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    Radius: int
        The radius of the patches.
    LatticeSurfSep: int
        The separation between surfaces of the patches.
    LatticeType: int
        A (currently redundant) orientation of the lattice desired
    BottomLeftxdisp: int
        The bottommost, leftmost patch surface separation from the left edge
        of the system.
    BottomLeftydisp: int
        The bottommost, leftmost patch surface separation from the bottom edge
        of the system.
    mutprob: float
        The probability that a given occupied site is a M.
    fitness: float
        The fitness of M: the ratio to which a M is chosen relative to a WT.
    EndCondition: int
        This is the option corresponding to the condition for which we will
        cease executing the rough front Eden Model.
    Repeat: int
        The enumerated repeat value of this repeat.

    RETURNS:
    int:
        The separation between surfaces of the patches.
    float:
        The measured ratio of the patch sites to the whole system.
    bool:
        If the System ended due to M dominating the front
    bool:
        If the process ended due to no M on the front after mutation occures.
    int:
        The number of WT in the end-state
    int: 
        The number of M in the end-state
    """

    #Initialisation
 
    print("ObstacleList: " + str(ObsCenterList))
 
    (InitialPopulationDict,
    InitialUnevolvedList,
    InitialEvolvedList) = CoreSim.Initialise(
        0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)

    #####################
    ##For Single Mutant##
    #####################

    #A repeat is only valid of these many mutants are on the front at a time
    MUTANTNUMBERTHRESHOLD = 2.*Radius       
    
    #Is this MUTANTNUMBERTHRESHOLD met yet?
    MUTANT_THRESHOLD_MET = False           
    
    #Variable for the mutation probability
    current_MutationProb = 0              
    
    #The height at which mutation shall be turned
    MUTATION_HEIGHT = 100

    Extinction = False
    Domination = False


    REPEAT_COUNT = 50

    while (not MUTANT_THRESHOLD_MET) and (REPEAT_COUNT >= 0):   
        print("Repeat %d, sub-repeat %d"%(Repeat,REPEAT_COUNT))
        REPEAT_COUNT -= 1

        PopulationDict = copy.deepcopy(InitialPopulationDict)
        UnevolvedList = copy.deepcopy(InitialUnevolvedList)
        EvolvedList = copy.deepcopy(InitialEvolvedList)


        print("Initial Length of Evolved List: " + str(len(EvolvedList)))
        print("Initial Length of Unevolved List: " + str(len(UnevolvedList)))

        current_MutationProb = 0
        MUTATION_HEIGHT_MET = False
        Mutant_Occurred = False
        MaxMutNum = 0

    
        #Main Sim
        EndConditionMet = False
        SimStep = 0
        while EndConditionMet == False:
            ##########################################
            ###CORE###################################
            ##########################################
            if (SimStep % 10000) == 0:
                print("SimStep = " + str(SimStep))
                print("#Wt, #M: %d, %d"%(len(UnevolvedList),len(EvolvedList)))
            (PopulationDict,
            EvolvedList, 
            UnevolvedList, 
            EndConditionMet, 
            ChildPest) = CoreSim.RoughFront(
                DomainHeight,DomainWidth,current_MutationProb,fitness,
                PopulationDict,UnevolvedList,EvolvedList,Radius,ObsCenterList,
                ObsCenterDict,EndCondition)
            SimStep += 1


            #And end condition commented out, as generally I don't want any extension
            """
            #Extend System if Necessary
            if (DomainHeight - ChildPest.GetPos()[1]) < (np.sqrt(3)/2)*Radius:
                (ObsCenterList,
                ObsCenterDict,
                DomainHeight) = CoreSim.ExtendSystem(
                    Radius,AreaFraction,DomainHeight,DomainWidth,
                    ObsCenterList,ObsCenterDict,ExtensionHeight)
                print("SYSTEM EXTENDED TO HEIGHT: " + str(DomainHeight))
            """
            ##########################################
            ##########################################
            ##########################################


            ##########################################
            ###EndConditions##########################
            ##########################################
            #if both lists become empty, then there's been a premature end 
            #due to WT being headed off by obsstacles.
            if (len(UnevolvedList) + len(EvolvedList)) == 0:
                EndConditionMet = True
                print("Ended Without mutation: WT cut off by pesticide.")

            if ChildPest.GetPos()[1] >= 10000-1:
                EndConditionMet = True
                print("Ended because system got way too big.")


            if not EndConditionMet:
                #If mutation has't occurred yet, check if you're at the right
                # height to turn it on
                if not Mutant_Occurred:
                    if not MUTATION_HEIGHT_MET:
                        if UnevolvedList[-1][0] > MUTATION_HEIGHT:
                            current_MutationProb = copy.deepcopy(mutprob)
                            MUTATION_HEIGHT_MET = True

                    if MUTATION_HEIGHT_MET:
                        if len(EvolvedList) > 0:
                            Mutant_Occurred = True
                            current_MutationProb = 0
            
                #If mutation occurred, check if the mutants stay in existence
                # and if the righ tnumber form before extinction/domination
                else:
                    MutNum = len(EvolvedList)
                    if MutNum > MaxMutNum:
                        MaxMutNum = MutNum
            
                    if MutNum == 0:
                        if MaxMutNum < MUTANTNUMBERTHRESHOLD:
                            EndConditionMet = True

                        else:
                            EndConditionMet = True
                            MUTANT_THRESHOLD_MET = True
                            Extinction = True

                    elif len(UnevolvedList) == 0:
                        EndConditionMet = True
                        MUTANT_THRESHOLD_MET = True
                        Domination = True
        print("Domination, Extinction: %d, %d"%(Domination,Extinction))
    print("Repeat %d, subrepeat %d, Dom %d"%(Repeat,REPEAT_COUNT,Domination))

    #Finally, measure endstate area fraction
    MeasuredAreaFraction = CoreSim.MeasureAreaFractionLattice(
        DomainWidth,DomainHeight,Radius,ObsCenterList,ObsCenterDict,
        BottomLeftydisp)

    print("Measured Area Fraction %0.4f"%(MeasuredAreaFraction))

    return (LatticeSurfSep,MeasuredAreaFraction,int(Domination),
            int(Extinction),len(UnevolvedList),len(EvolvedList))





#############################################################################
#############################################################################
#############################################################################
###Args for repeats
#Pull out the parameters from the files based on the argument to the script,
# then create a suitable savefile directory

ParamPairsFileName = str(sys.argv[1])
ParamPairID = int(sys.argv[2])

with np.load(os.path.join(ParamPairsFileName)) as data:
    DataMatrix = data['DataMatrix']
    SaveFile = data['SaveFileDirName']
    randomseed = int(data['randomseed'])

random.seed(randomseed + ParamPairID)

AreaRatio = float(DataMatrix[ParamPairID][0])
fitness = float(DataMatrix[ParamPairID][1])

LatticeSurfSep = np.ceil(Radius * np.sqrt(2*np.pi/(np.sqrt(3) * AreaRatio))
     * (1. - 2 * np.sqrt(np.sqrt(3) * AreaRatio/(2.* np.pi))))

print("LatticeSurfSep: %0.5f, Mutant Fitness: %0.3f"%(LatticeSurfSep,fitness))

SpecificSaveFile = (str(SaveFile) + 
    '/LatticeSurfSep_%0.5f_Fitness_%0.5f'%(LatticeSurfSep,fitness))

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

#Initialisation
ObsCenterList,ObsCenterDict,DomainWidth = CoreSim.LatticeObstacles(
    Radius,DomainHeight,DomainWidth,LatticeType,LatticeSurfSep,
    BottomLeftxdisp,BottomLeftydisp)


DataList = []
for R in range(Repeats):
    DataList.append(
    DominationProb(ObsCenterList,ObsCenterDict,DomainHeight,DomainWidth,
    Radius,LatticeSurfSep,LatticeType,BottomLeftxdisp,BottomLeftydisp,mutprob,
    fitness,EndCondition,R))


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
    AreaRatio=AreaRatio,
    LatticeSurfSep=LatticeSurfSep,
    DataList=DataList,
    BottomLeftydisp=BottomLeftydisp,
    BottomLeftxdisp=BottomLeftxdisp,
    ObsCenterList=ObsCenterList,
    rseed=randomseed + ParamPairID)

