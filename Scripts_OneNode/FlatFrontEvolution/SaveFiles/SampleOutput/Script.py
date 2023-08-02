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
for x in range(DomainWidth):
    newoutputrow = [0,0,0,0]
    newoutputrow[0] = 0                     #time
    newoutputrow[1] = x                     #x
    newoutputrow[2] = 0                     #y
    newoutputrow[3] = 0                     #Species/colour
    OutputData.append(newoutputrow)


#Domain
System = np.zeros((DomainHeight,DomainWidth))

simstarttime = time.time()


#Below is the main loop executed amongst al the scripts. The system will loop
#util a predefined end condition is met, but by including it in every script
#there is the option to define specialist end conditions, such as stopping
#when there are no more mutants etc. We could also define things such as 
#turning the mutation rate on and off given certain conditions.

EndConditionMet = False
for y in range(1,DomainHeight):
    for x in range(DomainWidth):

        #Determine the possible parents of the selected site
        leftparent = 0
        rightparent = 0

        if x%2 == 1:
            leftparent = System[y-1][x]
            rightparent = System[y-1][(x+1)%DomainWidth]

        else:
            leftparent = System[y-1][x-1]
            rightparent = System[y-1][x]

        #Given the possible parents, decide how may infection spread
        if leftparent == rightparent == 1:
            System[y][x] = 1

        elif leftparent == rightparent == 0:
            randomnum = random.uniform(0,1)
    
            if randomnum <= mutprob:
                System[y][x] = 1

        else:
            randomnum = random.uniform(0,1)

            if randomnum <= fitness/(1+fitness) + mutprob/(1+fitness):
                System[y][x] = 1
        

        #Add output to the list
        newoutputrow = [0,0,0,0]
        newoutputrow[0] = SimStep                     #time
        newoutputrow[1] = x                           #x
        newoutputrow[2] = y                           #y
        newoutputrow[3] = System[y][x]                #color
        OutputData.append(newoutputrow)

    SimStep +=1


simendtime = time.time()

print("Time taken for sim: " + str(simendtime-simstarttime))






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
    DataSizeList=DataSizeList,
    randomseed=randomseed)


with open(SaveFileDir+'/randomseed.txt', 'w') as f:
    f.write('seed:'+str(randomseed))


