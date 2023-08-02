from Params import *
import numpy as np
import random
import sys
import shutil
#############################################################################
#############################################################################
#############################################################################
#Individual Parameters
AreaFractionList = np.linspace(0,0.85,18)
FitnessList = np.linspace(0.7,1.02,17)

PairNum = len(AreaFractionList) * len(FitnessList)
print("Number of Parameter Pairs (therefore, last number in array)" + 
    str(PairNum))




#Create the random seed
#################################
###Argparse
#################################
from argparse import ArgumentParser
parser = ArgumentParser(description='Define the (optional) random seed')
parser.add_argument('-r','--randomseed',type=int,required=False,help='The random seed')
args = parser.parse_args()

if args.randomseed is None:
    print("No seed provided")
    randomseed = random.randrange(sys.maxsize-10000)

else:
    randomseed = int(args.randomseed)
    print("seed provided:",randomseed)



#Create Directory for Saving all of the data:
#############################################################################
#############################################################################
#############################################################################
#Do naming now
#SaveFile Name
ClassicVariables = ('DomainHeight_%d_DomainWidth_%d_Repeats_%d'
    %(DomainHeight,DomainWidth,Repeats))

#Evolution
EvolutionVariables = '_RoughFront'

#Infection
InfectionVariables = '_FlatInfection'

#Obstacles
ObstacleVariables = '_Random_Radius_%d'%(Radius)

SaveFileDirName = (
    'SaveFiles/' + 
    ClassicVariables + 
    InfectionVariables + 
    EvolutionVariables + 
    ObstacleVariables)

if not path.exists(SaveFileDirName):
    os.mkdir(SaveFileDirName)
shutil.copy("Params.py",SaveFileDirName)


with open(SaveFileDirName+'/randomseed.txt', 'w') as f:
    f.write('seed:'+str(randomseed))
#############################################################################
#############################################################################
#############################################################################



#Save Parameter Pairs
if os.path.exists("ParamPairs.npz"):
    os.remove("ParamPairs.npz")
    print("Deleted Previous")
else:
    print("The file did not exist in the first place")



DataMatrix = []
for A in AreaFractionList:
    for F in FitnessList:
        DataMatrix.append([A,F])
OutputDatafilename = 'ParamPairs.npz'
np.savez(OutputDatafilename,
    DataMatrix=DataMatrix,
    SaveFileDirName=SaveFileDirName,
    randomseed=randomseed)

