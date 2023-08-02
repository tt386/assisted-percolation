from Params import *
import numpy as np
import random
import sys
import shutil
#############################################################################
#############################################################################
#############################################################################

#Create the list of sites next to the obstacle
FitnessList = [0.1,0.8,0.9]


MinMutRate = 1e-6
MaxMutRate = 1e-2

MutationRateList = np.linspace(np.log(MinMutRate),np.log(MaxMutRate),50)
                   
print("FitnesList",FitnessList)
print("MutRateList",MutationRateList)
ParamNum = len(FitnessList) * len(MutationRateList)
print("Number of Parameters (therefore, last number in array)" + 
    str(ParamNum))


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
ClassicVariables = ('Height_%d_Width_%d_Repeats_%d_Radius_%d_BottomBuffer_%d'
    %(DomainHeight,DomainWidth,Repeats,Radius,BottomBuffer))

#Evolution
EvolutionVariables = '_RoughFront'

#Infection
InfectionVariables = '_FlatInfection'

#Obstacles
ObstacleVariables = '_Random_Radius_%d'%(Radius)

SaveFileDirName = ('SaveFiles/' + 
    ClassicVariables + 
    InfectionVariables + 
    EvolutionVariables)

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
for F in FitnessList:
    for M in MutationRateList:
        DataMatrix.append([F,M])
OutputDatafilename = 'ParamPairs.npz'
np.savez(OutputDatafilename,
    DataMatrix=DataMatrix,
    SaveFileDirName=SaveFileDirName,
    randomseed=randomseed)

