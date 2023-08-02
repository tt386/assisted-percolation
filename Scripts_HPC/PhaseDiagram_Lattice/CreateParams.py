from Params import *
import numpy as np
import random
import sys
import shutil
#############################################################################
#############################################################################
#############################################################################
#Individual Parameters
UpperAreaRatio = 0.9
LowerAreaRatio = 0.2

AreaRatioList = np.linspace(LowerAreaRatio,UpperAreaRatio,20)
FitnessList = np.linspace(0.7,1.04,18)
#FitnessList = np.arange(0.7,1.04,0.02)

PairNum = len(AreaRatioList) * len(FitnessList)
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
ClassicVariables = ('DomainHeight_%d_DomainWidth_%d_Repeats_%d'%
    (DomainHeight,DomainWidth,Repeats))

#Evolution
EvolutionVariables = '_RoughFront'

#Infection
InfectionVariables = '_FlatInfection'

#Obstacles
ObstacleVariables = '_Lattice_Radius_%d'%(Radius)

SaveFileDirName = ('SaveFiles/' + 
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


"""
######
#Method for doing the quickest ones first!

def TheoryLine(A,R):

    eta = np.sqrt((2*np.pi)/(np.sqrt(3)*A)) - 2
    
    F = ((np.sqrt(3))/(2) * (2+eta)/(np.sqrt(eta**2 + 4*eta) + 2*(np.pi/3 - 
        np.arccos(2/(2+eta)))))
    
    return F

SortList = []

for F in reversed(FitnessList):
    for A in reversed(AreaRatioList):
        SortList.append([F-TheoryLine(A,Radius),[A,F]])
        
#print(SortList)

SortedList = sorted(SortList, key=lambda x: x[0])[::-1]


DataMatrix = []
for i in range(len(SortedList)):
    DataMatrix.append(SortedList[i][1])
"""
#############################################################################
DataMatrix = []
for A in AreaRatioList:
    for F in FitnessList:
        DataMatrix.append([A,F])
OutputDatafilename = 'ParamPairs.npz'
np.savez(OutputDatafilename,
    DataMatrix=DataMatrix,
    SaveFileDirName=SaveFileDirName,
    randomseed=randomseed)

#############################################################################
#############################################################################
#############################################################################

"""
OutputDatafilename = 'ParamPairs.npz'
np.savez(OutputDatafilename,
    DataMatrix=DataMatrix,
    SaveFileDirName=SaveFileDirName,
    randomseed=randomseed)
"""
