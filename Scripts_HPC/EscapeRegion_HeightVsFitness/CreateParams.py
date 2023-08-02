from Params import *
import numpy as np
import random
import sys
import shutil
#############################################################################
#############################################################################
#############################################################################
"""
#Create the list of sites next to the obstacle
Centerposy = BottomBuffer + Radius
Centerposx = int(DomainWidth/2)
NextToResistanceList = []
for X in range(Centerposx - Radius - 10, Centerposx + Radius + 10):
    for Y in range(Centerposy - int(Radius*2/np.sqrt(3)) - 10,
             Centerposy + int(Radius*2/np.sqrt(3)) + 10):

        if (Centerposx - X)**2 + 0.75*(Centerposy - Y) **2 >= Radius**2:
            NeighbourList = []
            for y in range(Y-1,Y+2):
                if (y % 2) == 0:
                    ax = cx = ex = (X-1)%(DomainWidth)
                    bx = fx = X
                    dx = (X+1)%(DomainWidth)
                else:
                    cx = (X-1)%(DomainWidth)
                    ax = ex = X
                    bx = dx = fx = (X+1)%(DomainWidth) 
                
                ay = by = y+1
                cy = dy = y
                ey = fy = y-1
                    
            
                NeighbourList.append((ax,ay))
                NeighbourList.append((bx,by))
                NeighbourList.append((cx,cy))
                NeighbourList.append((dx,dy))
                NeighbourList.append((ex,ey))
                NeighbourList.append((fx,fy)) 

            for n in NeighbourList:
                x = n[0]
                y = n[1]
                if (((Centerposx - x)**2 + 0.75*(Centerposy - y)**2) < 
                        Radius**2):
                    if [X,Y] not in NextToResistanceList:                    
                        NextToResistanceList.append([X,Y])
                    
print("NextToResistanceList",NextToResistanceList)
ParamNum = len(NextToResistanceList)
print("Number of Parameters (therefore, last number in array)" + 
    str(ParamNum))
"""
fitnesslist = np.linspace(Min_fitness,Max_fitness,num_fitness)
ParamNum = len(fitnesslist)
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
ClassicVariables = ("Height_%d_Width_%d_Repeats_%d_Radius_%d_"
    "Min_Fitness_%0.3f_"
    "Max_Fitness_%0.3f_"
    "num_Fitness_%d"
    %(DomainHeight,DomainWidth,Repeats,Radius,Min_fitness,
    Max_fitness,num_fitness))

#Evolution
EvolutionVariables = '_RoughFront'

#Infection
InfectionVariables = '_FlatInfection'

#Obstacles
ObstacleVariables = '_Random_Radius_%d'%(Radius)

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



DataMatrix = []
for f in fitnesslist:
    DataMatrix.append(f)
OutputDatafilename = 'ParamPairs.npz'
np.savez(OutputDatafilename,
    DataMatrix=DataMatrix,
    SaveFileDirName=SaveFileDirName,
    randomseed=randomseed)

