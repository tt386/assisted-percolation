from Params import *

import numpy as np
import copy

import matplotlib.pyplot as plt
import matplotlib as mpl
from textwrap import wrap

import time

import os

import random

import sys

from argparse import ArgumentParser
parser = ArgumentParser(description='Flat Front')
parser.add_argument('-f','--F',help='The fitness of the pest')
parser.add_argument('-R','--R',help='Repeats')
parser.add_argument('-d','--directory',
    help='The directory where we save the data and images')
parser.add_argument('-r','--randomseed',type=int,required=False,help='The random seed')
args = parser.parse_args()

#Set variables
F = fitness
if args.F is not None:
    F = float(args.F)

R = Repeats
if args.R is not None:
    R = int(args.R)

SaveDirectory = SaveFileDir#"SaveFiles/Fitness_%0.3f_Repeats_%d"%(F,R)
if args.directory is not None:
    SaveDirectory = str(args.directory)

if not os.path.isdir(SaveDirectory):
    os.mkdir(SaveDirectory)

#Set random seed
if args.randomseed is None:
    print("No seed provided")
    randomseed = random.randrange(sys.maxsize)

else:
    randomseed = int(args.randomseed)
    print("seed provided:",randomseed)

random.seed(randomseed)


with open(SaveDirectory+'/FLAT_randomseed.txt', 'w') as f:
    f.write('seed:'+str(randomseed))

##############################################################################
##############################################################################
##############################################################################



HeightList = []
WidthList = []
"""
R = int(args.R)#1000000

F = float(args.F)#0.95
"""
p_grow = (F/(1+F))**2
p_shrink = 1./(1+F)**2
p_same = 2*F/(1+F)**2

starttime = time.time()

for r in range(R):
    Height = 1
    Width = 1 
    MaxWidth = 1
    
    while Width > 0:
        randomnum = random.uniform(0,1)
        
        if randomnum < p_grow:
            Height += 1
            Width += 1
            if Width > MaxWidth:
                MaxWidth = copy.deepcopy(Width)
        elif (randomnum>p_grow) and (randomnum<(p_grow+p_shrink)):
            Width -= 1
            if Width > 0:
                Height += 1
        else:
            Height += 1
            
    HeightList.append(Height)
    WidthList.append(MaxWidth)
    
endtime = time.time()
print(endtime-starttime)


OutputDatafilename = SaveDirectory + '/FlatFrontDataFile.npz'
np.savez(OutputDatafilename,
            F=F,
            R=R,
            HeightList=HeightList,
            WidthList=WidthList)






#Plotting
Heightxlist = np.arange(max(HeightList)) + 1
Heightfreqlist = np.zeros(len(Heightxlist))

for i in range(len(HeightList)):
   Heightfreqlist[HeightList[i]-1] += 1

Heightfreqlist = Heightfreqlist/sum(Heightfreqlist)
 
plt.figure(1)
plt.loglog(Heightxlist,Heightfreqlist)
plt.grid(True)
plt.xlim(0,200)
plt.savefig(SaveDirectory + '/FlatFrontHeightDistOnly.png')
plt.close()


