import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt


import sys
sys.path.insert(0,'../../CoreFunctions')

import TheoreticalPredictionFuns


import math
import numpy as np
import random

from numpy.linalg import matrix_power

from scipy import stats
from scipy.stats import chisquare

#Plotting aspects
from matplotlib.pyplot import cm


#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
#mpl.rcParams['agg.path.chunksize'] = 10000

from textwrap import wrap

import copy
from copy import deepcopy

import matplotlib.ticker as mtick

def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))

from scipy import stats


#############################################################################
##ArgParse###################################################################
#############################################################################
import os.path

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The File %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

from argparse import ArgumentParser

parser = ArgumentParser(description='Monte-Carlo Clusters')
parser.add_argument('-d','--directory',help='The directory of the data')

args = parser.parse_args()

#Find all the directories
templist = os.listdir(args.directory)
print(templist)

dirlist = []

for i in range(len(templist)):
    if os.path.isdir(args.directory + '/' + templist[i]):
        print("Is directory!")
        npzlist = os.listdir(args.directory + '/' + templist[i])
        for j in range(len(npzlist)):
            if npzlist[j].endswith(".npz"):
                dirlist.append(templist[i])

FitnessList = []
AreaFractionList = []
MeanTimeList = []
MeanDominationProb = []

TheoryAreaFractionList = []
TheoryFitnessList = []

qList = []

for d in dirlist:
    filelist = os.listdir(args.directory + "/" + d)
    for names in filelist:
        if names.endswith(".npz"):
            filename = names
            print("Found File!")

    with np.load(os.path.join(args.directory,d,filename)) as data:
        #rownum = data['rownum']
        #colnum = data['colnum']
        Repeats = data['Repeats']
        fitness = data['fitness']
        mutprob = data['mutprob']
        AreaFraction = data['AreaFraction']
        Radius = data['Radius']
        DataList = data['DataList'] #AreaFraction,MeasuredAreaFraction,
                                    #int(Domination),int(Extinction),Time,
                                    #SavedUnevolvedList,SavedEvolvedList

    print("AreaFraction %0.3f, Fitness %0.6f"%(AreaFraction,fitness))

    if AreaFraction not in TheoryAreaFractionList:
        TheoryAreaFractionList.append(AreaFraction)

    if fitness not in TheoryFitnessList:
        TheoryFitnessList.append(fitness)

    MeanAreaFraction = 0
    MeanDomination = 0
    MeanTime = 0
    TotalRepeats = 0
    for R in range(Repeats):
        Domination = DataList[R][2]
        Extinction = DataList[R][3]
        #Time = DataList[R][4]
        if Domination or Extinction:

            MeanAreaFraction += DataList[R][0]
            MeanDomination += float(DataList[R][2])
            TotalRepeats += 1
 
    FitnessList.append(float(fitness))
    if TotalRepeats == 0:
        TotalRepeats = 1
    AreaFractionList.append(MeanAreaFraction/TotalRepeats)
    MeanTimeList.append(MeanTime/TotalRepeats)
    MeanDominationProb.append(MeanDomination/TotalRepeats)
    print("Total_Repeats: " + str(TotalRepeats))

#################################
###Theory########################
#################################
TheoryAreaFractionList = []
TheoryFitnessList = sorted(TheoryFitnessList)
TheoryFitnessList = np.linspace(min(TheoryFitnessList),1,100)
for Fitness in TheoryFitnessList:
    a = TheoreticalPredictionFuns.WaveObstacleCapHeight(
        -np.pi/2,Fitness,Radius) 
    #Space contraction consideratons
    
    AreaFraction_LIMIT = 0.676#8

    TheoryAreaFrac=1.-(1.-AreaFraction_LIMIT)**((2*Radius)/(2*Radius+a))

    TheoryAreaFractionList.append(TheoryAreaFrac)
    
    print("Fitness, a, AreaFraction: %0.3f, %0.3f, %0.3f"%
        (Fitness,a,TheoryAreaFractionList[-1]))

#############################################################################
###Plotting##################################################################
#############################################################################


#######
TimeVar = False
#######

A = np.array(AreaFractionList)
F = np.array(FitnessList)
if not TimeVar:
    D = np.array(MeanDominationProb)
else:
    D = np.array(MeanTimeList)
#####################
###pcolor approach

#Create the matrix
AreaSet = set(AreaFractionList)
print(FitnessList)
FitnessSet = set(FitnessList)

print("AreaSet:",AreaSet)
print("Fitness Set:",FitnessSet)

SortedAreaList = list(AreaSet)
SortedAreaList.sort()
dA = (SortedAreaList[1]-SortedAreaList[0])/2

SortedFitnessList = list(FitnessSet)
SortedFitnessList.sort()
dF = (SortedFitnessList[1]-SortedFitnessList[0])/2

print("AreaList:",SortedAreaList)
print("FitnessList:",SortedFitnessList)

ProbMatrix = np.zeros((len(SortedFitnessList),len(SortedAreaList)))

num = 0
for row in range(len(ProbMatrix)):
    for col in range(len(ProbMatrix[row])):
        targetArea = SortedAreaList[col]
        targetFitness = SortedFitnessList[row]

        for i in range(len(MeanDominationProb)):
            if (abs(AreaFractionList[i] - targetArea)<dA) and (abs(FitnessList[i]-targetFitness)<dF):
                ProbMatrix[row][col] = MeanDominationProb[i]
                #print(i)
                num += 1
                break

            if i == len(MeanDominationProb):
                print("Shouldn't happen")
#np.flip(ProbMatrix,0)
#np.flip(ProbMatrix,1)
print("Number of times it was occupied:",num)
print("Matrix total element number:",ProbMatrix.size)
print(ProbMatrix)

SortedFitnessList.append(SortedFitnessList[-1]+SortedFitnessList[-1]-SortedFitnessList[-2])
SortedAreaList.append(SortedAreaList[-1]+SortedAreaList[-1]-SortedAreaList[-2])

fig=plt.figure(2)
ax = fig.add_subplot(1,1,1)

cp = ax.tricontourf(A.ravel(), F.ravel(), D.ravel(),10,cmap='coolwarm')


ax.pcolor(np.asarray(SortedAreaList)-dA,np.asarray(SortedFitnessList)-dF,np.array(ProbMatrix),cmap='coolwarm')



ax.set_xticks([0,0.25,0.5,0.75])
ax.set_xticklabels(
    ['$0.00$',r'$0.25$',r'$0.5$',r'$0.75$'],
    fontsize=7)

ax.set_yticks([0.7,0.8,0.9,1])
ax.set_yticklabels(
    ['$0.7$',r'$0.8$',r'$0.9$',r'$1.0$'],
    fontsize=7)



ax.set_ylabel(r'mutant fitness $F$',fontsize=8)
ax.set_xlabel(r'patch area ratio $\phi$',fontsize=8)

print("cp",cp)
cbar = plt.colorbar(cp)
cbar.set_label(label=r'prob. of M dominating $P_D$',size=8)
cbar.ax.tick_params(labelsize=7)

ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.savefig(str(args.directory) + '/testa.png')

plt.figure(2).set_size_inches(9.00000/2.54, 6.000000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
plt.figure(2).axes[0].set_position([0.1, 0.1, 0.85, 0.85])


#Adding coordinates of interest
"""
plt.scatter([0.4],[0.9],s=70,color='black',zorder=10,marker='x')
plt.scatter([0.4],[0.98],s=70,color='black',zorder=10)
plt.scatter([0.7],[0.9],s=70,color='black',zorder=10,marker='+')
"""

plt.scatter([0.4],[0.9],s=50,color='white',zorder=10,marker='X',edgecolors='black')
plt.scatter([0.4],[0.98],s=50,color='white',zorder=10,edgecolors='black')
plt.scatter([0.7],[0.9],s=50,color='white',zorder=10,marker='P',edgecolors='black')

plt.tight_layout()

plt.savefig(str(args.directory) + '/pcolor_Domination_NoTheory.png')
ax.plot(TheoryAreaFractionList,TheoryFitnessList,color='black',linewidth=2)

ax.plot([0.676,0.676],[0.6,1.1],'--',color='black',linewidth=1)
ax.plot([-0.5,0.9],[1,1],'--',color='black',linewidth=1)

plt.ylim(min(SortedFitnessList)-dF,max(SortedFitnessList)-dF)
plt.xlim(min(SortedAreaList)-dA,max(SortedAreaList)-dA)

print("TheoryAreaFractionList",TheoryAreaFractionList)
print("TheoryFitnessList",TheoryFitnessList)

plt.savefig(str(args.directory) + '/pcolor_Domination.png')
plt.savefig(str(args.directory) + '/pcolor_Domination.pdf')
plt.savefig(str(args.directory) + '/pcolor_Domination.eps')
plt.close()



print("Finished")
