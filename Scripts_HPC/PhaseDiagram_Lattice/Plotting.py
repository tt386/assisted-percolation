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

from textwrap import wrap

import copy
from copy import deepcopy

import matplotlib.ticker as mtick

def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))

from scipy import stats

from scipy import optimize

#############################################################################
#############################################################################
#############################################################################



def WaveObstacleCapHeightROWCONSIDERATION(
        InfectionAngle,b,Fitness,ObstacleRadius,SurfSep):
    """
    Use the principle of least time to calculate the height of the emergence
    region of the M above the patch, for M of a given Fitness, invading
    a patch of radius ObstacleRadius at an angle InfectionAngle. The original
    front of the population is perturbed by the previous layers of patches.


    ARGS:
    InfectionAngle: float
        The angle, as measured in polar coordinates, at which the patch
         becomes invaded.
    b: int
        The offset from the patch center we estimate the peak of the escape
        region to be. We use 0 most of the time.
    Fitness: float
        The fitness of the M sites.
    ObstacleRadius: int
        The radius of the patch.
    SurfSep: int
        The surface separation between nearest neighbour patches.

    RETURNS:
    float:
        Height of the escape region.
    """    

    def function(x,InfectionAngle,b,Fitness,Radius,SurfSep):
        if InfectionAngle > 0:
            return ((1.0/Fitness) * np.sqrt((Radius**2 ) + (x**2) - 
                2*x*Radius*np.sin(InfectionAngle)) + Radius * InfectionAngle-
                 Radius * np.arcsin(Radius/x) - np.sqrt((x)**2 - (Radius)**2))
        else:
            if b!= 0:
                return ((((np.sqrt(3.)/2)*(2*Radius+SurfSep)) -Radius) + 
                    (1.0/Fitness) * np.sqrt((Radius**2 ) + (b**2) + (x**2) - 
                    2*Radius*(x*np.sin(InfectionAngle) + 
                    b*np.cos(InfectionAngle))) - (np.sqrt(SurfSep**2 +
                    4*Radius*SurfSep)+2*Radius*(np.pi/3. - 
                    np.arccos(2.*Radius/(2.*Radius+SurfSep)))) - Radius * 
                    (np.arctan(x/b) - np.arccos(Radius/np.sqrt(x*x+b*b))) - 
                    np.sqrt((x)**2 + b*b - (Radius)**2))
            else:
                return ((((np.sqrt(3.)/2)*(2*Radius+SurfSep)) -Radius) + 
                    (1.0/Fitness) * np.sqrt((Radius**2 ) + (b**2) + (x**2) - 
                    2*Radius*(x*np.sin(InfectionAngle) + 
                    b*np.cos(InfectionAngle))) - (np.sqrt(SurfSep**2 +
                    4*Radius*SurfSep)+2*Radius*(np.pi/3. - 
                    np.arccos(2.*Radius/(2.*Radius+SurfSep)))) - Radius * 
                    (np.pi/2 - np.arccos(Radius/np.sqrt(x*x+b*b))) - 
                    np.sqrt((x)**2 + b*b - (Radius)**2))
    InitialGuess = ObstacleRadius
    try:
        root = optimize.newton(
            function,InitialGuess,args=(InfectionAngle,b,Fitness,
            ObstacleRadius,SurfSep))
    except:
        
        root = 0
    return root

#############################################################################
#############################################################################
#############################################################################



def WaveObstacleExitAngleROWCONSIDERATION(
        InfectionAngle,Fitness,ObstacleRadius,SurfSep):
    """
    Use the principle of least time to calculate the angle around the patch
    perimeter at which the M and WT strains would reach that point at the
    same time, for M of a given Fitness, invading a patch of radius
    ObstacleRadius at an angle InfectionAngle. The original
    front of the population is perturbed by the previous layers of patches.


    ARGS:
    InfectionAngle: float
        The angle, as measured in polar coordinates, at which the patch
         becomes invaded
    Fitness: float
        The fitness of the M sites
    ObstacleRadius: int
        The radius of the patch
    SurfSep: int
        The surface separation between nearest neighbour patches.

    RETURNS:
    float:
        Angle at which WT and M would meat around the patch perimeter
    """


    def function(phi,InfectionAngle,Fitness,Radius,SurfSep):
        return ((((np.sqrt(3.)/2)*(2*Radius+SurfSep)) -
            Radius*np.sin(abs(InfectionAngle))) + (1.0/Fitness) * 
            np.sqrt((Radius*np.cos(phi) - Radius*np.cos(InfectionAngle))**2 +
            (Radius*np.sin(abs(InfectionAngle)) + Radius*np.sin(phi))**2) - 
            (np.sqrt(SurfSep**2 +4*Radius*SurfSep)+2*Radius*(np.pi/3. - 
            np.arccos(2.*Radius/(2.*Radius+SurfSep)))) - Radius*phi)
    
    InitialGuess = np.pi/4
    try:
        root = optimize.newton(
            function,InitialGuess,
            args=(InfectionAngle,Fitness,ObstacleRadius,SurfSep))
    except:
        root = 0
    return root


#############################################################################
#############################################################################
#############################################################################

def LatticeTheoryROW(SepList,Radius):
    """
    Use the least time escape region height and the least time escape from
    the patch height, to find if the straight line connecting these points 
    intersects with an obstacle above. This is repeated for each surface
    separation, exploring through fitness values to output a list of fitness
    corresponding to each separation value.

    ARGS:
    SepList: list of ints
        The list of surface separation values between patches
    Radius: int
        The radius of the patch

    RETURNS:
    list of floats:
        Fitness values corresponding the the seplist values that outline the
        theoretical results
    """
    InfectionAngle = -np.pi/2
    FitnessList = []
    
    
    for Sep in SepList:
        FoundFitness = False
        Fitness = 0.7   #Work up from this value
        while not FoundFitness:
    
            EscapeAngle = WaveObstacleExitAngleROWCONSIDERATION(
                InfectionAngle,Fitness,Radius,Sep)
            
            v = Radius * np.cos(EscapeAngle)
            w = Radius * np.sin(EscapeAngle)
                
            alpha = WaveObstacleCapHeightROWCONSIDERATION(
                InfectionAngle,0,Fitness,Radius,Sep)
            
            
            Cx = (2*Radius + Sep) * np.cos(np.pi/3)
            Cy = (2*Radius + Sep) * np.sin(np.pi/3)
            
            a = (w-alpha)/v
            b = -1
            c = alpha
            
            Cx = (2*Radius + Sep) * np.cos(np.pi/3)
            Cy = (2*Radius + Sep) * np.sin(np.pi/3)
            
            d = abs(a*Cx + b*Cy + c)/np.sqrt(a**2 + b**2)
            
            
            
            if d <= Radius:
                FoundFitness = True
                FitnessList.append(Fitness)
                
   
            elif Fitness > 1:
                FoundFitness = True
                FitnessList.append(1)
            Fitness += 0.001
    return FitnessList


#############################################################################
#############################################################################
#############################################################################


def LatticeTheoryRACE(SepList,Radius):
    """
    For each value in the SepList, find the corresponding fitness value for
    which M and WT would travel from one equator to the next in the same time

    ARGS:
    SepList: list of ints
        The list of surface separation values between patches
    Radius: int
        The radius of the patch
    
    RETURNS:
    list of floats:
        Fitness values corresponding the the seplist values that outline the
        theoretical results
    """
    FitnessList = []
    for Sep in SepList:
        f = (np.sqrt(Sep**2 + 4*Radius*Sep) + 2*Radius * (np.pi/3 - 
            np.arccos(2.*Radius/(2.*Radius + Sep))))
        FitnessList.append((np.sqrt(3)/2.) * (2.*Radius + Sep)/f)

    return FitnessList




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

TheorySepList = []

qList = []


#Lists for checking for holes due to the long time taken to 
#create systems that probably never dominate
UniqueFList = []
UniqueAList = []
CoordList = []

for d in dirlist:
    filelist = os.listdir(args.directory + "/" + d)
    for names in filelist:
        if names.endswith(".npz"):
            filename = names
            print("Found File!")

    with np.load(os.path.join(args.directory,d,filename)) as data:
        Repeats = data['Repeats']
        fitness = data['fitness']
        mutprob = data['mutprob']
        LatticeSurfSep = data['LatticeSurfSep']
        AreaFraction = round(float(data['AreaRatio']),3)
        Radius = data['Radius']
        DataList = data['DataList'] #LatticeSurfSep,MeasuredAreaFraction,
                                    #int(Domination),int(Extinction),
                                    #SavedUnevolvedList,SavedEvolvedList

    print("LatticeSurfSep %0.3f, Fitness %0.6f"%(LatticeSurfSep,fitness))

    #AreaFraction = DataList[0][1]
    if AreaFraction not in TheoryAreaFractionList:
        TheoryAreaFractionList.append(AreaFraction)

    if fitness not in TheoryFitnessList:
        TheoryFitnessList.append(fitness)

    if LatticeSurfSep not in TheorySepList:
        TheorySepList.append(int(LatticeSurfSep))



    #Checks for empty spaces due to slow simulation
    if fitness not in UniqueFList:
        UniqueFList.append(fitness)
    if AreaFraction not in UniqueAList:
        UniqueAList.append(AreaFraction)
    if (AreaFraction,fitness) not in CoordList:
        CoordList.append((AreaFraction,fitness))



    MeanAreaFraction = 0
    MeanDomination = 0
    MeanTime = 0
    TotalRepeats = 0

    INCLUDEUNFINISHED = False

    for R in range(Repeats):
        Domination = DataList[R][2]
        Extinction = DataList[R][3]
        #Time = DataList[R][4]
        if Domination or Extinction:

            MeanAreaFraction += AreaFraction#DataList[R][1]
            MeanDomination += float(DataList[R][2])
            #MeanTime += DataList[R][4]
            TotalRepeats += 1
        
        else:
            if INCLUDEUNFINISHED:

                #Argue that if surivived this long, must be successful

                #Using the given area fraction is close enough generally.
                #Using given makes sure all data points aligned
                MeanAreaFraction += DataList[R][0]
                MeanDomination += 1
                TotalRepeats += 1


    FitnessList.append(float(fitness))
    AreaFractionList.append(MeanAreaFraction/TotalRepeats)
    MeanTimeList.append(MeanTime/TotalRepeats)
    MeanDominationProb.append(MeanDomination/TotalRepeats)
    print("Total_Repeats: " + str(TotalRepeats))




#Plugging holes in the results:
for A in UniqueAList:
    for F in UniqueFList:
        if (A,F) not in CoordList:
            AreaFractionList.append(A)
            FitnessList.append(F)
            MeanDominationProb.append(0)



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

    TheoryAreaFractionList.append(
        1. -(1.-AreaFraction_LIMIT)**((2*Radius)/(2*Radius+a)))
    
    print("Fitness, a, AreaFraction: %0.3f, %0.3f, %0.3f"%
        (Fitness,a,TheoryAreaFractionList[-1]))






TheorySepList = sorted(TheorySepList)
print("THEORYSEPLIST",TheorySepList)
for s in range(len(TheorySepList)):
    TheorySepList[s] = TheorySepList[s]/Radius
print("THEORYSEPLIST",TheorySepList)
#TheoryFitnessList = LatticeTheory(TheorySepList,1.)
TheoryROWFitnessList = LatticeTheoryROW(TheorySepList,1.)
print("THEORYFITNESSLIST:",TheoryROWFitnessList)


TheoryRaceFitnessList = LatticeTheoryRACE(TheorySepList,1.)

TheoryAreaFractionList = []
for i in TheorySepList:
    TheoryAreaFractionList.append(2*np.pi * ((1./(2. + i))**2) / np.sqrt(3.))






#############################################################################
###Plotting##################################################################
#############################################################################


#######
TimeVar = False
#######
print("AreaFractionList:",AreaFractionList)
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
FitnessList = [x if not isinstance(x, np.ndarray) else x.item() for x in FitnessList]
print(FitnessList)
FitnessSet = set(FitnessList)

print("AreaSet:",AreaSet)
print("Fitness Set:",FitnessSet)

SortedAreaList = list(AreaSet)
SortedAreaList.sort()
dA = (SortedAreaList[-1]-SortedAreaList[-2])/2

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
#print("Number of times it was occupied:",num)
#print("Matrix total element number:",ProbMatrix.size)
#print(ProbMatrix)

#Add final elements (A[-1] +dA, F[-1]+dF)
SortedFitnessList.append(SortedFitnessList[-1]+SortedFitnessList[-1]-SortedFitnessList[-2])
SortedAreaList.append(SortedAreaList[-1]+SortedAreaList[-1]-SortedAreaList[-2])


#for i in range(1,len(SortedAreaList)):
    #print(SortedAreaList[i],SortedAreaList[i]-SortedAreaList[i-1])


fig=plt.figure(2)
ax = fig.add_subplot(111)

cp = ax.tricontourf(A.ravel(), F.ravel(), D.ravel(),10,cmap='coolwarm')

ax.pcolor(np.asarray(SortedAreaList)-dA,np.asarray(SortedFitnessList)-dF,np.array(ProbMatrix),cmap='coolwarm')

#ax.pcolor([AreaFractionList,FitnessList],np.array(MeanDominationProb),cmap='coolwarm')


#ax.scatter(
#    AreaFractionList,FitnessList,c=np.array(MeanDominationProb),s=150,
#    cmap='coolwarm')

#ax.plot(TheoryAreaFractionList,TheoryFitnessList,color='black')


ax.set_xticks([0.25,0.5,0.75])
ax.set_xticklabels(
    [r'$0.25$',r'$0.5$',r'$0.75$'],
    fontsize=7)

ax.set_yticks([0.7,0.8,0.9,1])
ax.set_yticklabels(
    ['$0.7$',r'$0.8$',r'$0.9$',r'$1.0$'],
    fontsize=7)




ax.set_ylabel(r'mutant fitness $F$',fontsize=8)
ax.set_xlabel(r'patch area ratio $\phi$',fontsize=8)

#plt.colorbar(cp,label='Probability of Domination',fontsize=15)

cbar = fig.colorbar(cp)
cbar.set_label(label='prob. of M dominating $P_D$',size=8)
cbar.ax.tick_params(labelsize=7)
ax.tick_params(axis='both', which='major', labelsize=20)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.figure(2).set_size_inches(9.00000/2.54, 6.000000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
plt.figure(2).axes[0].set_position([0.1, 0.1, 0.85, 0.85])#[0.18, 0.24, 0.8, 0.7])

#Adding coordinates of interest
plt.scatter([0.4],[0.9],s=50,color='white',zorder=10,marker='X',edgecolors='black')
plt.scatter([0.4],[0.98],s=50,color='white',zorder=10,edgecolors='black')
plt.scatter([0.8],[0.9],s=50,color='white',zorder=10,marker='P',edgecolors='black')

#plt.scatter([0.4],[0.9],s=30,color='#ff007f',zorder=10)
#plt.scatter([0.4],[0.98],s=30,color='#50C878',zorder=10)
#plt.scatter([0.8],[0.9],s=30,color='white',zorder=10)



plt.tight_layout()
#plt.grid(True)

plt.savefig(str(args.directory) + '/pcolor_Domination_NoTheory.png')

ax.plot(
    TheoryAreaFractionList,TheoryROWFitnessList,':k',linewidth=2,
    label="short term theory")


ax.plot(
    TheoryAreaFractionList,TheoryRaceFitnessList,'k',linewidth=2,
    label="long term theory")


upperlim = 0.907#np.pi/2 * Radius**2 / (np.sqrt(3)/4 * (2*Radius+1)**2)
ax.plot([upperlim,upperlim],[0.6,1.1],'--',color='black',linewidth=1)
ax.plot([-0.5,0.9],[1,1],'--',color='black',linewidth=1)

plt.ylim(min(SortedFitnessList)-dF,max(SortedFitnessList)-dF)
plt.xlim(min(SortedAreaList)-dA,max(SortedAreaList)-dA)


#plt.legend(loc='lower left',fontsize=7)

plt.savefig(str(args.directory) + '/pcolor_Domination.png')
plt.savefig(str(args.directory) + '/pcolor_Domination.eps')
plt.savefig(str(args.directory) + '/pcolor_Domination.pdf')

print("Finished")
