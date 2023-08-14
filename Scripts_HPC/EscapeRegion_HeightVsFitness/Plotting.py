import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
#import pylustrator

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
from scipy.stats import gaussian_kde

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

CapHeightList = []
FitnessList = []

MeanFitnessList = []
MeanCapHeightList = []
MedianCapHeightList = []

UniqueAngleList = []

TheoryMeanCapHeight = []
TheoryRadiusList = []
TheoryFitnessList = []

for d in dirlist:
    filelist = os.listdir(args.directory + "/" + d)
    for names in filelist:
        if names.endswith(".npz"):
            filename = names
            print("Found File!")

    with np.load(os.path.join(args.directory,d,filename),allow_pickle=True) as data:
        Repeats = data['Repeats']
        fitness = data['fitness']
        mutprob = data['mutprob']
        Radius = data['Radius']
        SeededPos = data['SeededPos']
        DataList = data['DataList'] #ObsCenterList,MaxMutantHeight

    print("Radius %0.3f, Fitness %0.6f"%(Radius,fitness))


    TheoryRadiusList.append(Radius)
    TheoryFitnessList.append(fitness)
    TheoryMeanCapHeight.append(
        TheoreticalPredictionFuns.WaveObstacleCapHeight(
            -np.pi/2,fitness,Radius)/Radius)

    mean = 0
    medianlist = []

    HistogramList = []

    angle = -np.pi/2
    """
    ydisp = (np.sqrt(3.)/2.) *(SeededPos[1] - DataList[0][0][0][1])
    xdisp = (SeededPos[0] - DataList[0][0][0][0])

    
    if ydisp == 0:
        if xdisp > 0:
            angle = np.pi/2
        else:
            angle = -np.pi/2
    elif ydisp < 0:
        if xdisp > 0:
            angle = abs(np.arctan(xdisp/ydisp))
        else:
                        angle = -np.arctan(xdisp/ydisp)
    else:
        if xdisp > 0:
            angle = np.arctan(ydisp/xdisp) + np.pi/2
        else:
            angle = -np.pi/2 - abs(np.arctan(ydisp/xdisp))

    if angle not in UniqueAngleList:
        UniqueAngleList.append(angle)
    """
    meanlist = []
    for R in range(Repeats):
        print("Obs Center:",DataList[R][0])
        print("Max height:",DataList[R][1])

        capheight = ((DataList[R][1] - DataList[R][0][0][1]
             - int(Radius*2/np.sqrt(3)))*(np.sqrt(3)/2.))
        capratio = capheight/Radius
        print("HEIGHT:",capheight)
        if capheight >= 0:
            print("HEIGHT:",capheight)
            FitnessList.append(fitness)
            CapHeightList.append(capratio)
            HistogramList.append(capheight)
            meanlist.append(capratio)

    MeanFitnessList.append(fitness)
    MeanCapHeightList.append(np.mean(meanlist))
    MedianCapHeightList.append(np.median(meanlist))

    print(angle)
    print(HistogramList)
    """
    #Histogram of cap heights
    plt.figure(1)
    upperlim = 300
    binstep = 10
    bins = np.arange(0,upperlim+binstep,binstep)
    

    plt.hist(HistogramList,bins=bins)

    plt.axvline(np.mean(HistogramList),color='r',linewidth=3,label='Mean')
    plt.axvline(np.median(HistogramList),color='k',linewidth=3,label='Median')


    effectiveangle=0
    if angle > 0:
        #Function takes Angle as measured from the right hand side
        effectiveangle = angle - np.pi/2
    else:
        effectiveangle = abs(angle) - np.pi/2

    plt.axvline(
        TheoreticalPredictionFuns.WaveObstacleCapHeight(
            effectiveangle,fitness,Radius),color='orange',linewidth=3,
            label='Theory')

    plt.xlim(0,300)

    plt.xlabel("Cap Height")
    plt.ylabel("Frequency")

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    #plt.legend(loc='upper right',fontsize=30)

    plt.grid()

    plt.savefig(args.directory + '/Hist_fitness_%0.3f.png'%(fitness))
    plt.close()
    """


#Sort the mean:
(MeanFitnessList,
MeanCapHeightList,
MedianCapHeightList) = zip(*sorted(zip(MeanFitnessList,
     MeanCapHeightList, MedianCapHeightList)))
print("Length of capheightlist",len(CapHeightList))

TheoryFitnessList = sorted(FitnessList)
TheoryCapHeight = []
for fitness in TheoryFitnessList:
    effectiveangle = 0
    TheoryCapHeight.append(
        TheoreticalPredictionFuns.WaveObstacleCapHeight(
            -np.pi/2,fitness,Radius)/Radius)


#Colour map based on density
FitnessList = np.array(FitnessList)
CapHeightList = np.log(np.array(CapHeightList))

#pylustrator.start()

fig = plt.figure(2)
ax = fig.add_subplot(1, 1, 1)

plt.plot(MeanFitnessList,MedianCapHeightList,'k',label="Median")
plt.plot(TheoryFitnessList,TheoryCapHeight,'--c',label="Theory")

#plt.legend(loc='lower center',fontsize = 20)
#plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0,fontsize=20)

plt.xlim(0.7,0.95)
plt.ylim(-0.1,2)
#plt.ylim(bottom = 0)

ax.set_xticks([0.7,0.8,0.9,1.0])
ax.set_xticklabels(
    [r'$0.7$',r'$0.8$',r'$0.9$',r'$1.0$'],
    fontsize=7)

plt.yticks(fontsize=7)

#plt.ylim(0,np.ceil(max(MeanCapHeightList)))

plt.xlabel(r'fitness $F$',fontsize=8)
plt.ylabel(r'escape reg. height $\Lambda$',fontsize=8)

plt.figure(2).set_size_inches(4.000000/2.54, 4.500000/2.54, forward=True)
plt.figure(2).axes[0].set_position([0.45, 0.24, 0.5, 0.7])#[0.18, 0.24, 0.8, 0.7])

#plt.tight_layout()
plt.savefig(str(args.directory) + '/Median_CapHeight.pdf')
plt.savefig(str(args.directory) + '/Median_CapHeight.eps')
#% start: automatic generated code from pylustrator
#% end: automatic generated code from pylustrator
#plt.show()

#############################################################################
############################################################################
#############################################################################
