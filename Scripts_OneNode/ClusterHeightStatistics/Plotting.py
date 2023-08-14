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

from collections import Counter

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
#mpl.rcParams['agg.path.chunksize'] = 10000

from textwrap import wrap

import copy
from copy import deepcopy

import matplotlib.ticker as mtick

from matplotlib.ticker import MaxNLocator

def ticks(y, pos):
    return r'$e^{:.0f}$'.format(np.log(y))

from scipy import stats
from scipy.stats import gaussian_kde
import scipy
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


#Find all the npz files within the directory
templist = os.listdir(args.directory)
print(templist)

npzlist = []

for i in range(len(templist)):
    if templist[i].endswith(".npz"):
        npzlist.append(templist[i])


print(npzlist)


#############################################################################
###Data extraction###########################################################
#############################################################################

RoughClusterHeightList = []
FlatClusterHeightList = []
FlatClusterWidthList = []
F = 0

WidthDirectoryName = args.directory + '/WidthGivenHeightPlots'
if not os.path.isdir(WidthDirectoryName):
    os.mkdir(WidthDirectoryName)

for n in npzlist:
    with np.load(os.path.join(args.directory,n)) as data:

        if n.startswith("Rough"):
            Output = data['DataList']
            F = data['fitness']
            for i in range(len(Output)):
                RoughClusterHeightList.append(
                    Output[i][0] - Output[i][1] + 1)
        
        elif  n.startswith("Flat"):
            FlatClusterHeightList = data['HeightList']
            F = data['F']
            FlatClusterWidthList = data['WidthList']

#############################################################################
###Height list shenanigans###################################################
#############################################################################
#pylustrator.start()

plt.figure(1)

#Check for if there was a rough directory involved
if len(RoughClusterHeightList) > 0:
    Roughxlist = np.arange(max(RoughClusterHeightList))
    Roughfreqlist = np.zeros(len(Roughxlist))

    for i in range(len(RoughClusterHeightList)):
        Roughfreqlist[RoughClusterHeightList[i]-1] += 1
    Roughfreqlist /= sum(Roughfreqlist)

    plt.loglog(Roughxlist,Roughfreqlist,'k',label='Sim: Rough')


#Check for if there was a flat directory involved
if len(FlatClusterHeightList) > 0:
    Flatxlist = np.arange(max(FlatClusterHeightList))
    Flatfreqlist = np.zeros(len(Flatxlist))
    for i in range(len(FlatClusterHeightList)):
        Flatfreqlist[FlatClusterHeightList[i]-1] += 1
    Flatfreqlist /= sum(Flatfreqlist)



    Flattheorylist = []
    for i in Flatxlist:
        H = i+1
        Flattheorylist.append(
            ((F**(H-1))/((1+F)**(2*H))) * (1/(H+1)) * 
            scipy.special.binom(2*H,H))

    plt.loglog(Flatxlist,Flatfreqlist,'orange',label='Sim: Flat')
    plt.loglog(Flatxlist,Flattheorylist,'--b',label='Theory: Flat')

plt.legend(loc='upper right',fontsize=5,frameon=False)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.xlabel(r"cluster height $H$", fontsize=8)
plt.ylabel("probability density", fontsize=8)

plt.figure(1).axes[0].set_xticks([1.0, 10.0, 100.0, 1000.0])
plt.figure(1).axes[0].set_yticks([1e-05, 0.0001, 0.001, 0.01, 0.1])

plt.figure(1).set_size_inches(4.500000/2.54, 4.500000/2.54, forward=True)
plt.figure(1).axes[0].set_position([0.35, 0.25, 0.6, 0.7])

#plt.legend(loc="upper right",fontsize=8, frameon=False)


plt.savefig(args.directory + '/Rough_Vs_Flat.pdf')
plt.savefig(args.directory + '/Rough_Vs_Flat.png')
plt.show()

#############################################################################
###Width List Shenanigans####################################################
#############################################################################

zipped_lists = zip(FlatClusterHeightList, FlatClusterWidthList)
sorted_pairs = sorted(zipped_lists)

tuples = zip(*sorted_pairs)
HeightList, WidthList = [ list(tup) for tup in  tuples]

WidthGivenHeightList = [[HeightList[0]]]
for i in range(1,len(HeightList)):
    if HeightList[i] != HeightList[i-1]:
        WidthGivenHeightList.append([])

    WidthGivenHeightList[-1].append(WidthList[i])



for H in range(1,41):
    print(H)
    TheoryWidthProbList = TheoreticalPredictionFuns.ClusterWidthDist(F,H)
    
    plt.figure(2)
    bins = np.arange(max(WidthGivenHeightList[H-1])+1) - 0.5

    countdict = Counter(WidthGivenHeightList[H-1])

    countdictkeys = countdict.keys()

    countdictvals = []

    for i in countdictkeys:
        countdictvals.append(countdict[i])

    countdictvals = np.asarray(countdictvals)/sum(countdictvals)

    #plt.hist(WidthGivenHeightList[H-1],bins,density=True)
    plt.scatter(countdictkeys,countdictvals,s=100,marker="x",color='k',zorder=5)
    plt.plot(np.arange(len(TheoryWidthProbList))+1,TheoryWidthProbList,'.c-',
        label='Analytical',linewidth=3,markersize=20)

    #plt.legend(loc='upper right',fontsize=30,frameon=False)

    plt.xticks(fontsize=30)

    xint = []
    locs, labels = plt.xticks()
    for each in locs:
        xint.append(int(each))
    plt.xticks(xint)
    
    plt.yticks(fontsize=30)
   
    #plt.figure(2).axes[0].set_xticks(list(np.arange(0,max(WidthGivenHeightList[H-1])+1)))

    plt.figure(2).axes[0].set_yticks([0, 0.2,0.4,0.6])
    plt.figure(2).axes[0].set_xticks([0,2,4,6,8,10,12])

 
    plt.xlabel("cluster width",fontsize=30)
    plt.ylabel("probability",fontsize=30)

    #plt.figure(1).set_size_inches(4.500000/2.54, 4.500000/2.54, forward=True)
    plt.figure(2).axes[0].set_position([0.20, 0.25, 0.75, 0.65])

    plt.title("cluster height: %d"%(H),fontsize=30)
    plt.xlim(0,12)
    plt.ylim(0,0.6)
    plt.savefig(WidthDirectoryName + "/Height_%0.4d.png"%(H))
    plt.savefig(WidthDirectoryName + "/Height_%0.4d.eps"%(H))
    plt.close()



