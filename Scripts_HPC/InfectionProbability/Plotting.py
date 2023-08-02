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
FMPList = []
InfectionProbList = []
MutationRateList = []


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
        DataList = data['DataList']


    meanprob = 0
    num = 0

    for i in range(len(DataList)):
        if DataList[i][1] == 1:
            meanprob += 1
            num += 1

        elif DataList[i][1] == 2:
            num += 1

    InfectionProbList.append(meanprob/num)
    FitnessList.append(float(fitness))
    MutationRateList.append(mutprob)
    

print(FitnessList)
fitnessset = set(FitnessList)

orderedfitness = list(sorted(fitnessset))


MutAndProbMatrix = []
for i in range(len(orderedfitness)):
    MutAndProbMatrix.append([[],[]])

for i in range(len(FitnessList)):
    for j in range(len(orderedfitness)):
        if orderedfitness[j] == FitnessList[i]:
            MutAndProbMatrix[j][0].append(MutationRateList[i])
            MutAndProbMatrix[j][1].append(InfectionProbList[i])

for i in range(len(MutAndProbMatrix)):
    (MutAndProbMatrix[i][0], 
    MutAndProbMatrix[i][1]) = zip(*sorted(zip(
        MutAndProbMatrix[i][0], MutAndProbMatrix[i][1])))



#Theories
ZerothOrder = []    #Surface Mutations
FirstOrder = []     #Mean cluster sizes below, and upper surface mutations
SecondOrder = []    #Flat Front cluster height distribution employed


def probless(fitness,d):
    """Probability that the size of a cluster of M fitness F is less
    than d"""
    prob = 0

    for i in range(1,d+1):
        prob += ((((fitness)**(i-1))/((1+fitness)**(2*i))) * 
            (1/(i+1))*scipy.special.binom(2*i,i))
    return prob

for i in FitnessList:
    ZerothOrder.append(1-(1-mutprob)**(2*np.pi*Radius))
    FirstOrder.append(1-((1-mutprob)**(np.pi*Radius*np.ceil(
        (1+i)/(1-i)))) * (1-mutprob)**(np.pi*Radius))

    probnot = 1
    for h in range(2,100):   #distance from the patch
        prob_less = probless(i,h)
        print(h,prob_less)
        probnot *= ((1 - mutprob) + mutprob*prob_less)**(np.pi*Radius)

    print(probnot)
    SecondOrder.append(1 - probnot*(1-mutprob)**(2*np.pi*Radius))




#############################################################################
#############################################################################
#############################################################################

#Sight line: for highlighting the linear relation
guidexlist = np.linspace(1e-6,1e-4,200)
guideylist = guidexlist*1e4


plt.figure(1)
for i in range(len(orderedfitness)):
    plt.loglog(
        (MutAndProbMatrix[i][0]),
        (MutAndProbMatrix[i][1]),
        label='%0.1f'%(orderedfitness[i]))

plt.loglog(guidexlist,guideylist,'k',linewidth=1)

plt.annotate(r"$\sim\mu$",[1e-6,3e-1],fontsize=8)


plt.legend(loc='lower right',title=r'$F$',fontsize=7,frameon=False)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.xlabel("mutation prob. "+r"$\mu$", fontsize=8)
plt.ylabel("prob. patch invasion", fontsize=8)

plt.figure(1).axes[0].set_xticks([1.e-6, 1e-5, 1.e-4, 1.e-3,1.e-2])
plt.figure(1).axes[0].set_yticks([1e-04, 1e-3, 1e-2, 1e-1, 1])

plt.figure(1).set_size_inches(4.500000/2.54, 4.500000/2.54, forward=True)
plt.figure(1).axes[0].set_position([0.35, 0.25, 0.6, 0.7])

plt.minorticks_off()


"""
plt.grid()
plt.legend(loc='lower right',fontsize=20)

plt.xlabel("Mutation Probability")
plt.ylabel("Probability of Patch Invasion")

plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
"""
"""
yint = []
locs, labels = plt.yticks()
for each in locs:
    yint.append(int(each))
plt.yticks(yint)

xint = []
locs, labels = plt.xticks()
for each in locs:
    xint.append(int(each))
plt.xticks(xint)
"""
plt.savefig(str(args.directory) + '/InfectionProb.png')
plt.savefig(str(args.directory) + '/InfectionProb.pdf')
