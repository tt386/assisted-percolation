"""
The argument required after this script name is the directory of the stored data.
"""
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
from matplotlib import colors
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

from matplotlib import animation

from copy import copy
from copy import deepcopy

#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
#mpl.rcParams['agg.path.chunksize'] = 10000

from textwrap import wrap

from scipy import stats

import copy

from collections import defaultdict

################################
##ArgParse######################
################################
import os.path

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The Directory %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

from argparse import ArgumentParser

parser = ArgumentParser(description='FieldPlot')
parser.add_argument('-d','--directory',help='The directory of the data')
args = parser.parse_args()

#############################################################################
#############################################################################
###NPZ FILES
filename = 'timecourse.npz'
with np.load(os.path.join(args.directory, filename)) as data:
    OutputData = data['OutputData']
    DomainHeight = data['DomainHeight']
    DomainWidth = data['DomainWidth']
    fitness = data['fitness']
    mutprob = data['mutprob']
    DataSizeList = data['DataSizeList']


#Finished dat stuff
#############################################################################
#############################################################################
#Datasize plotting
#There was a time I cared immensly about the data usage, so I plotted it so
#I could compare various techniques

simsteplist = []
memorylist = []

for i in DataSizeList:
    simsteplist.append(i[0])
    memorylist.append(i[1])

plt.figure(1)
plt.plot(simsteplist,memorylist)
plt.title("Memory Consumed as SimStep increases.")
plt.xlabel("SimStep")
plt.ylabel("Memory Consumed")
plt.grid()
plt.savefig(args.directory + '/MemoryPlot.png')

#############################################################################
#############################################################################
#############################################################################
















#Construct the CellList (a matrix of sites we colour to represent the system)
MaxHeight = 0
for O in OutputData:
    if O[2] > MaxHeight:
        MaxHeight = int(copy.deepcopy(O[2]))




OneSiteAtATime = False


#Construct the world, empty so far:
CellList = np.zeros(shape=(MaxHeight,DomainWidth),dtype=int)
print("Length of CellList: " + str(len(CellList)))


#Let's have 5 second differences ('second' is subjective) betwwen frames
tdiff = 5
maxtime = OutputData[-1][0]


#Cut the OutputData into slices: Create list of indices which represent time
#bounds so each frame is advanced by that amount

TimeBoundsList = []
TimeBoundsList.append(0)
#Start for t = 0
counter = 0

while OutputData[counter][0] == 0:
    counter += 1

TimeBoundsList.append(counter)

#For each 5-second gap:
ledgecounter = 1

while counter < len(OutputData):
    if OutputData[counter][0] > ledgecounter*tdiff:
        ledgecounter += 1
        TimeBoundsList.append(counter)
    counter += 1

#Add last element
TimeBoundsList.append(counter)

if OneSiteAtATime:
    TimeBoundsList = []
    for i in range(len(OutputData)):
        TimeBoundsList.append(colnum + i)






#############################################################################
#############################################################################
#############################################################################
SYSTEMPLOT = True
rownum = len(CellList)
colnum = len(CellList[0])
if SYSTEMPLOT:
    #########################################################################
    #########################################################################
    #########################################################################
    #Create the bacground patches patches of the entire domain -
    fig, ax = plt.subplots()

    #background
    patches = []

    x = np.linspace(0,colnum-1,colnum)
    y = np.linspace(0,rownum-1,rownum)
    X,Y = np.meshgrid(x,y)

    dx = np.diff(x)[0]
    dy = np.diff(y)[0]
    ds = np.sqrt(dx**2 + dy**2)

    hexsize = 0.55


    PlottedSites = 0

    for i in x:
        #for n,j in enumerate(y):
        for j in y:
            if not(j%2):    #if on even rows
                xpos = i-dx/2.
            else:
                xpos = i
            ypos = j * np.sqrt(3.)/2.0

    

            fcol = '#CCCCCC'
            polygon = mpatches.RegularPolygon(
                (xpos, ypos), 6, hexsize*dx,facecolor=fcol)
            PlottedSites += 1
            if PlottedSites%DomainWidth == 0:
                print("Plotted Sites: " + str(PlottedSites))
            

            patches.append(polygon)
    print("Finished Background Obstacles")
    background = PatchCollection(patches, match_original=True)

    ax.add_collection(background)
    ax.set_axis_off()

    plt.xlim(-1,rownum)
    plt.ylim(-1,colnum)
    plt.axis('equal')
    
    fig.savefig(str(args.directory) + "/Background.png")
    plt.title('\n'.join(wrap(
        'Mutation Prob: %0.3f, Mutant Fitness: %0.3f'%(mutprob,fitness),60)))

    #########################################################################
    #########################################################################
    #########################################################################
    #Now, iterate through the TimeBoundsList, plotting an image at each point
    
    for timebound in range(len(TimeBoundsList) - 1):

        newpatches = []
        #print(1)
        CurrentLower = TimeBoundsList[timebound]
        CurrentUpper = TimeBoundsList[timebound + 1]
        for timevals in range(CurrentLower, CurrentUpper):
            col = int(OutputData[timevals][1])
            row = int(OutputData[timevals][2])
            #print(2)
            if not(row%2):      #if on even rows
                xpos = col-dx/2.
            else:               
                xpos = col      #if on odd rows

            ypos = row * np.sqrt(3.)/2.0

            if OutputData[timevals][3] == 0:
                polygon = mpatches.RegularPolygon(
                    (xpos, ypos), 6, hexsize*dx,facecolor='#33CCFF',
                    alpha=0.5)
            elif OutputData[timevals][3] == 1:
                polygon = mpatches.RegularPolygon(
                    (xpos, ypos), 6, hexsize*dx,facecolor='#CC0033',
                    alpha=0.5)
            newpatches.append(polygon)
    

        newcollection = PatchCollection(newpatches, match_original=True)

        ax.add_collection(newcollection)
        
        plt.xlim(-1,rownum)
        plt.ylim(-1,colnum)
        plt.axis('equal')

        Savename = 'AField_%0.5d'%(timebound) + '.png'

        plt.savefig(
            str(args.directory) + "/" + Savename,dpi=300,bbox_inches="tight",
            pad_inches=0)

        print("Produced plot: "+str(timebound)+'/'+str(len(TimeBoundsList)))

    plt.close()

