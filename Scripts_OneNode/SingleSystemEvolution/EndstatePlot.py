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

import time
starttime=time.time()
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
    ObsCenterList = data['ObsCenterList']
    ObstacleRadius = data['Radius']
    AreaFraction = data['AreaFraction']
    LatticeSurfSep = data['LatticeSurfSep']
    DataSizeList = data['DataSizeList']


#####
maxheight = 0
for i in range(len(OutputData)):
    if maxheight < int(OutputData[i][2]):
        maxheight = int(OutputData[i][2])+1


#Create EndState
print("Starting Endstate")
EndState = np.ones((maxheight,DomainWidth))*-1
for i in range(len(OutputData)):
    print(OutputData[i][1],OutputData[i][2])
    EndState[int(OutputData[i][2])][int(OutputData[i][1])] = OutputData[i][3]
print("Created Endstate")

#########################################################################
#########################################################################
#########################################################################
 #Create the bacground patches patches of the entire domain -
colnum = DomainWidth
rownum = maxheight
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



ObsCenterDict = defaultdict(lambda: defaultdict(list))
xsquaresnum = (np.ceil(float(DomainWidth)/ObstacleRadius))
ysquaresnum =(np.ceil(float(DomainHeight)/ObstacleRadius))

#xsquaresnum = 0

for O in ObsCenterList:
    x_square = int((O[0] - (O[0]%ObstacleRadius))/ObstacleRadius)
    y_square = int((O[1] - (O[1]%ObstacleRadius))/ObstacleRadius)
    ObsCenterDict[str(x_square)][str(y_square)].append([O[0],O[1]])

PlottedSites = 0

for i in x:
    #for n,j in enumerate(y):
    for j in y:
        if not(j%2):    #if on even rows
            xpos = i-dx/2.
        else:
            xpos = i
        ypos = j * np.sqrt(3.)/2.0

        x_square = int((i - (i%ObstacleRadius))/ObstacleRadius)
        y_square = int((j - (j%ObstacleRadius))/ObstacleRadius)


        fcol = '#CCCCCC'
        xdictsquare = -2
        ObstaclePlaced = False
        while xdictsquare < 2:
            #for xdictsquare in range(-1,2):
            xval = int((xdictsquare + x_square)%xsquaresnum)

            ydictsquare = -2
            while ydictsquare < 3:
                #for ydictsquare in range(-2,3):
                yval = int(ydictsquare + y_square)
                if (yval < ysquaresnum) and (yval >= 0):
                    for c in ObsCenterDict[str(xval)][str(yval)]:
                        centerx = c[0]
                        xdist = min(abs(centerx-xpos),
                            DomainWidth-abs(centerx-xpos))
                        if (((xdist)**2 + 0.75*(c[1] - j)**2) <
                                 ObstacleRadius**2):
                            fcol = '#666666'
                            xdictsquare += 3
                            ydictsquare += 5
                            ObstaclePlaced = True
                            break
                ydictsquare += 1
            xdictsquare += 1
        Opolygon = mpatches.RegularPolygon(
            (xpos, ypos), 6, hexsize*dx,facecolor=fcol)
        PlottedSites += 1
        if PlottedSites%DomainWidth == 0:
            print("Plotted Sites: " + str(PlottedSites))

        patches.append(Opolygon)

        
        addpoly = False

        if EndState[int(j)][int(i)] ==-1:
            if not ObstaclePlaced:
                polygon = mpatches.RegularPolygon((xpos, ypos), 6, hexsize*dx,facecolor='#CCCCCC',alpha=0.5)
                addpoly = True
        elif EndState[int(j)][int(i)] ==0:
            polygon = mpatches.RegularPolygon((xpos, ypos), 6, hexsize*dx,facecolor='#33CCFF',alpha=0.5)
            addpoly = True
        elif EndState[int(j)][int(i)] ==1:
            polygon = mpatches.RegularPolygon((xpos, ypos), 6, hexsize*dx,facecolor='#CC0033',alpha=0.5)
            addpoly = True

        if addpoly:
            patches.append(polygon)

print("Finished Background Obstacles")
background = PatchCollection(patches, match_original=True)

ax.add_collection(background)
ax.set_axis_off()

plt.xlim(-1,rownum)
plt.ylim(-1,colnum)
plt.axis('equal')

plt.title('\n'.join(wrap('MutProb: %0.5f, Fitness: %0.3f, AreaFraction: %0.3f'%(mutprob,fitness,AreaFraction),60)))

fig.savefig(str(args.directory) + "/ENDSTATE.png", dpi=300,bbox_inches="tight", pad_inches=0)

endtime = time.time()

print("Time taken:",endtime-starttime)
