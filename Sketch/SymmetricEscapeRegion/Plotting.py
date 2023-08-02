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

from matplotlib.transforms import Bbox


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

blue = '#33CCFF'
red = '#CC0033'

####################################
fig, ax = plt.subplots()

#The Cirlce
circle1 = plt.Circle((0, 0), 1, facecolor='w',edgecolor='k',linewidth=3)
ax.add_patch(circle1)

#The mutant Path
plt.plot([0,0],[-1,2],color = red,linewidth=3)

#The WT Path
WTPathx = [1,1] 
WTPathy = [-1,0]

thetalist = np.linspace(0,np.pi/2-np.arccos(0.5),1000)

for theta in thetalist:
    WTPathx.append(np.cos(theta))
    WTPathy.append(np.sin(theta))

WTPathx.append(0)
WTPathy.append(2)

plt.plot(WTPathx,WTPathy,color=blue,linewidth=3)

#Arrows
#Mutant Arrows
plt.arrow(0,-0.5,0,0.1,head_width=0.1,color=red)
plt.arrow(0,1,0,0.1,head_width=0.1,color=red)

#WT Arrows
plt.arrow(1,-0.5,0,0.1,head_width=0.1,color=blue)

slope = -(2-np.sin(np.pi/2-np.arccos(0.5)))/(
        np.cos(np.pi/2-np.arccos(0.5)))

dx = 0.1/np.sqrt(1+slope**2)
dy = dx * slope

y = 1
x = -1/slope

plt.arrow(x,y,-dx,-dy,head_width=0.1,color=blue)


#Formatting
plt.xlim(-1,1)
plt.ylim(-1,2)

plt.axis('off')
plt.axis('equal')
plt.tight_layout()


bbox = Bbox([[-1.1,-1.1],[1.1,2.1]]) #[x0,y0],[x1,y1]
bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())
#plt.figure().set_size_inches(1.500000/2.54, 4.500000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
#plt.figure().axes[0].set_position([0.05, 0.24, 0.90, 0.7])

plt.savefig("SymmetricPatch.pdf",bbox_inches=bbox)
