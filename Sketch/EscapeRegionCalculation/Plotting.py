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
circle1 = plt.Circle((0, 0), 1, facecolor='w',edgecolor='k',linewidth=2)
ax.add_patch(circle1)

#The mutant Path
#Have it start at pi/4
plt.plot([np.cos(np.pi/4),0],[-np.sin(np.pi/4),2],zorder=5,color = red,linewidth=3)

#The WT Path
WTPathx = [1,1] 
WTPathy = [-np.sin(np.pi/4),0]

thetalist = np.linspace(0,np.pi/2-np.arccos(0.5),1000)

for theta in thetalist:
    WTPathx.append(np.cos(theta))
    WTPathy.append(np.sin(theta))

WTPathx.append(0)
WTPathy.append(2)

plt.plot(WTPathx,WTPathy,zorder=5,color=blue,linewidth=3)

"""
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
"""
##################
#Dotted Aid lines

OBx = [0,1]
OBy = [0,0]
plt.plot(OBx,OBy,'--k',linewidth=3)


OCx = [0,np.cos(np.pi/2-np.arccos(0.5))]
OCy = [0,np.sin(np.pi/2-np.arccos(0.5))]
plt.plot(OCx,OCy,'--k',linewidth=3)

ODx = [0,0]
ODy = [0,2]
plt.plot(ODx,ODy,'--k',linewidth=3)

OIx = [0,np.cos(np.pi/4)]
OIy = [0,-np.sin(np.pi/4)]
plt.plot(OIx,OIy,'--k',linewidth=3)

##################
#Adding letters
fs = 15

plt.text(-0.15,0.,"O",fontsize=fs,ha='center', va='center')

plt.text(-0.15,2,"D",fontsize=fs,ha='center', va='center')

plt.text(np.cos(np.pi/2-np.arccos(0.5))+0.1,np.sin(np.pi/2-np.arccos(0.5))+0.1,
    "C",fontsize=fs,ha='center', va='center')

plt.text(1.15,0,"B",fontsize=fs,ha='center', va='center')

plt.text(np.cos(np.pi/4),-np.sin(np.pi/4)-0.2,
    "I",fontsize=fs,ha='center', va='center')

plt.text(1,-np.sin(np.pi/4)-0.2,
    "A",fontsize=fs,ha='center', va='center')

plt.text(0.1,0.2,r"$\gamma$",fontsize=fs,ha='center', va='center')

plt.text(0.35,0.1,r"$\psi$",fontsize=fs,ha='center', va='center')

plt.text(0.3,-0.15,r"$\theta$",fontsize=fs,ha='center', va='center')
##################
#Formatting
plt.xlim(-1,1.2)
plt.ylim(-1,2.2)

plt.axis('off')
plt.axis('equal')
plt.tight_layout()


bbox = Bbox([[-1.1,-1.1],[1.3,2.3]]) #[x0,y0],[x1,y1]
bbox = bbox.transformed(ax.transData).transformed(fig.dpi_scale_trans.inverted())
#plt.figure().set_size_inches(1.500000/2.54, 4.500000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
#plt.figure().axes[0].set_position([0.05, 0.24, 0.90, 0.7])

plt.savefig("EscapeCalculation.pdf",bbox_inches=bbox)
