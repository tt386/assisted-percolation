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
r = 1
s = r/4

mid = r + s/2

#Left
leftcircle = plt.Circle((0, 0), r, facecolor='w',edgecolor='k',linewidth=2)
ax.add_patch(leftcircle)

#right
rightcircle = plt.Circle((2*r+s,0),r, facecolor='w',edgecolor='k',
                        linewidth=2)
ax.add_patch(rightcircle)

#top middle
topx = (2*r+s)*np.cos(np.pi/3)
topy = (2*r+s)*np.sin(np.pi/3)
topcircle = plt.Circle((topx,topy),r, facecolor='w',edgecolor='k',
                        linewidth=2)
ax.add_patch(topcircle)

#above
highercircle = plt.Circle((0,2*topy),r,facecolor='w',edgecolor='k',
                        linewidth=2)
ax.add_patch(highercircle)



#The Mutant path
plt.plot([mid,mid],[0,topy-r],color=blue,linewidth=3)
plt.plot([mid,mid],[topy-r,topy+2*r],color=red,linewidth=3)

#M path for escape angle
plt.plot([mid,topx-np.cos(np.pi/2-np.arccos(0.5))],
        [topy-r,topy+np.sin(np.pi/2-np.arccos(0.5))],
        color=red,linewidth=3)


#M path for the escape envelope
plt.plot([topx-np.cos(np.pi/2-np.arccos(0.5)),topx],
        [topy+np.sin(np.pi/2-np.arccos(0.5)),topy+2*r],
        '--',color=red,linewidth=3)

#The Wt path
tangentangle = np.pi/3 - np.arccos(r/(r+s/2))
print(tangentangle)

wtx = []
wty = []

for i in np.linspace(0,tangentangle,1000):
    wtx.append(r * np.cos(i))
    wty.append(r * np.sin(i))

wtx.append(topx - r*np.cos(tangentangle))
wty.append(topy - r*np.sin(tangentangle))


for i in np.linspace(-tangentangle,np.pi/2-np.arccos(0.5),1000):
    wtx.append(topx - r*np.cos(i))
    wty.append(topy + r*np.sin(i))

"""
wtx.append(topx)
wty.append(topy+2*r)
"""
plt.plot(wtx,wty,color=blue,linewidth=3)


#M path for escape angle




#####Guide lines
#Distance between adjacent edges
plt.plot([r,r+s],[0,0],'--',color='k',linewidth=3)

#Distance between adjacent centers O-O'
plt.plot([2*r+s,topx],[0,topy],'--',color='k',linewidth=3)


#Angle DOE
#plt.plot([0.5*(2*r+s)-r,0.5*(2*r+s),0.5*(2*r+s)])

#Guide points
m = 30
fs = 15 


plt.scatter([0],[0],color='k',s=m,zorder=3)
plt.text(0.3,0,
    "O''",fontsize=fs,ha='center', va='center')

plt.scatter([2*r+s],[0],color='k',s=m,zorder=3)
plt.text(2*r+s+0.2,0,
    "O'",fontsize=fs,ha='center', va='center')

plt.scatter([topx],[topy],color='k',s=m,zorder=3)
plt.text(topx+0.2,topy,
    "O",fontsize=fs,ha='center', va='center')

plt.scatter([topx],[topy+2*r],color='k',s=m,zorder=3)
plt.text(topx+0.2,topy+2*r,
    "F",fontsize=fs,ha='center', va='center')

plt.scatter([r],[0],color='k',s=m,zorder=3)
plt.text(r-0.2,0,
    "A",fontsize=fs,ha='center', va='center')

plt.scatter([r+s/2],[0],color='k',s=m,zorder=3)
plt.text(r+s/2,0-0.3,
    "G",fontsize=fs,ha='center', va='center')

plt.scatter([r*np.cos(tangentangle)],
            [r*np.sin(tangentangle)],color='k',s=m,zorder=3)
plt.text(r*np.cos(tangentangle)-0.2,r*np.sin(tangentangle),
    "B",fontsize=fs,ha='center', va='center')

plt.scatter([r+s/2],[topy-r],color='k',s=m,zorder=3)
plt.text(r+s/2,topy-r-0.3,
    "H",fontsize=fs,ha='center', va='center')

plt.scatter([topx - r*np.cos(tangentangle)],
            [topy - r*np.sin(tangentangle)],color='k',s=m,zorder=3)
plt.text(topx - r*np.cos(tangentangle)-0.2,topy - r*np.sin(tangentangle),
    "C",fontsize=fs,ha='center', va='center')

plt.scatter([topx-r],[topy],color='k',s=m,zorder=3)
plt.text(topx-r-0.2,topy,
    "D",fontsize=fs,ha='center', va='center')

plt.scatter([topx-r*np.cos(np.pi/2-np.arccos(0.5))],
            [topy + r*np.sin(np.pi/2-np.arccos(0.5))],color='k',s=m,zorder=3)
plt.text(topx-r*np.cos(np.pi/2-np.arccos(0.5))-0.2,
        topy + r*np.sin(np.pi/2-np.arccos(0.5)),
        "E",fontsize=fs,ha='center', va='center')

#Angle DOE
plt.plot([topx-r/2,topx,topx-(r*np.cos(np.pi/2-np.arccos(0.5)))],
        [topy,topy,topy + (r*np.sin(np.pi/2-np.arccos(0.5)))],
        color='k',zorder=1)
plt.text(topx-0.5,topy+0.1,r"$\beta$",fontsize=fs,ha='center',va='center')

#Formatting
plt.xlim(-1,1)
plt.ylim(-1,2)

plt.axis('off')
plt.axis('equal')
plt.tight_layout()


bbox = Bbox([[-1.1,-1.1],[3.1,4.2]]) #[x0,y0],[x1,y1]
bbox = bbox.transformed(ax.transData).transformed(
        fig.dpi_scale_trans.inverted())
#plt.figure().set_size_inches(1.500000/2.54, 4.500000/2.54, forward=True)#8.200000/2.54, 4.500000/2.54, forward=True)
#plt.figure().axes[0].set_position([0.05, 0.24, 0.90, 0.7])

plt.savefig("ShortTerm.pdf",bbox_inches=bbox)
