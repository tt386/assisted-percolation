"""
Includes all the predictions from theory
"""
import numpy as np
import scipy

import copy
from copy import deepcopy

#############################################################################
#############################################################################
#############################################################################
def ClusterHeightDist(Fitness,MaxClusterHeight):
    """
    Calculate the height distribution of a cluster of mutants of Fitness,
    on a flat front, up to a cutoff heigth MaxClusterHeight, and output
     this as a list. Uses a matrix approach.

    ARGS:
    Fitness: float
        The fitness of the M sites
    MaxClusterHeight: int
        The largest clustr height included in the height distribution

    RETURNS:
    List
        The height distribution of clusters

    """


    a = (Fitness/(1.+Fitness))**2
    b = 2.*Fitness/((1.+Fitness)**2)
    c = 1./((1.+Fitness)**2)

    TransitionMatrix = np.zeros(
        (MaxClusterHeight,MaxClusterHeight),
        dtype=float)


    #Populate the transition matrix:
    TransitionMatrix[0,0] = 1   #Remember [row,col]


    for i in range(1,MaxClusterHeight):
        TransitionMatrix[i,i-1] = c
        TransitionMatrix[i,i] = b
        if i < MaxClusterHeight-1:
            TransitionMatrix[i,i+1] = a

    
    #Multiply Calculate the (row = 1,column=0) site for T^n-T^(n-1) 
    #to find probability of ending at that time given 
    
    ProbList = np.zeros((MaxClusterHeight),dtype=float)
    ProbList[0] = c

    TM1 = deepcopy(TransitionMatrix)
    TM2 = np.matmul(TransitionMatrix,TransitionMatrix)

    for i in range(1,MaxClusterHeight):
        print(i)
        ProbList[i] = (TM2-TM1).item(MaxClusterHeight)
        TMTemp = deepcopy(TM2)

        TM2 = np.matmul(TM2,TransitionMatrix)
        TM1 = TMTemp


    return ProbList


#############################################################################
#############################################################################
#############################################################################
def ClusterWidthDist(Fitness,ClusterHeight):
    """
    Calculate the width distribution of a flat front cluster of Fitness,
    for a flat front cluster of height ClusterHeight

    ARGS:
    Fitness: float
        The fitness of the M sites
    ClusterHeight: int
        The Height of the flat front cluster

    RETURNS:
    List
        The width distribution of the clusters

    """


    def Tstar(n,k):
        """Used to find the number of ways a dyck path length 2n has 
        height k
        """
        tempsum = 0.0
        for i in range(1,k+2):
            tempsum += (np.sin(np.pi*i/(k+2))*(np.cos(np.pi*i/(k+2)))**n)**2

        return ((2**(2*n+1))/(k+2))*tempsum

    def T(n,k):
        """The number of ways a dyck path length 2n has height k"""
        return np.round(Tstar(n,k) - Tstar(n,k-1))



    def Cn(n):                              
        """The nth catalan number"""
        return (1/(n+1)) * scipy.special.binom(2*n,n)


    ProbList = np.zeros(int(ClusterHeight/2),dtype=float)


    #Iterate through possible widths
    for W in range(1,int(np.ceil(float(ClusterHeight)/2.0))):
        ProbWidthAndHeightSum = 0

        for x in range(W-1,int((ClusterHeight-1)/(2)) + 1):
            ProbWidthAndHeightSum += (2**(ClusterHeight-1)/
                 Cn(ClusterHeight)* T(x,W-1) * 
                 scipy.special.binom(ClusterHeight-1,2*x) * 1/(2**(2*x)))

        ProbList[W-1] = ProbWidthAndHeightSum#/ProbClusterHeight
    

    return ProbList




#############################################################################
#############################################################################
#############################################################################
def WaveObstacleCapHeight(InfectionAngle,Fitness,ObstacleRadius):
    """
    Use the principle of least time to calculate the height of the emergence
    region of the M above the patch, for M of a given Fitness, invading
    a patch of radius ObstacleRadius at an angle InfectionAngle


    ARGS:
    InfectionAngle: float
        The angle, as measured in polar coordinates, at which the patch
         becomes invaded
    Fitness: float
        The fitness of the M sites
    ObstacleRadius: int
        The radius of the patch

    RETURNS:
    float:
        Height of the escape region

    """

    from scipy import optimize

    def function(x,InfectionAngle,Fitness,Radius):
        if InfectionAngle > 0:
            return ((1.0/Fitness) * np.sqrt((Radius**2 ) + (x**2) - 
                2*x*Radius*np.sin(InfectionAngle)) + Radius * 
                InfectionAngle - Radius * np.arcsin(Radius/x) - 
                np.sqrt((x)**2 - (Radius)**2))

        else:
            return ((1.0/Fitness) * np.sqrt((Radius**2 ) + (x**2) + 
                2*x*Radius*np.sin(-InfectionAngle)) - Radius * 
                np.sin(-InfectionAngle) - Radius * np.arcsin(Radius/x) - 
                np.sqrt((x)**2 - (Radius)**2))

    if Fitness >= 1:
        return np.inf

    #Equations only valid for -np.pi/2 < angle < np.pi/2
    if InfectionAngle > np.pi/2:
        InfectionAngle = np.pi-InfectionAngle
    elif InfectionAngle < -np.pi/2:
        InfectionAngle = -np.pi - InfectionAngle

    #Newtwo initial guess
    InitialGuess = ObstacleRadius

    try:
        root = (optimize.newton(
            function,
            InitialGuess,
            args=(InfectionAngle,Fitness,ObstacleRadius)) - ObstacleRadius)

    except:
        root = 0

    return root


#############################################################################
#############################################################################
#############################################################################
def WaveMutantWTAngleMeet(InfectionAngle,Fitness,ObstacleRadius): 
    """
    Use the principle of least time to calculate the angle around the patch
    perimeter at which the M and WT strains would reach that point at the
    same time, for M of a given Fitness, invading a patch of radius 
    ObstacleRadius at an angle InfectionAngle


    ARGS:
    InfectionAngle: float
        The angle, as measured in polar coordinates, at which the patch
         becomes invaded
    Fitness: float
        The fitness of the M sites
    ObstacleRadius: int
        The radius of the patch

    RETURNS:
    float:
        Angle at which WT and M would meat around the patch perimeter
    """


    from scipy import optimize
    
    def function(x,InfectionAngle,Fitness,Radius):
        if InfectionAngle > 0:
            return ((1.0/Fitness) * np.sqrt(2.0 - 2.0 * 
                np.cos(x - InfectionAngle)) - x + InfectionAngle)
        else:
            return ((1.0/Fitness) * np.sqrt(2.0 - 2.0 * 
                np.cos(x - InfectionAngle)) - x - np.sin(abs(InfectionAngle)))


    if InfectionAngle < 0:
        InitialGuess = 0

    else:
        InitialGuess = InfectionAngle


    try:
        root = optimize.newton(
            function,
            InitialGuess,
            args=(InfectionAngle,Fitness,ObstacleRadius))

    except:
        root = 2 * np.pi

    return root
