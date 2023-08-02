import numpy as np

from collections import defaultdict

import copy

import random

import sys


#############################################################################
#############################################################################
#############################################################################


"""
Class encapsulates everything about the pests.
"""

class Pest:
    """
    A class to represent a single pest, AKA occupied site

    ATTRIBUTES:
    pos: tuple of ints
        The coordinates of the pest.
    parentpos: tuple of ints
        The coordinates of the parent of the pest.
    mut: bool
        The species of the pest: 1=M, 0=WT.
    neighbourlist: list of tuples of ints
        The list of nearest neighbours coordinates.
    emptyneighbourlist: list of tuples of ints
        The list of empty nearest neighbours coordinates.
    active: bool
        Whether the pest is on the population front or not.
    
    METHODS
    Mutation(self,parentmut,MutationProb):
        Sets the mut value of the pest.

    def NeighboursPosList(
            self,ObsRadius,ObsCenterList,DomainHeight,DomainWidth):
        Construct the neighbourlist.

    def EmptyNeighbourPosList(
            self,DomainWidth,Radius,ObsCenterList,ObsCenterDict,
            PopulationDict):
        Construct the emptyneighbourlist.

    def Reproduce(
            self,mutprob,ObsRadius,ObsCenterList,ObsCenterDict,PopulationDict,
            DomainHeight,DomainWidth):
        Replicate the pest into a neighbouring, empty site.

    def UpdateEmptyNeighbours(self,PopulationDict):
        Checks to see if the nearest neighbours, and itself, are still active.
    """
    def __init__(
            self,pos,parentpos,parentmut,parentgenotype,mutprob,ObsRadius,
            ObsCenterList,ObsCenterDict,DomainHeight,DomainWidth,
            PopulationDict):
        """
        Set all the necessary attributes for a Pest.

        ARGS:
        pos: tuple of ints.
            The coordinates of the pest.
        parentpos: tuple of ints
            The coordinates of the parent of the pest.
        parentmut: int
            The type of the parent: 0=WT, 1=M
        mutprob: float
            Mutation probability upon generation of this pest
        ObsRadius: int
            The radius of the patches
        ObsCenterList: list
            A list of the centers of patches
        ObsCenterDict: dict
            A dictionary of regions of the system: within each region we store
            occupied site origins. This removes the need to loop through each.
        DomainHeight: int
            The height of the system
        DomainWidth: int
            The width of the system    
        PopulationDict: dict
            Dict:
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.    
        """
        

        self.pos = pos                          #Position of the pest

        self.parentpos = parentpos              #Position of the pest's parent

        self.mut = self.Mutation(parentmut,mutprob) #If self is a mutant

        self.genotype = parentgenotype

        #Find the neighbours of the pasr
        self.neighbourlist = self.NeighboursPosList(
            ObsRadius,ObsCenterList,DomainHeight,DomainWidth)

        #find he empty neghbours of the pest
        self.emptyneighbourlist = self.EmptyNeighbourPosList(
            DomainWidth,ObsRadius,ObsCenterList,ObsCenterDict,PopulationDict)

        if len(self.emptyneighbourlist) > 0:
            self.active = True
        else:
            self.active = False

    #########################################################################
    #########################################################################

    def Mutation(self,parentmut,MutationProb):
        """
        Determine the type of the pest based on the mutation rtae. If the 
        parent is a Mutant, the child MUST be a mutant

        ARGS:
        parentmut: int
            The type of the parent: 0=WT, 1=M
        MutationProb: float
            Mutation probability upon generation of this pest

        RETRUNS:
        int:
            The 'mut' type of this pest: 0=WT, 1=M
        """
        if parentmut == 1:
            return 1
        else:
            if random.uniform(0,1) < MutationProb:
                return 1
            else:
                return 0
    #########################################################################
    #########################################################################
    def NeighboursPosList(
            self,ObsRadius,ObsCenterList,DomainHeight,DomainWidth):
        """
        Generate all the neighbour sites of the site, on the hexagonal
        lattice.

        ARGS:
        ObsRadius: int
            The radius of the patches.
        ObsCenterList: list
            A list of the centers of patches.
        DomainHeight: int
            The height of the system.
        DomainWidth: int
            The width of the system.

        RETURNS:
        list of tuples of ints
            The coordinates of nearest neighbours coordinates.
        """
        NeighbourList = []
        x = self.pos[0]
        y = self.pos[1]

        """
         . . . . . . . . . 
        . . . . . . . . . . 
         . . . . . . . . . 
        . . . . . . . . . .

         a b
        c ! d
         e f
        """

        ax=bx=cx=dx=ex=fx=0  #The x and y positions of each neighbour cell
        ay=by=cy=dy=ey=fy=0

        if (y % 2) == 0:
            ax = cx = ex = (x-1)%(DomainWidth)
            bx = fx = x
            dx = (x+1)%(DomainWidth)
        else:
            cx = (x-1)%(DomainWidth)
            ax = ex = x
            bx = dx = fx = (x+1)%(DomainWidth)
        

        ay = by = y+1
        cy = dy = y
        ey = fy = y-1

        #GO THROUGH OBSCENTERLIST AND CHECK IF YOU CAN LAND
        if y != (DomainHeight-1):
            NeighbourList.append((ax,ay))
            NeighbourList.append((bx,by))
        NeighbourList.append((cx,cy))
        NeighbourList.append((dx,dy))
        if y != 0:
            NeighbourList.append((ex,ey))
            NeighbourList.append((fx,fy))

        
        return NeighbourList

    #########################################################################
    #########################################################################

    def EmptyNeighbourPosList(
            self,DomainWidth,Radius,ObsCenterList,ObsCenterDict,
            PopulationDict):
        """
        Iterate through the Neighbourlist and identify sites that are
        unoccupied, outputting this as a list of tuples of the coordinates
        of these sites.

        ARGS:
        DomainWidth: int
            The width of the system.
        Radius: int
            The radius of the patches.
        ObsCenterList: list
            A list of the centers of patches.
        ObsCenterDict: dict
            A dictionary of regions of the system: within each region we store
            occupied site origins. This removes the need to loop through each.
        PopulationDict: dict
            A dictionary of pests: the keys are their locations within the 
            system, and the values are the Pest class objects. 

        RETURNS:
        list of tuples of ints:
           The list of empty nearest nrighbours coordinates. 
        """

        NeighbourList = self.neighbourlist.copy()


        xsquaresnum =  ObsCenterDict[("xsquares","ysquares")][0]
        ysquaresnum =  ObsCenterDict[("xsquares","ysquares")][1]

        #Dict Lookup
        removelist = []
        for n in NeighbourList:
            #Obstacles: check if a neighbour point is within the radius of
            #an obstacle. Only care about checking for obstacles if not
            #a mutant 
            if self.mut == 0:
                y = n[1]
                x = n[0]

                if not y%2:
                    x -= 0.5

                x_square = int((x - (x%Radius))/Radius)
                y_square = int((y - (y%Radius))/Radius)

                #loop through nearby squares in the x direction
                for xdictsquare in range(-2,2):
                    xval = int((xdictsquare + x_square)%xsquaresnum)
                    #loop through nearby squares in the y direction
                    for ydictsquare in range(-2,3):
                        yval = int(ydictsquare + y_square)
                        if (yval < ysquaresnum) and (yval >= 0): 
                            for c in ObsCenterDict[(xval,yval)]:
                                centerx = 0
                                if abs(c[0] - x) > DomainWidth/2:
                                    if c[0] - x > 0:
                                        centerx = c[0] - DomainWidth
                                    else:
                                        centerx = c[0] + DomainWidth
                                        
                                else:
                                    centerx = c[0]

                                #This takes the hexagonal offset into account
                                if (((centerx - x)**2 + 0.75*(c[1] - y)**2) <
                                         Radius**2):
                                    removelist.append(n)
                                    xdictsquare += 3
                                    ydictsquare += 5
                                    break
            #Occupied Sites
            if n in PopulationDict:
                removelist.append(n)                     
                        
        
        for r in removelist:
            if r in NeighbourList:
                NeighbourList.remove(r)

        return NeighbourList
    #########################################################################
    #########################################################################

    def Reproduce(
            self,mutprob,ObsRadius,ObsCenterList,ObsCenterDict,PopulationDict,
            DomainHeight,DomainWidth):
        """
        Choose a random empty neighbour, and create a pest at that site.

        ARGS:
        mutprob: float
            Mutation probability upon generation of this pest
        ObsRadius: int
            The radius of the patches.
        ObsCenterList: list
            A list of the centers of patches.
        ObsCenterDict: dict
            A dictionary of regions of the system: within each region we store
            occupied site origins. This removes the need to loop through each.
        PopulationDict: dict
            A dictionary of pests: the keys are their locations within the
            system, and the values are the Pest class objects.
        DomainHeight: int
            The height of the system.
        DomainWidth: int
            The width of the system.


        RETURNS:
        pest:
            The newly generated pest
        """
        neighboursite = random.choice(self.emptyneighbourlist)
        ChildPest = Pest(
            neighboursite,self.pos,self.mut,self.genotype,mutprob,ObsRadius,
            ObsCenterList,ObsCenterDict,DomainHeight,DomainWidth,
            PopulationDict)

        ChildPest.UpdateEmptyNeighbours(PopulationDict)
        self.emptyneighbourlist.remove(neighboursite)
        return ChildPest

    #########################################################################
    #########################################################################

    def UpdateEmptyNeighbours(self,PopulationDict):
        """
        Update the empty neighbour list of pests that have already been 
        placed. If their empty neighbour list is empty, de-active
        the pest, so it is no longer on the front.

        ARGS:
        PopulationDict: dict
            A dictionary of pests: the keys are their locations within the
            system, and the values are the Pest class objects.

        RETURNS:
            NONE
        """
        removelist = []
        for n in self.emptyneighbourlist:
            if n in PopulationDict:
                removelist.append(n)

        for r in removelist:
            self.emptyneighbourlist.remove(r)

        if len(self.emptyneighbourlist) == 0:
            self.active = False

    #########################################################################
    #########################################################################

    #Setters
    def SetMut(self,newmut):
        self.mut = newmut
        return 0
    #Getters
    def GetPos(self):
        return self.pos
    def GetMut(self):
        return self.mut
    def GetNeighbours(self):
        return self.neighbourlist
    def GetEmptyNeighbours(self):
        return self.emptyneighbourlist
    def GetActive(self):
        return self.active
    def GetParentPos(self):
        return self.parentpos
    def GetGenotype(self):
        return self.genotype
##############################################################################
##############################################################################
##############################################################################
##############################################################################
def Initialise(
        mutprob,DomainHeight,DomainWidth,ObsRadius,ObsCenterList,
        ObsCenterDict):
    """
    Creates a population dictionary to contain M and WT occupied sites, and
    two lists to store M and WT sites still on the population front.

    ARGS:
    mutprob: float
        The probability that a given occupied site is a M
    DomainHeight: int
        The height of the system
    DomainWidth: int
        The width of the system
    ObsRadius: int
        The radius of the patches
    ObsCenterList: list
        A list of the centers of patches
    ObsCenterDict: dict
        A dictionary of regions of the system: within each region we store
        occupied site origins. This removes the need to loop through each.

    RETURNS:
    dict
        A dictionary of pests: the keys are their locations within the system,
        and the values are the Pest class objects.
    list
        The coordinates of WT pests on the front.
    list
        The coordinates of M pests on the front.
    """


    UnevolvedList = []
    EvolvedList = []
    PopulationDict = defaultdict(dict)
    for x in range(DomainWidth):
        PopulationDict[(x,0)] = Pest(
            (x,0),(x,0),0,x,mutprob,ObsRadius,ObsCenterList,ObsCenterDict,
            DomainHeight,DomainWidth,{})

    #Making sure they all have correct neighbourlist
    for xy in PopulationDict.keys():
        PopulationDict[xy].UpdateEmptyNeighbours(PopulationDict)
        UnevolvedList.append(xy)        #Because all pests start off as WT

    return PopulationDict,UnevolvedList,EvolvedList


##############################################################################
##############################################################################
##############################################################################
##############################################################################
def Obstacles(Radius,AreaFraction,DomainHeight,DomainWidth):
    """
    Creates a list of patch centers and a dictionary of packaged regions
    of patch centers, for randomly distributed patches about a system of 
    a given height and width, for some anticipated patch area ratio of 
    patch sites.

    ARGS:
    Radius: int
        The radius of the patches.
    AreaFraction: float
        The fraction of the system to be covered in patch sites.
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    

    RETURNS:
    List:
        List of origins of the patches in the system.
    Dict:
        The domain has been split into regions: these regions are keys to the 
        dict. The values of the dict are the patch centers within the regions.
    """
    

    numberdensity = -(1./(np.pi*((Radius**2) * 2./np.sqrt(3.)))) * np.log(
        1.-AreaFraction)     #Units: obstacles per site. hence the 2/root(3)

    ObsNum = numberdensity * (DomainHeight - Radius*2./np.sqrt(3))*DomainWidth

    ObsCenterList = []
    ObsCenterDict = defaultdict(list)


    #Create the list of obstacle center positions, using uniform random values
    for n in range(int(ObsNum)):
        x = random.randint(0,DomainWidth-1)
        y = random.randint(int(Radius*2./np.sqrt(3)) + 2,DomainHeight-1)

        ObsCenterList.append((x,y))
   
    
    #ObsCenterDict is a dictionary that stores the positions of the obstacle
    #centers, but in chuncks so we needn't iterate through every single 
    #obstacle center: we will only need to search through nearby chunks.
    
    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainWidth)/Radius))
    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainHeight)/Radius))


    for O in ObsCenterList:
        x_square = int((O[0] - (O[0]%Radius))/Radius)
        y_square = int((O[1] - (O[1]%Radius))/Radius)
        ObsCenterDict[(x_square,y_square)].append((O[0],O[1]))



    #Output Dict
    for key in ObsCenterDict:
        print("Keys ",str(key[0]),", ",str(key[1]),", leads to list",
            str(ObsCenterDict[key]))

 
    return ObsCenterList, ObsCenterDict

##############################################################################
##############################################################################
##############################################################################
##############################################################################
def LatticeObstacles(
        Radius,DomainHeight,DomainWidth,LatticeType,LatticeSurfSep,
        BottomLeftxdisp,BottomLeftydisp):
    """
    Creates a list of patch centers and a dictionary of packaged regions
    of patch centers, for lattice distributed patches about a system of
    a given height and width, for some lattice type with surface separation.

    ARGS:
    Radius: int
        The radius of the patches.
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    LatticeType: int
        A (currently redundant) orientation of the lattice desired
    LatticeSurfSep: int
        The separation between surfaces of the patches.
    BottomLeftxdisp: int
        The bottommost, leftmost patch surface separation from the left edge
        of the system.
    BottomLeftydisp: int 
        The bottommost, leftmost patch surface separation from the bottom edge
        of the system.

    RETURNS:
    List:
        List of origins of the patches in the system.
    Dict:
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    Int:
        The recalculated width of the domain, to ensure we have a properly
        periodic lattice of patches.
    """

    ObsCenterList = []
    ObsCenterDict = defaultdict(list)

    #The separation between obstacle centers is more useful
    LatticeCentreSep = LatticeSurfSep + 2.0*Radius

    if LatticeType == 0:    #Horizontally adjacent
        DomainWidth = int((2*Radius+LatticeSurfSep)*np.ceil(
            float(DomainWidth)/float(2*Radius+LatticeSurfSep)))

    if LatticeType == 1:    #Vertically adjacent
        DomainWidth = int((np.sqrt(3)*(2*Radius+LatticeSurfSep))*np.ceil(
            float(DomainWidth)/float(np.sqrt(3)*(2*Radius+LatticeSurfSep)))) 

    


    if LatticeType == 0:    #Shortest dist is between horizontal neighbours
        xObstaclesNum = int(np.floor(
            float(DomainWidth-BottomLeftxdisp)/LatticeCentreSep) )
        yObstaclesNum = int(np.floor(
            float(DomainHeight-BottomLeftydisp)/LatticeCentreSep))

        
        for obsrow in range(yObstaclesNum):
            ypos = int(np.ceil(
                BottomLeftydisp+(2/np.sqrt(3))*Radius+obsrow*
                LatticeCentreSep))
            for obscol in range(xObstaclesNum):
                xpos = int(np.ceil(
                    BottomLeftxdisp + Radius + (obscol+0.5*(obsrow%2))*
                    LatticeCentreSep))
                ObsCenterList.append((xpos,ypos))
    
    else:
        print("Have not implemented this lattice type yet: please use type 0")
        sys.exit()


    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainWidth)/Radius))
    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainHeight)/Radius))


    for O in ObsCenterList:
        x_square = int((O[0] - (O[0]%Radius))/Radius)
        y_square = int((O[1] - (O[1]%Radius))/Radius)
        ObsCenterDict[(x_square,y_square)].append((O[0],O[1]))



    #Output Dict
    for key in ObsCenterDict:
        print("Keys ",str(key[0]),", ",str(key[1]),", leads to list",
            str(ObsCenterDict[key]))


    return ObsCenterList, ObsCenterDict, DomainWidth

##############################################################################
##############################################################################
##############################################################################
##############################################################################
def SingleObstacle(Radius,DomainHeight,DomainWidth,xpos,ypos):
    """
    Creates a list of patch centers and a dictionary of packaged regions
    of patch centers, for a single patch within a system of a given height
    and width, for some lattice type with surface separation.

    ARGS:
    Radius: int
        The radius of the patches.
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    xpos: int
        The x-coordinate of the center of the patch.
    ypos: int:
        the y-coordinate of the center of the patch.

    RETURNS:
    List:
        List of origins of the patches in the system.
    Dict:
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.

    """

    ObsCenterList = []
    ObsCenterDict = defaultdict(list)

    ObsCenterList.append((xpos,ypos))

    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainWidth)/Radius))
    ObsCenterDict[("xsquares","ysquares")].append(
        np.ceil(float(DomainHeight)/Radius))

    for O in ObsCenterList:
        x_square = int((O[0] - (O[0]%Radius))/Radius)
        y_square = int((O[1] - (O[1]%Radius))/Radius)
        ObsCenterDict[(x_square,y_square)].append((O[0],O[1]))

    return ObsCenterList, ObsCenterDict



##############################################################################
##############################################################################
##############################################################################
##############################################################################
def RoughFront(
        DomainHeight,DomainWidth,mutprob,fitness,PopulationDict,UnevolvedList,
        EvolvedList,ObsRadius,ObsCenterList,ObsCenterDict,EndCondition):
    """
    Execture a single step of the Rough Front Eden model, taking as input,
    updating then outputting various lists and dictionary regarding occupied
    sites within the system. A member of the fron tis selected to infect a
    nearest unoccupied neighbour: during this process it may mutate if it is
    a WT. The lists and dicts are updated to ascertain if any local sites on
    the front are now not on the front.

    ARGS:
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    mutprob: float
        The probability that a given occupied site is a M.
    fitness: float
        The fitness of M: the ratio to which a M is chosen relative to a WT.
    PopulationDict: dict
        A dictionary of pests: the keys are their locations within the system,
        and the values are the Pest class objects.
    UnevolvedList: list of tuples
        The coordinates of WT pests on the front.
    EvolvedList: list of tuples
        The coordinated of M pests on the front.
    ObsRadius: int
        The radius of the patches.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.    
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    EndCondition: int
        This is the option corresponding to the condition for which we will
        cease executing the rough front Eden Model.


    RETURNS:
    dict:
        Updated dictionary of pests: the keys are their locations within the
         system, and the values are the Pest class objects.
    list of tuples:
        Updated list of coordinates of M pests on the front. 
    list of tuples:
        Updates list of coordinates of WT pests on the front.
    bool:
        Whether the Rough Front Eden model step should be executed again
    Pest:
        The child pest that was created during this step.
    """
    


    #Randomly Select
    #Weight prob of being chosen against the fitness
    ran = random.uniform(0,float(len(UnevolvedList)+fitness*len(EvolvedList)))

    ParentPest = 0  #Assume WT by default

    #WT Chosen
    if ran < len(UnevolvedList):

        #Select random member keys:
        rand = (random.choice((UnevolvedList)))

        ParentPest = PopulationDict[rand]  

    else:
        #Select random member keys:
        rand = (random.choice((EvolvedList)))
        ParentPest = PopulationDict[rand]

    #Create New
    ChildPest = ParentPest.Reproduce(
        mutprob,ObsRadius,ObsCenterList,ObsCenterDict,PopulationDict,
        DomainHeight,DomainWidth)

    #Add to population
    ChildPos = ChildPest.GetPos()
    PopulationDict[ChildPos] = ChildPest

    if ChildPest.GetActive() == True:
        if ChildPest.GetMut() == 0:
            UnevolvedList.append(ChildPos)
        else:
            EvolvedList.append(ChildPos)


    #Update Lists
    #Cycle through neighbours of the child to see if any neighbours are empty
    ChildNeighbours = ChildPest.GetNeighbours()

    #Check if neighbours are still on the front
    for n in ChildNeighbours:
        if n in PopulationDict:
            NeighbourPest = PopulationDict[n]       #Choose site
            if not (type(NeighbourPest) is int):    #If it is pest-typ object
                if NeighbourPest.GetActive() == True:

                    #Update number of neighbours
                    NeighbourPest.UpdateEmptyNeighbours(PopulationDict)
                    
                    #If the pest is no longer active, delete them
                    if NeighbourPest.GetActive() == False:
                        NeighbourPos = NeighbourPest.GetPos()
                        if NeighbourPest.GetMut() == 0: 
                            UnevolvedList.remove(NeighbourPos)
                        else:
                            EvolvedList.remove(NeighbourPos)
    
                        PopulationDict[NeighbourPos] = NeighbourPest.GetMut()

    if not ChildPest.GetActive():
        PopulationDict[ChildPos] = ChildPest.GetMut()


    #Check for end-conditions
    EndConditionMet = False

    #If there are no more sites on the front
    if EndCondition == 0:   
        if (len(EvolvedList) + len(UnevolvedList)) == 0:
            EndConditionMet = True

    #Child reaches domain height. Also stop if there can be no more expansion.
    elif EndCondition == 1: 
        if ChildPos[1] == DomainHeight-2:
            EndConditionMet = True
        
        if (len(EvolvedList) + len(UnevolvedList)) == 0:
            EndConditionMet = True

    #When the first child reaches the ehight of 'EndCondition'
    else:
        if ChildPos[1] == EndCondition:
            EndConditionMet = True
    
        if (len(EvolvedList) + len(UnevolvedList)) == 0:
            EndConditionMet = True

    return PopulationDict,EvolvedList,UnevolvedList,EndConditionMet,ChildPest

##############################################################################
##############################################################################
##############################################################################
##############################################################################
def ExtendSystem(
        Radius,AreaFraction,DomainHeight,DomainWidth,ObsCenterList,
        ObsCenterDict,ExtensionHeight):
    """
    Increase the height of the system, adding more randomly distributed
    patches in the newly generated space.

    ARGS:
    Radius: int
        The radius of the patches.
    AreaFraction: float
        The fraction of the system to be covered in patch sites.
    DomainHeight: int
        The height of the system.
    DomainWidth: int
        The width of the system.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    ExtensionHeight: int:
        The extra height we wish to add to the domain.

    RETURNS:
    list:
        Updated list of coordinates of the patches in the system.
    dict:
        Updated dict in which domain has been split into regions: these
        regions are keys to the dict. The values of the dict are the 
        patch centers within the regions.
    int:
        Updated height of the system.
    """


    numberdensity = -(1./(np.pi*(Radius**2) * 2./np.sqrt(3.))) * np.log(
        1.-AreaFraction)

    ObsNum = numberdensity * ExtensionHeight*DomainWidth

    for n in range(int(ObsNum)):
        x = random.randint(0,DomainWidth-1)
        y = random.randint(DomainHeight,DomainHeight + ExtensionHeight)
        ObsCenterList.append((x,y))    

        x_square = int((x - (x%Radius))/Radius)
        y_square = int((y - (y%Radius))/Radius)
        ObsCenterDict[(x_square,y_square)].append((x,y))    

    DomainHeight += ExtensionHeight

    ObsCenterDict[("xsquares","ysquares")][1] = (np.ceil(
        float(DomainHeight)/Radius))


    return ObsCenterList,ObsCenterDict,DomainHeight

##############################################################################
##############################################################################
##############################################################################
##############################################################################
def MeasureAreaFraction(
        DomainWidth,DomainHeight,ObstacleRadius,ObsCenterList,ObsCenterDict):
    """
    Measure the area fraction of the system, with randomly distributed
    patches taking the ratio of the patch area and the total area of the 
    system.

    ARGS:
    DomainWidth: int
        The width of the system.
    DomainHeight: int
        The height of the system.
    ObstacleRadius: int
        The radius of the patches.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.

    RETURNS:
    float:
        The ratio of the patch area and the total area of the system.
    """


    TreatedArea = 0

    xsquaresnum =  ObsCenterDict[("xsquares","ysquares")][0]
    ysquaresnum =  ObsCenterDict[("xsquares","ysquares")][1]


    for i in range(DomainWidth):
        #for n,j in enumerate(y):
        for j in range(int(ObstacleRadius*2./np.sqrt(3)),DomainHeight):
            if not(j%2):    #if on even rows
                xpos = i-0.5
            else:
                xpos = i
            ypos = j * np.sqrt(3.)/2.0

            x_square = int((i - (i%ObstacleRadius))/ObstacleRadius)
            y_square = int((j - (j%ObstacleRadius))/ObstacleRadius)


            xdictsquare = -2
            
            while (xdictsquare < 2):
                xval = int((xdictsquare + x_square)%xsquaresnum)

                ydictsquare = -2
                while (ydictsquare < 3):
                    yval = int(ydictsquare + y_square)
                    if (yval < ysquaresnum) and (yval >= 0):
                        for c in ObsCenterDict[(xval,yval)]:
                            centerx = 0
                            if abs(c[0] - xpos) > DomainWidth/2:
                                if c[0] - xpos > 0:
                                    centerx = c[0] - DomainWidth
                                else:
                                    centerx = c[0] + DomainWidth

                            else:
                                centerx = c[0]
                            if (((centerx - xpos)**2 + 0.75*(c[1] - j)**2) < 
                                    ObstacleRadius**2):
                                TreatedArea += 1
                                xdictsquare += 3
                                ydictsquare += 5
                                break
                    ydictsquare += 1
                xdictsquare += 1


    return TreatedArea/(
        (DomainHeight - ObstacleRadius*2./np.sqrt(3))*DomainWidth)


##############################################################################
##############################################################################
##############################################################################
##############################################################################
def MeasureAreaFractionLattice(
        DomainWidth,DomainHeight,ObstacleRadius,ObsCenterList,ObsCenterDict,
        BottomLeftydisp):
    """
    Measure the area fraction of the system, with randomly distributed
    patches taking the ratio of the patch area and the total area of the
    system.

    ARGS:
    DomainHeight: int
        The width of the system.
    DomainWidth: int
        The height of the system.
    ObstacleRadius: int
        The radius of the patches.
    ObsCenterList: list of tuples:
        List of coordinates of the patches in the system.
    ObsCenterDict: dict
        The domain has been split into regions: these regions are keys to the
        dict. The values of the dict are the patch centers within the regions.
    BottomLeftydisp: int
        The bottommost, leftmost patch surface separation from the bottom edge
        of the system.

    RETURNS:
    float:
        The ratio of the patch area and the total area of the system.
    """


    TreatedArea = 0

    xsquaresnum =  ObsCenterDict[("xsquares","ysquares")][0]
    ysquaresnum =  ObsCenterDict[("xsquares","ysquares")][1]


    for i in range(DomainWidth):
        #for n,j in enumerate(y):
        for j in range(BottomLeftydisp,
                ObsCenterList[-1][1]+int(ObstacleRadius*2./np.sqrt(3))):
            if not(j%2):    #if on even rows
                xpos = i-0.5
            else:
                xpos = i
            ypos = j * np.sqrt(3.)/2.0

            x_square = int((i - (i%ObstacleRadius))/ObstacleRadius)
            y_square = int((j - (j%ObstacleRadius))/ObstacleRadius)


            xdictsquare = -2

            while (xdictsquare < 2):
                xval = int((xdictsquare + x_square)%xsquaresnum)

                ydictsquare = -2
                while (ydictsquare < 3):
                    yval = int(ydictsquare + y_square)
                    if (yval < ysquaresnum) and (yval >= 0):
                        for c in ObsCenterDict[(xval,yval)]:
                            centerx = 0
                            if abs(c[0] - xpos) > DomainWidth/2:
                                if c[0] - xpos > 0:
                                    centerx = c[0] - DomainWidth
                                else:
                                    centerx = c[0] + DomainWidth

                            else:
                                centerx = c[0]
                            if (((centerx - xpos)**2 + 0.75*(c[1] - j)**2) <
                                     ObstacleRadius**2):
                                TreatedArea += 1
                                xdictsquare += 3
                                ydictsquare += 5
                                break
                    ydictsquare += 1
                xdictsquare += 1


    return TreatedArea/(
        (ObsCenterList[-1][1]+int(ObstacleRadius*2./np.sqrt(
        3))-BottomLeftydisp)*DomainWidth)
