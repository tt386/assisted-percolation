# Code for *Assisted percolation of slow-spreading mutants in heterogeneous environments*


## Overview

[Publication](#publication)

[Brief description of the generalized Eden Model](#brief-description-of-the-generalized-eden-model)

[Structure of simulation code](#structure-of-simulation-code)

[Directory structure and exectuing code](#directory-structure-and-executing-code)

[Reproducing figures](#reproducing-figures)


## Publication

Manuscript published in XXX in 2023 by Thomas Tunstall, Tim Rogers, and Wolfram Möbius. The corresponding preprint can be found on arXiv XXX.

We have included all data files and results pertinent to the manuscript (except those larger than 100MB), hence why this repository is so large.
## Brief description of the generalized Eden Model

A full description can be found in the paper and supplementary information.

A hexagonal lattice of sites that come in one of 4 types:

* Unoccupied, 'Background'&nbsp;&nbsp;&nbsp;(<span style="color:#CCCCCC">light grey</span>)
* Unoccupied, 'Patch'&nbsp;&nbsp;&nbsp;(<span style="color:#666666">Dark grey</span>)
* Occupied, Wild-Type&nbsp;&nbsp;&nbsp;(<span style="color:#33CCFF">Blue</span>)
* Occupied, Mutant&nbsp;&nbsp;&nbsp;(<span style="color:#CC0033">Red</span>)

### Initial conditions

Patch sites are arranged into circles of radius $R$ called patches. They are placed either individually, randomly (allowed to overlap) or onto a fixed lattice. The total area ratio of patch sites is defined to be $\phi$. A number of patches, $n$, are placed such that a target patch area ratio is achieved, related by:
$$\phi = 1- e^{-n\pi R^2 }.$$
The bottom of the lattice is initialised as a line of WT sites.

### Update Rules

At each iteration of the model, the system is updated according to the following rules:

* A random occupied site, which has neighbours that it can occupy, is chosen
    * WT can only occupy background sites, whilst M can occupy both background and patch sites.
    * The probability of a M being selected is weighted by a 'fitness' factor $F$, where $F<1$ throughout this work.
* For this chosen site, a random nearest neighbour site (which can be occupied by the chosen site) is selected.
* The selected neighbouring site is then updated to be of the same type as the intially selected site. However, if the initial site was WT, where is a fixed mutation probability $\mu$ that the neighbour site becomes of type M.

A simulation terminates when some end condition is satisfied as detailed in the publication.


## Structure of simulation code

The code structure here is similar across all simulations with the most general example being `Scripts_OneNode/SystemEvolution/Script.py`. The following gives an overview of that script.

### Initialisation

    ObsCenterList,ObsCenterDict = [],{}

    if Random:
        ObsCenterList,ObsCenterDict = CoreSim.Obstacles(Radius,AreaFraction,DomainHeight,DomainWidth)

    elif Lattice:
        ObsCenterList,ObsCenterDict,DomainWidth = CoreSim.LatticeObstacles(Radius,DomainHeight,DomainWidth,LatticeType,LatticeSurfSep,BottomLeftxdisp,BottomLeftydisp)

    else:
        ObsCenterList,ObsCenterDict = CoreSim.Obstacles(Radius,0,DomainHeight,DomainWidth)

    PopulationDict,UnevolvedList,EvolvedList = CoreSim.Initialise(0,DomainHeight,DomainWidth,Radius,ObsCenterList,ObsCenterDict)

`ObsCenterList` corresponds to a list of all patch centers. (Note: When developing the code, patches were called obstacles, thus the choice of variable names and folder names.) `ObsCenterDict` corresponds to a dictionary of squares (part of domain) of the system, within which the patch center coordinates within that square are stored. This auxiliary variable increases computational efficiency. Patch centers are initialised by `CoreFunctions/CoreSim.py`, with the options of `Obstacles` for random patches, `LatticeObstacles` for patches on a lattice, and `SingleObstacle` for a single patch.

`PopulationDict` is a record of all occupied sites: it is a dictionary with keys equivalent to the spatial position of the sites and stores either an object (if the occupied site is active - i.e., it can still invade a neighbouring site) or an integer (if the occupied site is inactive). This is useful for interating over sites to see if an occupied site can invade another site. `UnevolvedList` and `EvolvedList` are lists that store the positions of WT and M occupied sites that are active, respectively. These are the lists that are sampled to determine the next member of the front that will replicate.

### Main loop


    EndConditionMet = False
    while EndConditionMet == False:
        PopulationDict, EvolvedList, UnevolvedList, EndConditionMet, ChildPest = CoreSim.RoughFront(DomainHeight,DomainWidth,mutprob,fitness,PopulationDict,UnevolvedList,EvolvedList,Radius,ObsCenterList,ObsCenterDict,EndCondition)


        #Extend System if Necessary
        if (DomainHeight - ChildPest.GetPos()[1]) < (np.sqrt(3)/2)*Radius:
            ObsCenterList,ObsCenterDict,DomainHeight = CoreSim.ExtendSystem(Radius,AreaFraction,DomainHeight,DomainWidth,ObsCenterList,ObsCenterDict,ExtensionHeight)

        #Emergency End Condition
        if ChildPest.GetPos()[1] >= UpperDomainLimit-1:
            EndConditionMet = True

`EndConditionMet`, a Boolean, is constructed to keep track of when to stop the evolution of the system.

The most important aspect is `RoughFront`. Within a single iteration it selects a member of the population front to reproduce. Has that member infected a random nearest neighbour, then it checks to make sure the nearby member of the front are still indeed on the front. Therefore, `PopulationDict`, `EvolvedList`, and `UnevolvedList` are updated every loop, and the inputs regarding M behaviour (fitness and mutation rate) and patches (`ObsCenterList`, `ObsCenterDict`) are required inputs.

The `EndCondition` variable decides the end condition of the system. For example, it can stop the system when the entire front becomes M, or when there are no more members of the population front. This triggers `EndConditionMet` to become true. A full list:
-   0: There is no more expansion due to there being no more occupied
    sites on the population frontier
-   1: The first occupied site reaches the height of the domain. This
    has found reduced utility since the system has been capable of
    expanding.
-   else: when the height of the most recent occupied site produced is
    at this height, the system stops. A natural upper limit to the
    system.

Finally in our main loop we have extra functions that can be script-specific. In `ExtendSystem`, there is a way of generating new patches if our most recent occupied site is too close to the maximum current height of the system.

There are definable end conditions that are oftentime more specific than those previously listed. In the code given, it is just a check that out system doesn’t exceed some UpperDomainLimit.

### Output

Data are saved to an NPZ file named `timecourse.npz` (in directory `SystemEvolution`) or `DataFile.npz` (elsewhere).

    OutputDatafilename = SaveFileDirName + '/timecourse.npz'
    np.savez(OutputDatafilename,
        OutputData=OutputData,
        DomainHeight=DomainHeight,
        DomainWidth=DomainWidth,
        fitness=fitness,
        mutprob=mutprob,
        Radius=Radius,
        ObsCenterList=ObsCenterList,
        DataSizeList=DataSizeList)

There is a corresponding plotting file within the same directory reading in those data.


## Directory structure and executing code

### Dependencies

Create the environment using ```conda env create -f EdenEnv.yml```. [[Manual on setting up environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)]

### Directory structure

Below is a tree respresenting the structure of directories. Scripts that run using a single node are found in ```Scripts_OneNode```. Scripts that are meant to run on HPC (using SLURM) are found in ```Scripts_HPC```. Scripts which only create visualisations of geometrical aspects of the theory are found in ```Sketch```.

```
.
├── CoreFunctions
├── RawFigures
├── Scripts_HPC
│   ├── EscapeRegion_HeightVsAngle
│   ├── EscapeRegion_HeightVsFitness
│   ├── InfectionProbability
│   ├── PhaseDiagram_Lattice
│   └── PhaseDiagram_Random
├── Scripts_OneNode
│   ├── ClusterHeightStatistics
│   ├── FlatFrontEvolution
│   └── SingleSystemEvolution
└── Sketch
    ├── EscapeRegionCalculation
    ├── EscapeRegionCalculation_Over
    ├── Lattice_LongTerm
    ├── Lattice_Orientation
    ├── Lattice_ShortTerm
    └── SymmetricEscapeRegion
```

### Directory `Scripts_OneNode`

Directories contain code that runs on a single core. To execute them within their home directories execute

    $ python Script.py
    $ python Plotting.py -d SaveFiles/Dirname

When executing `Script.py`, there is an optional argument in which one can add a random seed:

    $ python Script.py -r INTEGER

If no random seed is called, one is automatically generated. The random seed can be read in the output directory in the `randomseed.txt` file and used for rerunning the code.

#### Directory `Scripts_OneNode/SingleSystemEvolution`

This produces an .npz file that can be used to create many plots of the system, for each of a predefined timestep defined in `Plotting.py`. An example is given below.

![SampleOutput for `SystemEvolution`](RawFigures/Fig_1_C_ii.png)

These images can then be concatenated into a video by entering the directory with the images:

    ffmpeg -r 20 -pattern_type glob -i 'AField_*.png' -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf scale=-2:1240 SystemEvolution.mp4

In this directory, the `Plotting.py` has two optional arguments:

- `python Plotting.py -b 1` means that only the patches will be plotted ('b' stands for 'background', and the '1' means `True`)
- `python Plotting.py -n 2 4 105` means that only the 2nd, 4th and 105th frames will be plotted. These numbers can be any integer.  


Note: In this repository we have endeavoured to include all data files used to generate figures. The exception to this are the Random and Lattice patch distributions with a (height,width) of (2000,2000): these `.npz` files exceeded 100MB, the limit Github will accept.

#### Directory `Scripts_OneNode/ClusterHeightStatistics`

This directory is a bit different from the rest, in such that it does not have a single `Script.py`: instead it has a script for both rough and flat fronts: `RoughFrontClusterHeights.py` and `FlatFrontClusterHeights.py`, respectively.

`RoughFrontClusterHeights.py` takes input parameters entirely from the `Params.py` file, creates a directory and executes many repeats of a small rough front system, in which M is seeded. The maximum and minimum height extents of the resulting cluster is monitored until there are no more M on the front. The data is saved to a created directory for plotting a cluster height distribution. Executing for the sample output will take a few hours.

`FlatFrontClusterHeights.py` is run separately, and has three optional arguments: the M fitness, the number of repeats, and the name of the directory to save data (with argument keys `-f`, `-R` and `-d` respectively). If only the first two of these arguments are included, a unique `SaveFiles` directory is created, and the corresponding data is saved here. If any argument is missing, it will use the values recorded in the `Params.py` file. Executing for the sample output will take roughly a minute.

The `Plotting.py` file takes a directory argument and plots both the rough and flat front outputs on the same plot, along with the analytical result for the flat front. Furthermore, if there is a datafile present for the flat front, another directory is created and, for a number of cluster heights, the distribution of cluster widths for the flat front is determined; each of these are outputted as a separate png.

[[PDF for a cluster height distribution](./RawFigures/Fig_2_A.pdf)]

[[PDF for a cluster width distribution for a given height](./RawFigures/SuppFig_2_B_i.eps)]


## Directory `Scripts_HPC`

Within this directory are all the directories and scripts making use of the Exeter HPC ISCA (using Slurm Workload Manager) to typically run many repeats per node across many nodes, each with different parameters for the system. Therefore, there is additional complexity within these directories:

The start of `CreateParams.py` the parameter values and combinations thereof to be used. For 2D maps, there are two lists of parameter values.

    \$ python CreateParams.py

creates a `ParamPairs.npz` file, which includes a list of all combinations of each parameter - this is to be iterated over! Any previous list is deleted and a directory in `SaveFiles` is created to store results when running the actual simulation. The number $n$ in the output corresponds to the number of parameters pairs, which corresponds to the number of jobs to run.

There is an optional argument which defines the random seed of the simulation set:

    $ python CreateParams.py -r INTEGER

If a seed is not defined, a seed is generated and can be found in ```randomseed.txt```.

The number of parameter pairs $n$ is used to run the
`Parallel_Bash.sh` script to send off the jobs:

    $ sbatch --array=0-n Parallel_Bash.sh

The random seed of each parameter pair is taken to be the sum of the simulation set random seed and the index of the parameter pair.

An example of the full process:

    $ python CreateParams.py
    Number of Parameter Pairs (therefore, last number in array)360
    Deleted Previous
    $ sbatch --array=0-360 Parallel_Bash.sh

Data are plotted as before:

    $ python Plotting.py -d SaveFiles/Dirname

### Directory `Scripts_HPC/EscapeRegion_HeightVsAngle`

The system is executed for a system with a single patch. Mutation is disabled; instead, a site at a certain angle (given by the parameter passed to the job) is seeded to become M upon being infected. The system is stopped when there are no mutants on the population front. We keep track of the highest M site: This corresponds to the maximum height of the escape region.

[[PDF for a distribution fo escape region heights with changing invasion angle4](./RawFigures/Fig_3_D.pdf)]

### Directory `Scripts_HPC/InfectionProbability`

Taking the mutation rate as a parameter, the probability that a patch of fixed radius becomes invaded by mutants is computed for a fixed fitness. The system is rerun for different fitnesses.

[[PDF for a distribution of probability of patch invasion with changing mutation prob](./RawFigures/Fig_2_B.pdf)]

### Directory `Scripts_HPC/PhaseDiagram_Lattice`

The parameters are fitness of M and patch area fraction. The aim is to find the probability M domination given a single mutation from WT to M.

The mutation rate is set to 0 at the beginning to allow for the front to roughen, then set to a finite value for a single mutation to occur randomly, and then set to 0 again. For each repeat of the system, there are a number of ‘sub-repeats’. This is to ensure that at least one patch was invaded. In addition, repeats might be discarded to keep simulation time within reasonable bounds. See supplementary information to the manuscript for details.

[[PDF for the lattice patch distribution phase diagram](./RawFigures/Fig_5_B.eps)]

### Directory `Scripts_HPC/PhaseDiagram_Random`

This is very similar to the case of patches organised on a lattice, except that for each repeat a new distribution of patches is used. The same location of patches is used across sub-repeats.

[[PDF for the random patch distribution phase diagram](./RawFigures/Fig_4_B.eps)]


# Reproducing Figures 

The subfigures or panels used in the paper can be found in the directory ```RawFigures```. Script ```Bash_Copying.sh``` can be used to collect these from the individual folders where they were generated:

    bash Bash_Copying.sh RawFigures/

The original locations of the subfigures are listed in ```Bash_Copying.sh```. These source directories also include the files ```Params.py``` and ```randomseed.txt``` which contain parameters and random seed used.

To recreate a given subfigure, move into the appropriate folder and execute the commands detailed below. Note that the directory created to record results will not be identical to the folder from which the parameter file is taken. The name of the directory is found in the plotting part of the instructions.

Also note that in the folder `Scripts_OneNode` cases, overwriting a directory is prevented to ensure no data is accidentally overwritten.

## `Scripts_OneNode/SingleSystemEvolution`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [1B_i](./RawFigures/Fig_1_B_i.png) | `cp SaveFiles/Fig_1B_i_Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.300/Params.py .` <br> `python Script.py -r 6195044168430747995` <br> `python Plotting.py -d SaveFiles/Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.300 -b 1` |
| [1B_ii](./RawFigures/Fig_1_B_ii.png) | `cp SaveFiles/Fig_1B_ii_Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.600/Params.py .` <br> `python Script.py -r 6844594242896927532` <br> `python Plotting.py -d SaveFiles/Height_1000_Width_1000_Fitness_0.700_MutProb_0.00100_Radius_50_Random_AreaFraction_0.600 -b 1` |
| [1C_i](./RawFigures/Fig_1_C_i.png), [1C_ii](./RawFigures/Fig_1_C_ii.png) | `cp SaveFiles/Fig_1C_Height_1500_Width_1000_Fitness_0.950_MutProb_0.00100_Radius_50_Random_AreaFraction_0.500/Params.py .` <br> `python Script.py -r 5404145099506580897` <br> `python Plotting.py -d SaveFiles/Height_1500_Width_1000_Fitness_0.950_MutProb_0.00100_Radius_50_Random_AreaFraction_0.500 -n 75 150` |
| [2C_i](./RawFigures/Fig_2_C_i_a.png), [3A_i](./RawFigures/Fig_3_A_i_a.png) | `cp SaveFiles/Fig_2C_i_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/Params.py .` <br> `python Script.py -r 1798679661298138938` <br> `python Plotting.py -d SaveFiles/Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle` |
| [2C_ii](./RawFigures/Fig_2_C_ii_a.png), [3A_ii](./RawFigures/Fig_3_A_ii_a.png) | `cp SaveFiles/Fig_2C_ii_Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle/Params.py .` <br> `python Script.py -r 8349812612310990994` <br> `python Plotting.py -d SaveFiles/Height_400_Width_400_Fitness_0.900_MutProb_0.00050_Radius_100_SingleObstacle` |
| [4A_i](./RawFigures/Fig_4_A_i.png), [S1B_i](./RawFigures/SuppFig_1_C_i.png) | `cp SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/Params.py .` <br> `python Script.py -r 1536085347793381594` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400 -n 200` |
| [S1B_ii](./RawFigures/SuppFig_1_C_ii.png) | `cp SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.02000_Radius_40_Random_AreaFraction_0.400/Params.py .` <br> `python Script.py -r 2701395967057991801` <br> `python Plotting.py -d Height_2000_Width_2000_Fitness_0.900_MutProb_0.02000_Radius_40_Random_AreaFraction_0.400 -n 200` |
| [4A_ii](./RawFigures/Fig_4_A_ii.png), [S1B_iii](./RawFigures/SuppFig_1_C_iii.png) | `cp SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400/Params.py .` <br> `python Script.py -r 5538001929365847712` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Random_AreaFraction_0.400 -n 200` |
| [4A_iii](./RawFigures/Fig_4_A_iii.png), [S1B_iiii](./RawFigures/SuppFig_1_C_iiii.png) | `cp SaveFiles/RandHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.700/Params.py .` <br> `python Script.py -r 530881470231357024` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Random_AreaFraction_0.700 -n 200` |
| [5A_i](./RawFigures/Fig_5_A_i.png) | `cp SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/Params.py .` <br> `python Script.py -r 6075686819584632846` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400 -n 200` |
| [5A_ii](./RawFigures/Fig_5_A_ii.png) | `cp SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400/Params.py .` <br> `python Script.py -r 6704889587093784354` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.980_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_41.000_AreaFraction_0.400 -n 200` |
| [5A_iii](./RawFigures/Fig_5_A_iii.png) | `cp SaveFiles/LattHeatMap_Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_6.000_AreaFraction_0.800/Params.py .` <br> `python Script.py -r 9063097149324546647` <br> `python Plotting.py -d SaveFiles/Height_2000_Width_2000_Fitness_0.900_MutProb_0.00100_Radius_40_Lattice_LatticeSurfSep_6.000_AreaFraction_0.800 -n 200` |
| [S1A](./RawFigures/SuppFig_1_B.png) | `cp SaveFiles/ZoomedRoughFront_Height_30_Width_30_Fitness_0.800_MutProb_0.01000_Radius_5_SingleObstacle/Params.py .` <br> `python Script.py -r 6819485700400171543` <br> `python Plotting.py -d SaveFiles/Height_30_Width_30_Fitness_0.800_MutProb_0.01000_Radius_5_SingleObstacle` |


## `Scripts_OneNode/ClusterHeightStatistics`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [2A](./RawFigures/Fig_2_A.pdf), [S2C_i](./RawFigures/SuppFig_2_B_i.eps), [S2C_ii](./RawFigures/SuppFig_2_B_ii.eps), [S2C_iii](./RawFigures/SuppFig_2_B_iii.eps) | `cp SaveFiles/SampleOutput/Params.py .` <br> `python RoughFrontClusterHeights.py -r 4757963275543371788` <br> `python FlatFrontClusterHeights.py -r 5943095296882043323 ` <br> `python Plotting.py -d Height_500_Width_100_Fitness_0.900_Repeats_100000_SeededHeight_50` |

## `Scripts_OneNode/FlatFrontEvolution`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [S2A](./RawFigures/SuppFig_2_A.png) | `cp SaveFiles/SampleOutput/Params.py .` <br> `python Script.py -r 8631765320194566459` <br> `python Plotting.py -d SaveFiles/Height_30_Width_30_Fitness_0.800_MutProb_0.01000_NoObstacles` |

## `Scripts_HPC/InfectionProbability`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [2B](./RawFigures/Fig_2_B.pdf) | `cp SaveFiles/Height_280_Width_180_Repeats_10000_Radius_60_BottomBuffer_100_FlatInfection_RoughFront/Params.py .` <br> `python CreateParams.py -r 788091261247258457` <br> `sbatch --array=0-150 Parallel_Bash.sh` <br> `python Plotting.py -d SaveFiles/Height_280_Width_180_Repeats_10000_Radius_60_BottomBuffer_100_FlatInfection_RoughFront` |

## `Scripts_HPC/EscapeRegion_HeightVsFitness`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [3C](./RawFigures/Fig_3_C.pdf) | `cp SaveFiles/Height_1000_Width_150_Repeats_1000_Radius_50_Min_Fitness_0.700_Max_Fitness_0.990_num_Fitness_100_FlatInfection_RoughFront_Random_Radius_50/Params.py .` <br> `python CreateParams.py -r 8392284138695379699` <br> `sbatch --array=0-100 Parallel_Bash.sh` <br> `python Plotting.py -d SaveFiles/Height_1000_Width_150_Repeats_1000_Radius_50_Min_Fitness_0.700_Max_Fitness_0.990_num_Fitness_100_FlatInfection_RoughFront_Random_Radius_50` |

## `Scripts_HPC/EscapeRegion_HeightVsAngle`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [3D](./RawFigures/Fig_3_D.pdf) | `cp SaveFiles/Height_1000_Width_150_Repeats_1000_Radius_50_Fitness_0.950_FlatInfection_RoughFront_Random_Radius_50/Params.py .` <br> `python CreateParams.py -r 6942065374998542844` <br> `sbatch --array=0-572 Parallel_Bash.sh` <br> `python Plotting.py -d SaveFiles/Height_1000_Width_150_Repeats_1000_Radius_50_Fitness_0.950_FlatInfection_RoughFront_Random_Radius_50` |

## `Scripts_HPC/PhaseDiagram_Random`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [4C](./RawFigures/Fig_4_B.eps) | `cp SaveFiles/SampleOutput/Params.py .` <br> `python CreateParams.py -r 1919149702158833669` <br> `sbatch --array=0-306 Parallel_Bash.sh` <br> `python Plotting.py -d SaveFiles/DomainHeight_10000_DomainWidth_2000_Repeats_50_FlatInfection_RoughFront_Random_Radius_50` |

## `Scripts_HPC/PhaseDiagram_Lattice`

| Figure(s) | Commands for simulation and creating figure |
| ------------- | ------------- |
| [5B](./RawFigures/Fig_5_B.eps) | `cp SaveFiles/SampleOutput/Params.py .` <br> `python CreateParams.py -r 1322269423698166977` <br> `sbatch --array=0-360 Parallel_Bash.sh` <br> `python Plotting.py -d SaveFiles/DomainHeight_10000_DomainWidth_2000_Repeats_50_FlatInfection_RoughFront_Lattice_Radius_40` |

## `Sketch/EscapeRegionCalculation`

| Figure(s) | Command for creating figure |
| ------------- | ------------- |
| [S3A](./RawFigures/SuppFig_3_A.pdf) | `python Plotting.py` |

## `Sketch/EscapeRegionCalculation_Over`

| Figure(s) | Command for creating figure |
| ------------- | ------------- |
| [S3B](./RawFigures/SuppFig_3_B.pdf) | `python Plotting.py` |

## `Sketch/Orientation`

| Figure(s) | Command for creating figure |
| ------------- | ------------- |
| [S4A_i](./RawFigures/SuppFig_5_Ai.eps),[S4A_ii](./RawFigures/SuppFig_5_Aii.eps) | `python Plotting.py` |
