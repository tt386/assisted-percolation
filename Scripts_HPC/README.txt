ISCA Specific:
To then run python:
module load Anaconda3

To check the squeue:
squeue --me




This will outline the uses of all the  files here!

Periphery:

Params.py:          This contains parameters endemic across all simulations we want to carry out, such as system size, end conditions, obstacle sizes, etc.

CreateParams.py:    This creates ParamPairs.npz (so names because usually there are two such parameters in each, for a 2d plot) from a generated list of parameter values: a list of parameter pairs values. The directory is named here too, so it takes input from Params.py
    
ParamPairs.npz:     This contains all combinations of parameters for use in the simulation: each repeat of the simulation will use a different row in the list of parameter pairs.

ContourPlotting.py:     This plots the data from simulations from an input directory within SaveFiles




################
###PROCEDURE####
################

1) Create the ParamPairs with CreateParamPairs.py. This also creates the directory that will save all of the results of the simulation, within SaveFiles. Make a note of the total number of parameter pairs that are produced: let's call this 'n'

2) Run the Parallel_Bash script:
    sbatch --array=0-n Parallel_Bash.sh





