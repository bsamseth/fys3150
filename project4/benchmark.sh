#!/bin/bash
# Note, if the script fails, a re-run has been effective in solving the problems
# This script will create the following structure (output of tree)
# .
# ├── benchmarks
# │   ├── output
# │   │   ├── gauss_laguerre_25.dat
# │   │   ├── gauss_legendre_25.dat
# │   │   ├── smartMonteCarlo_10000000.dat
# │   │   └── uniformmontecarlo_10000000.dat
# │   └── fig
# │       └── psiPlot.png
#
# These results should correspond to the data and figures referenced
# in the report, aswell as the benchmarks/

# define name of folders
PARENT=$(pwd)
DIR=$PARENT/benchmarks
SRC=src
DATA=output
FIG=fig
DIRSRC=$PARENT/$SRC
DIRDATA=$DIR/$DATA
DIRFIG=$DIR/$FIG
# create data and fig folder
echo "Creating folder structure"
mkdir -p $DIRDATA
mkdir -p $DIRFIG
mkdir -p $DIR/$SRC/data

cp $DIRSRC/*.cpp $DIR/$SRC
cp $DIRSRC/*.h $DIR/$SRC
cp $DIRSRC/*.py $DIR/$SRC
cp $DIRSRC/Makefile $DIR/$SRC
# enter $DIR, needed because some programs write to relative path
echo "Entering benchmark directory"
cd $DIR/$SRC

# compile all executables
echo "Compiling source files"
make
echo "Compilation complete"

# =================================
# run programs with selcted output.
# =================================
# This should reproduce the results reffered to in the report.
# Please note that this will take a long time to run
# you may change the parameters to get less acurate
# results in a shorter time
echo "Running MPIising.cpp, computing data refrenced in table 2 and 3"
mpirun -n 4 MPIising.x data/out 2 1000000 1 2.3 1.3 0
echo "Running MPIisingImproved.cpp, computing data refrenced i table 3"
mpirun -n 4 MPIisingImproved.x data/out 2 1000000 2.3 2.3 1 0

# this will genereate all the figures shown in the report
# you will need to close the figures in order to continue calculations. 
echo "Running expectationValues.py for ordered initial config. (will take some time)"
python expectationValues.py 20 1 2.4 1.4 1000 500000 35 0 
echo "Running expectationValues.py for random initial config. (will take some time)"
python expectationValues.py 20 1 2.4 1.4 1000 500000 35 1 
echo "Running probabilities.py"
python probabilities.py 20 1 2.4 0.7 1000000 0 
echo "Running critical.py (will take a long time)"
python critical.py [20,40,60,80,100] 2.2 2.35 0.01 1000000 0

echo "Cleaning up"
cp data/* ../$DATA
cd ..
rm -r $SRC/

# exit back to $PARENT
echo "Exiting benchmark directory"
cd $PARENT

