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

cp $DIRSRC/*.cpp $DIR/
cp $DIRSRC/*.h $DIR
cp $DIRSRC/*.py $DIR
cp $DIRSRC/Makefile $DIR
# enter $DIR, needed because some programs write to relative path
echo "Entering benchmark directory"
cd $DIR

# compile all executables
echo "Compiling source files"
make
echo "Compilation complete"

# =================================
# run programs with selcted output.
# =================================
#This should reproduce the results reffered to in the report.
echo "Running gauss_legendre.x for selected arguments"
./gauss_legendre.x -2 2 31 > $DATA/gauss_legendre.dat
echo "Running gauss_laguerre.x for selected arguments"
./gauss_laguerre.x 40 > $DATA/gauss_laguerre.dat
echo "Running uniformmonteCarlo.x for selected arguments"
./uniformmonteCarlo.x -2 2 10000000 > $DATA/uniformmontecarlo.dat
echo "Running smartMonteCarlo.x for selected arguments"
./smartMonteCarlo.x 10000000 > $DATA/smartMonteCarlo.dat

# this will genereate all the figures shown in the report (only one)
echo "Running plotFunc.py"
python plotFunc.py
# plotFunc.py stores the produced fig in ../fig/psiPlot.png
# so cp this benchmarks/fig/
cp $PARENT/fig/psiPlot.png $DIRFIG/

echo "Cleaning up"
rm *.cpp
rm *.h
rm *.py
rm Makefile
rm *.x
rm *.o

# exit back to $PARENT
echo "Exiting benchmark directory"
cd $PARENT
