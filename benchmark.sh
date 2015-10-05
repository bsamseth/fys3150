#!/bin/bash
# Note, if the script fails, a re-run has been effective in solving the problems
# This script will create the following structure (output of tree)
# .
# ├── benchmarks
# │   ├── data
# │   │   ├── lower_3_states_nstep_55_rho_3_5.dat
# │   │   ├── lower_3_states_nstep_85_rho_4.dat
# │   │   └── lower_3_states_nstep_90_rho_4_5.dat
# │   ├── fig
# │   │   ├── number_of_transformations_trend_plot.jpg
# │   │   └── theta.png
# │   └── unitTest_result.txt
#
# These results should correspond to the data and figures referenced
# in the report.


# define name of folders
PARENT=$(pwd)
DIR=benchmarks
DATA=data
FIG=fig
DIRDATA=$DIR/$DATA
DIRFIG=$DIR/$FIG
# create data and fig folder
echo "Creating folder structure"
mkdir -p $DIRDATA
mkdir -p $DIRFIG

cp $PARENT/*.cpp $DIR/
cp $PARENT/*.h $DIR
cp $PARENT/*.py $DIR
cp $PARENT/Makefile $DIR
# enter $DIR, needed because some programs write to relative path
echo "Entering benchmark directory"
cd $DIR

# compile all executables
echo "Compiling source files"
make all
echo "Compilation complete"

# =================================
# run programs with selcted output.
# =================================
#This should reproduce the results given in section 4
echo "Running solve_lower_3_states.cpp for selected arguments"
./solve_lower_3_states.x 55 3.5 > $DATA/lower_3_states_nstep_55_rho_3_5.dat
./solve_lower_3_states.x 55 3.5 > $DATA/lower_3_states_nstep_85_rho_4.dat
./solve_lower_3_states.x 55 3.5 > $DATA/lower_3_states_nstep_90_rho_4_5.dat

# this runs all tests. Output should report all completed without problem
echo "Running unitTest.cpp "
./unitTest.x > unitTest_result.txt

# this will genereate all the figures shown in the report
# they will aprear on screen one by one. When exited, they are
# saved in $DIRFIG. This does not include plot of eigenstates, because
# this plot should be adjusted in width and height before saved. Option
# to save is presented in plot GUI.
echo "Running tangensPlot.py"
python tangensPlot.py
echo "Running number_of_transformations.py"
python number_of_transformations.py 5 10 15 25 40 60 80 100 120
echo "Running plotEigenstates.py"
python plotEigenstates.py 100 2.25 5 100 5 1 100 7 0.5 100 50 0.01

echo "Cleaning up"
rm *.cpp
rm *.h
rm *.py
rm Makefile
rm *.x
rm *.o
rm $DATA/*interacting*.dat

# exit back to $PARENT
echo "Exiting benchmark directory"
cd $PARENT
