#!/bin/bash

# define name of benchmark folder
DIR=data
# create data folder
mkdir -p $DIR

# compile all executables
make all

# run programs with selcted output. This should reproduce the results
# given in the benchmark folder in the master branch of the Github-repo.
./solve_lower_3_states.x 55 3.5 > $DIR/lower_3_states_nstep_55_rho_3_5.dat
./solve_lower_3_states.x 55 3.5 > $DIR/lower_3_states_nstep_85_rho_4.dat
./solve_lower_3_states.x 55 3.5 > $DIR/lower_3_states_nstep_90_rho_4_5.dat

./unitTest.x > $DIR/unitTest_result.txt

./interactingElectrons.x 1 50 5 1
./interactingElectrons.x 0 50 5 1

mv $DIR benchmarks
