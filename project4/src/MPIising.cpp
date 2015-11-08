/* 
   Program to solve the two-dimensional Ising model 
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <sstream>
#include <string>

#include "functions.h"

using namespace  std;


// output file
ofstream ofile;

// Main program begins here

int main(int argc, char* argv[])
{
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs, my_rank, numprocs;
  double w[17], average[5], total_average[5], 
    initial_temp, final_temp, E, M, temp_step;
  int randomizer;
  int num_configurations = 0;
  

  if (argc < 8) {
    cout << "Usage: " << argv[0] << " outputFile n_spins MCCycles T0 T1 dT randomizer" << endl;
    exit(1);
  }
    
  n_spins = atoi(argv[2]); mcs = atoi(argv[3]);
  initial_temp = atof(argv[4]); final_temp = atof(argv[5]);
  temp_step = atof(argv[6]); randomizer = atoi(argv[7]);
  
  //  MPI initializations
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0 && argc <= 1) {
    cout << "Bad Usage: " << argv[0] << 
      " read outputfile n_spins MCs T0 T1 dT rand" << endl
	 << "rand is either zero or one. If zero initial spin" << endl 
	 << "configuration is all down. If =one initial spin is random." 
	 << endl;
    exit(1);
  }
  if (my_rank == 0) {
    outfilename=argv[1];
    stringstream sstm;
    sstm << outfilename << "_" << n_spins << "_" << mcs
	 << "_" << initial_temp << "_" << final_temp
	 << "_" << temp_step << "_" << randomizer << ".dat";
    string filename = sstm.str();
    ofile.open(filename.c_str());
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << setw(15) << "T";
      ofile << setw(15) << "E";
      ofile << setw(15) << "Evar/T^2";
      ofile << setw(15) << "M";
      ofile << setw(15) << "Mvar/T";
      ofile << setw(15) << "MvarAbs/T";
      ofile << setw(15) << "|M|" ;
      ofile << setw(15) << "Num_config" << endl;
  }

  /*
  Determine number of intervall which are used by all processes
  myloop_begin gives the starting point on process my_rank
  myloop_end gives the end point for summation on process my_rank
  */
  int no_intervalls = mcs/numprocs;
  int myloop_begin = my_rank*no_intervalls + 1;
  int myloop_end = (my_rank+1)*no_intervalls;
  if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

  // broadcast to all nodes common variables
  MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  //  Allocate memory for spin matrix
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
  // every node has its own seed for the random numbers, this is important else
  // if one starts with the same seed, one ends with the same random numbers
  idum = -1-my_rank; //-time(NULL)-my_rank;  // random starting point
  // Start Monte Carlo sampling by looping over T first
  for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){
    //    initialise energy and magnetization 
    E = M = 0.;
    // initialise array for expectation values
    initialize(n_spins, spin_matrix, E, M, randomizer, idum);
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 5; i++) total_average[i] = 0.;
    // start Monte Carlo computation
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
      Metropolis(n_spins, idum, spin_matrix, E, M, w, &num_configurations);
      // update expectation values  for local node
      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);
    }
    // Find total average
    for( int i =0; i < 5; i++){
      MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    // print results
    if ( my_rank == 0) {
      ofile << output(n_spins, mcs, temperature, total_average, &num_configurations);
    }
  }
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file
  // End MPI
  MPI_Finalize (); 
  return 0;
}

