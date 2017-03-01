
/*
 * Global structures for MPI communication
 */

#pragma once

#include "options.h"

typedef struct mpistate {

  /* Level of thread support provided by MPI implementation */
  int thread_support;

  /* Process rank */
  int rank;

  /* Number of ranks in MPI_COMM_WORLD */
  int nprocs;

  /* Number of ranks in 2D decomposition */
  int npey;
  int npez;

  /* Process rank in 2D decomposition */
  int y;
  int z;

  /* Neighbour ranks - 2D decomposition of 3D domain */
  int yhi;
  int ylo;
  int zhi;
  int zlo;

} mpistate;

void decompose(mpistate *mpi, options opt);

