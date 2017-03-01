
/*
 * Global structures for MPI communication
 */

#pragma once

typedef struct mpistate {

  /* Level of thread support provided by MPI implementation */
  int thread_support;

  /* Process rank */
  int rank;

  /* Number of ranks in MPI_COMM_WORLD */
  int nprocs;

} mpistate;

