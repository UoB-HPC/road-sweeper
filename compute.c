
#include <mpi.h>

#ifndef LOAD
#define LOAD 0.1
#endif

void compute(void) {
  /* Spin until enough time has passed */
  double start = MPI_Wtime();
  while (MPI_Wtime() - start < LOAD) {}
}

