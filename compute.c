
#include <mpi.h>

#ifndef LOAD
#define LOAD 0.1
#endif

void compute(void) {
  double start = MPI_Wtime();
  while (MPI_Wtime() - start < LOAD) {}
}

