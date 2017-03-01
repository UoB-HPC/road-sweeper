
#include <mpi.h>

#define LOAD 0.1

void compute(void) {
  double start = MPI_Wtime();
  while (MPI_Wtime() - start < LOAD) {}
}

