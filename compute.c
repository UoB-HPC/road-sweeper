
#include <mpi.h>

void compute(const double load) {
  /* Spin until enough time has passed */
  double start = MPI_Wtime();
  while (MPI_Wtime() - start < load) {}
}

