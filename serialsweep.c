
#include <mpi.h>
#include "serialsweep.h"

/* Perform a vanilla KBA sweep without using OpenMP threads */
void serial_sweep(mpistate mpi, options opt) {

  double * restrict ybuf = NULL;
  double * restrict zbuf = NULL;
  int ycount = 0;
  int zcount = 0;

  /* Send requests */
  MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  /* Receive payload from upwind neighbours */
  MPI_Recv(ybuf, ycount, MPI_DOUBLE, mpi.ylo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(zbuf, zcount, MPI_DOUBLE, mpi.zlo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  printf("Rank %d starting\n", mpi.rank);
  sleep(5);

  /* Send payload to downwind neighbours */
  MPI_Waitall(2, req, MPI_STATUS_IGNORE);
  MPI_Isend(ybuf, ycount, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req+0);
  MPI_Isend(zbuf, zcount, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req+1);
}

