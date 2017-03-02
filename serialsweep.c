
#include "compute.h"
#include <mpi.h>
#include "serialsweep.h"

/* Perform a vanilla KBA sweep without using OpenMP threads */
void serial_sweep(mpistate mpi, options opt, messages msg) {

  const int count = opt.pencil*opt.chunklen;

  /* Send requests */
  MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  /* Octant loops - 0 is stepping backwards, 1 is stepping forwards */
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {

        /* Loop over messages to send per octant */
        for (int c = 0; c < opt.nchunks; c++) {

          /* Receive payload from upwind neighbours */
          if (j == 0) {
            MPI_Recv(msg.ybuf, count, MPI_DOUBLE, mpi.yhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          else {
            MPI_Recv(msg.ybuf, count, MPI_DOUBLE, mpi.ylo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }

          if (k == 0) {
            MPI_Recv(msg.zbuf, count, MPI_DOUBLE, mpi.zhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }
          else {
            MPI_Recv(msg.zbuf, count, MPI_DOUBLE, mpi.zlo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          }

          /* Do "work" */
          compute(opt.work);

          /* Send payload to downwind neighbours */
          MPI_Waitall(2, req, MPI_STATUS_IGNORE);

          if (j == 0) {
            MPI_Isend(msg.ybuf, count, MPI_DOUBLE, mpi.ylo, 0, MPI_COMM_WORLD, req+0);
          }
          else {
            MPI_Isend(msg.ybuf, count, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req+0);
          }

          if (k == 0) {
            MPI_Isend(msg.zbuf, count, MPI_DOUBLE, mpi.zlo, 0, MPI_COMM_WORLD, req+1);
          }
          else {
            MPI_Isend(msg.zbuf, count, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req+1);
          }

        } /* End nchunks loop */
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */
}

