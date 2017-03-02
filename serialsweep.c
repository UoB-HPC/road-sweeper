
#include "comms.h"
#include "compute.h"
#include <mpi.h>
#include "options.h"
#include <stdlib.h>

void init_serial_sweep(int count, double **ybuf, double **zbuf);
void end_serial_sweep(double *ybuf, double *zbuf);

/* Perform a vanilla KBA sweep without using OpenMP threads */
double serial_sweep(mpistate mpi, options opt) {

  /* Message buffers */
  int count = opt.nang * opt.pencil * opt.chunklen;
  double *ybuf;
  double *zbuf;
  init_serial_sweep(count, &ybuf, &zbuf);

  /* Send requests */
  MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  /* Start the timer */
  double tick = MPI_Wtime();

  /* Octant loops - 0 is stepping backwards, 1 is stepping forwards */
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {

        /* Loop over energy groups in serial */
        for (int g = 0; g < opt.ng; g++) {

          /* Loop over messages to send per octant */
          for (int c = 0; c < opt.nchunks; c++) {

            /* Receive payload from upwind neighbours */
            if (j == 0) {
              MPI_Recv(ybuf, count, MPI_DOUBLE, mpi.yhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(ybuf, count, MPI_DOUBLE, mpi.ylo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (k == 0) {
              MPI_Recv(zbuf, count, MPI_DOUBLE, mpi.zhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(zbuf, count, MPI_DOUBLE, mpi.zlo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            /* Do "work" */
            compute(opt.work);

            /* Send payload to downwind neighbours */
            MPI_Waitall(2, req, MPI_STATUS_IGNORE);

            if (j == 0) {
              MPI_Isend(ybuf, count, MPI_DOUBLE, mpi.ylo, 0, MPI_COMM_WORLD, req+0);
            }
            else {
              MPI_Isend(ybuf, count, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req+0);
            }

            if (k == 0) {
              MPI_Isend(zbuf, count, MPI_DOUBLE, mpi.zlo, 0, MPI_COMM_WORLD, req+1);
            }
            else {
              MPI_Isend(zbuf, count, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req+1);
            }

          } /* End nchunks loop */
        } /* End ng loop */
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */

  /* End the timer */
  double tock = MPI_Wtime();

  end_serial_sweep(ybuf, zbuf);

  return tock-tick;
}

/* Allocate MPI message buffers */
void init_serial_sweep(int count, double **ybuf, double **zbuf) {
  (*ybuf) = malloc(sizeof(double)*count);
  (*zbuf) = malloc(sizeof(double)*count);
}

/* Free MPI message buffers */
void end_serial_sweep(double *ybuf, double *zbuf) {
  free(ybuf);
  free(zbuf);
}

