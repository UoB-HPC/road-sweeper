/*
 * This file is part of road-sweeper.
 *
 * road-sweeper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * road-sweeper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with road-sweeper.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "comms.h"
#include "compute.h"
#include <mpi.h>
#include "options.h"
#include <stdlib.h>
#include "sweep.h"

void init_serial_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf);
void end_serial_sweep(double *ybuf, double *zbuf);

/* Perform a vanilla KBA sweep without using OpenMP threads */
timings serial_sweep(mpistate mpi, options opt) {

  timings time = {
    .sweeping = 0.0,
    .setup = 0.0,
    .comms = 0.0
  };

  /* Message buffers */
  time.setup = MPI_Wtime();
  const int ycount = opt.nang * opt.nz * opt.chunklen;
  const int zcount = opt.nang * opt.ny * opt.chunklen;
  double *ybuf;
  double *zbuf;
  init_serial_sweep(ycount, zcount, &ybuf, &zbuf);
  time.setup = MPI_Wtime() - time.setup;

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
            double comtime = MPI_Wtime();
            if (j == 0) {
              MPI_Recv(ybuf, ycount, MPI_DOUBLE, mpi.yhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(ybuf, ycount, MPI_DOUBLE, mpi.ylo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (k == 0) {
              MPI_Recv(zbuf, zcount, MPI_DOUBLE, mpi.zhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(zbuf, zcount, MPI_DOUBLE, mpi.zlo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            time.comms += MPI_Wtime() - comtime;

            /* Do proportional "work" */
            for (int w = 0; w < opt.nang*opt.chunklen*opt.ny*opt.nz; w++) {
              compute();
            }

            /* Send payload to downwind neighbours */
            comtime = MPI_Wtime();
            MPI_Waitall(2, req, MPI_STATUS_IGNORE);

            if (j == 0) {
              MPI_Isend(ybuf, ycount, MPI_DOUBLE, mpi.ylo, 0, MPI_COMM_WORLD, req+0);
            }
            else {
              MPI_Isend(ybuf, ycount, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req+0);
            }

            if (k == 0) {
              MPI_Isend(zbuf, zcount, MPI_DOUBLE, mpi.zlo, 0, MPI_COMM_WORLD, req+1);
            }
            else {
              MPI_Isend(zbuf, zcount, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req+1);
            }
            time.comms += MPI_Wtime() - comtime;

          } /* End nchunks loop */
        } /* End ng loop */
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */

  /* End the timer */
  double tock = MPI_Wtime();

  time.sweeping = tock-tick;

  end_serial_sweep(ybuf, zbuf);

  time.setup += MPI_Wtime() - tock;

  return time;
}

/* Allocate MPI message buffers */
void init_serial_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf) {
  (*ybuf) = malloc(sizeof(double)*ycount);
  (*zbuf) = malloc(sizeof(double)*zcount);
}

/* Free MPI message buffers */
void end_serial_sweep(double *ybuf, double *zbuf) {
  free(ybuf);
  free(zbuf);
}

