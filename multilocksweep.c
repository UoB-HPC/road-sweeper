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
#include <omp.h>
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include "sweep.h"

void init_par_mpi_multi_lock_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf);
void end_par_mpi_multi_lock_sweep(double *ybuf, double *zbuf);

/* Perform a KBA sweep using OpenMP threads for concurrent group sweeps */
timings par_mpi_multi_lock_sweep(mpistate mpi, options opt) {

  timings time = {
    .sweeping = 0.0,
    .setup = 0.0,
    .comms = 0.0
  };

  time.setup = MPI_Wtime();

  /* Check MPI threading model is high enough */
  if (mpi.thread_support < MPI_THREAD_SERIALIZED) {
    if (mpi.rank == 0) {
      printf("MPI library must support MPI_THREAD_SERIALIZED\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  /*
   * Need to use OpenMP locks if we are only MPI_THREAD_SERIALIZED.
   * We must ensure only one thread calls at once, and we cannot
   * use the if clause on a critial region, so we implement this
   * manually with locks.
   * Only one lock should be required
   */
  int nthrds;
  #pragma omp parallel
  {
    nthrds = omp_get_num_threads();
  } 

  omp_lock_t *lock = malloc(sizeof(omp_lock_t)*nthrds);
  if (mpi.thread_support == MPI_THREAD_SERIALIZED) {
    for (int l = 0; l < nthrds; l++)
      omp_init_lock(lock+l);
  }

  /* Message buffers */
  const int ycount = opt.nang * opt.nz * opt.chunklen;
  const int zcount = opt.nang * opt.ny * opt.chunklen;
  double *ybuf;
  double *zbuf;
  init_par_mpi_multi_lock_sweep(opt.ng*ycount, opt.ng*zcount, &ybuf, &zbuf);
  time.setup = MPI_Wtime() - time.setup;

  /* Start the timer */
  double tick = MPI_Wtime();

  /* Start parallel region. Locks must be set within this region for spec conformance */
#pragma omp parallel
{
  /* Send requests - 2 per thread */
  MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  const int thrd = omp_get_thread_num();

  /* Lock all the locks, except the first */
  if (thrd > 0 && mpi.thread_support == MPI_THREAD_SERIALIZED) {
    omp_set_lock(lock+thrd);
  }


  /* Octant loops - 0 is stepping backwards, 1 is stepping forwards */
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {

        /* Loop over energy groups in parallel, setting up
         * one concurrent sweep per group
         */
        #pragma omp for
        for (int g = 0; g < opt.ng; g++) {


          /* Loop over messages to send per octant */
          for (int c = 0; c < opt.nchunks; c++) {

            /* Receive payload from upwind neighbours */
            double comtime = MPI_Wtime();

            /* Lock if necessary before comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_set_lock(lock+thrd);
            }

            if (j == 0) {
              MPI_Recv(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.yhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.ylo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (k == 0) {
              MPI_Recv(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zhi, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
              MPI_Recv(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zlo, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            /* Just time last thread */
            if (thrd == nthrds-1) {
              time.comms += MPI_Wtime() - comtime;
            }

            /* Unlock neighbour thread if necessary after comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              int nxt = (thrd+1 == nthrds) ? 0 : thrd+1;
              omp_unset_lock(lock+nxt);
            }

            /* Do proportional "work" */
            for (int w = 0; w < opt.nang*opt.chunklen*opt.ny*opt.nz; w++) {
              compute();
            }

            /* Send payload to downwind neighbours */
            comtime = MPI_Wtime();

            /* Lock if necessary before comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_set_lock(lock+thrd);
            }

            MPI_Waitall(2, req, MPI_STATUS_IGNORE);

            if (j == 0) {
              MPI_Isend(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.ylo, 0, MPI_COMM_WORLD, req+0);
            }
            else {
              MPI_Isend(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req+0);
            }

            if (k == 0) {
              MPI_Isend(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zlo, 0, MPI_COMM_WORLD, req+1);
            }
            else {
              MPI_Isend(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req+1);
            }

            /* Just time last thread */
            if (thrd == nthrds-1) {
              time.comms += MPI_Wtime() - comtime;
            }

            /* Unlock next thread if necessary after comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              int nxt = (thrd+1 == nthrds) ? 0 : thrd+1;
              omp_unset_lock(lock+nxt);
            }

          } /* End nchunks loop */
        } /* End ng loop */
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */

  /* Unset all the locks, except the first */
  if (thrd > 0 && mpi.thread_support == MPI_THREAD_SERIALIZED) {
    omp_unset_lock(lock+thrd);
  }

} /* End parallel region */

  /* End the timer */
  double tock = MPI_Wtime();

  time.sweeping = tock-tick;

  /* Destroy lock */
  if (mpi.thread_support == MPI_THREAD_SERIALIZED) {
    for (int l = 0; l < nthrds; l++) {
      omp_destroy_lock(lock+l);
    }
  }
  end_par_mpi_multi_lock_sweep(ybuf, zbuf);

  time.setup += MPI_Wtime() - tock;

  return time;
}

/* Allocate MPI message buffers */
void init_par_mpi_multi_lock_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf) {
  (*ybuf) = malloc(sizeof(double)*ycount);
  (*zbuf) = malloc(sizeof(double)*zcount);
}

/* Free MPI message buffers */
void end_par_mpi_multi_lock_sweep(double *ybuf, double *zbuf) {
  free(ybuf);
  free(zbuf);
}

