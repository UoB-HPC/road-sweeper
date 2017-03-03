
#include "comms.h"
#include "compute.h"
#include <mpi.h>
#include <omp.h>
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include "sweep.h"

void init_par_mpi_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf);
void end_par_mpi_sweep(double *ybuf, double *zbuf);

/* Perform a KBA sweep using OpenMP threads for concurrent group sweeps */
timings par_mpi_sweep(mpistate mpi, options opt) {

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
  omp_lock_t lock;
  if (mpi.thread_support == MPI_THREAD_SERIALIZED) {
    omp_init_lock(&lock);
  }

  /* Message buffers */
  const int ycount = opt.nang * opt.nz * opt.chunklen;
  const int zcount = opt.nang * opt.ny * opt.chunklen;
  double *ybuf;
  double *zbuf;
  init_par_mpi_sweep(opt.ng*ycount, opt.ng*zcount, &ybuf, &zbuf);
  time.setup = MPI_Wtime() - time.setup;

  /* Send requests - 2 per thread */
  int nthrds;
  #pragma omp parallel
  {
    nthrds = omp_get_num_threads();
  }
  MPI_Request req[nthrds][2];
#pragma omp parallel
  {
    req[omp_get_thread_num()][0] = MPI_REQUEST_NULL;
    req[omp_get_thread_num()][1] = MPI_REQUEST_NULL;
  }

  /* Start the timer */
  double tick = MPI_Wtime();

  /* Octant loops - 0 is stepping backwards, 1 is stepping forwards */
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {

        /* Loop over energy groups in parallel, setting up
         * one concurrent sweep per group
         */
        #pragma omp parallel for
        for (int g = 0; g < opt.ng; g++) {

          const int thrd = omp_get_thread_num();

          /* Loop over messages to send per octant */
          for (int c = 0; c < opt.nchunks; c++) {

            /* Receive payload from upwind neighbours */
            double comtime = MPI_Wtime();

            /* Lock if necessary before comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_set_lock(&lock);
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

            /* Unlock if necessary after comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_unset_lock(&lock);
            }

            /* Do proportional "work" */
            for (int w = 0; w < opt.nang*opt.chunklen*opt.ny*opt.nz; w++) {
              compute();
            }

            /* Send payload to downwind neighbours */
            comtime = MPI_Wtime();

            /* Lock if necessary before comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_set_lock(&lock);
            }

            MPI_Waitall(2, req[thrd], MPI_STATUS_IGNORE);

            if (j == 0) {
              MPI_Isend(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.ylo, 0, MPI_COMM_WORLD, req[thrd]+0);
            }
            else {
              MPI_Isend(ybuf+g*ycount, ycount, MPI_DOUBLE, mpi.yhi, 0, MPI_COMM_WORLD, req[thrd]+0);
            }

            if (k == 0) {
              MPI_Isend(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zlo, 0, MPI_COMM_WORLD, req[thrd]+1);
            }
            else {
              MPI_Isend(zbuf+g*zcount, zcount, MPI_DOUBLE, mpi.zhi, 0, MPI_COMM_WORLD, req[thrd]+1);
            }

            /* Just time last thread */
            if (thrd == nthrds-1) {
              time.comms += MPI_Wtime() - comtime;
            }

            /* Unlock if necessary after comms */
            if(mpi.thread_support == MPI_THREAD_SERIALIZED) {
              omp_unset_lock(&lock);
            }

          } /* End nchunks loop */
        } /* End ng loop */
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */

  /* End the timer */
  double tock = MPI_Wtime();

  time.sweeping = tock-tick;

  /* Destroy lock */
  if (mpi.thread_support == MPI_THREAD_SERIALIZED) {
    omp_destroy_lock(&lock);
  }
  end_par_mpi_sweep(ybuf, zbuf);

  time.setup += MPI_Wtime() - tock;

  return time;
}

/* Allocate MPI message buffers */
void init_par_mpi_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf) {
  (*ybuf) = malloc(sizeof(double)*ycount);
  (*zbuf) = malloc(sizeof(double)*zcount);
}

/* Free MPI message buffers */
void end_par_mpi_sweep(double *ybuf, double *zbuf) {
  free(ybuf);
  free(zbuf);
}

