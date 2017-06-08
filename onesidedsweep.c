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

void init_one_sided_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf, MPI_Win *ywin, MPI_Win *zwin, const int ylo, const int yhi, const int zlo, const int zhi);
void end_one_sided_sweep(MPI_Win *ywin, MPI_Win *zwin, const int ylo, const int yhi, const int zlo, const int zhi);

/* Perform a KBA sweep threading over groups inside the chunk
 *
 * This method uses the one-sided RMA MPI communication pattern,
 * with _passive_ synchronisation.
 * A window is created on each MPI rank.
 * Shared locks created on this window for neighbouring processors.
 * The sweep goes in all directions, and so therefore this lock is shared
 * rather than exclusive as each rank has two neighbours which may write
 * to the buffer dependent on the sweep direction.
 * Signals are used to synchronise between ranks.
 * The send protocol is therefore:
 *   1. Poll for safe signal
 *   2. MPI_Put payload + MPI_Flush
 *   3. Put done signal + MPI_Flush
 * The receive protocol is therefore:
 *   1. Put safe signal + MPI_Flush
 *   2. Poll for done signal
 *   3. Use data
 */

#define SAFE_SIGNAL 123456789.0
#define SENT_SIGNAL 987654321.0
#define NULL_SIGNAL -1.0
/* Offsets from the main array where the signals can be found */
#define SAFE_OFFSET 0
#define SENT_OFFSET 1

timings one_sided_sweep(mpistate mpi, options opt) {

  timings time = {
    .sweeping = 0.0,
    .setup = 0.0,
    .comms = 0.0
  };

  /* Message buffers */
  time.setup = MPI_Wtime();
  const int ycount = opt.nang * opt.nz * opt.chunklen * opt.ng;
  const int zcount = opt.nang * opt.ny * opt.chunklen * opt.ng;
  double *ybuf;
  double *zbuf;
  MPI_Win ywin, zwin;
  init_one_sided_sweep(ycount, zcount, &ybuf, &zbuf, &ywin, &zwin, mpi.ylo, mpi.yhi, mpi.zlo, mpi.zhi);
  time.setup = MPI_Wtime() - time.setup;

  /* Send requests */
  MPI_Request req[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  /* Start the timer */
  double tick = MPI_Wtime();

  /* Octant loops - 0 is stepping backwards, 1 is stepping forwards */
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 2; i++) {
        printf("%d: starting oct: %d\n", mpi.rank, i+2*j+4*k);

        /* Loop over messages to send per octant */
        for (int c = 0; c < opt.nchunks; c++) {

          /* Receive payload from upwind neighbours */
          double comtime = MPI_Wtime();
          if (j == 0) {
            /* Do comms if internal boundary */
            if (mpi.yhi != MPI_PROC_NULL) {
              printf("%d: recv j=0\n", mpi.rank);
              /* Send safe signal */
              ybuf[ycount+SAFE_OFFSET] = SAFE_SIGNAL;
              MPI_Put(ybuf+ycount+SAFE_OFFSET, 1, MPI_DOUBLE, mpi.yhi, ycount+SAFE_OFFSET, 1, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.yhi, ywin);

              /* Poll for sent signal - lock required around access */
              int sent = 0;
              while (!sent) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, ywin);
                if (ybuf[ycount+SENT_OFFSET] == SENT_SIGNAL)
                  sent = 1;
                MPI_Win_unlock(mpi.rank, ywin);
              }

              /* Reset signal */
              ybuf[ycount+SENT_OFFSET] = NULL_SIGNAL;
            }
          }
          else {
            /* Do comms if internal boundary */
            if (mpi.ylo != MPI_PROC_NULL) {
              printf("%d: recv j/=0\n", mpi.rank);
              /* Send safe signal */
              ybuf[ycount+SAFE_OFFSET] = SAFE_SIGNAL;
              MPI_Put(ybuf+ycount+SAFE_OFFSET, 1, MPI_DOUBLE, mpi.ylo, ycount+SAFE_OFFSET, 1, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.ylo, ywin);

              /* Poll for sent signal - lock required around access */
              int sent = 0;
              while (!sent) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, ywin);
                if (ybuf[ycount+SENT_OFFSET] == SENT_SIGNAL)
                  sent = 1;
                MPI_Win_unlock(mpi.rank, ywin);
              }

              /* Reset signal */
              ybuf[ycount+SENT_OFFSET] = NULL_SIGNAL;
            }
          }
          if (k == 0) {
            /* Do comms if internal boundary */
            if (mpi.zhi != MPI_PROC_NULL) {
              printf("%d: recv k=0\n", mpi.rank);
              /* Send safe signal */
              zbuf[zcount+SAFE_OFFSET] = SAFE_SIGNAL;
              MPI_Put(zbuf+zcount+SAFE_OFFSET, 1, MPI_DOUBLE, mpi.zhi, zcount+SAFE_OFFSET, 1, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zhi, zwin);

              /* Poll for sent signal - lock required around access */
              int sent = 0;
              while (!sent) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, zwin);
                if (zbuf[zcount+SENT_OFFSET] == SENT_SIGNAL)
                  sent = 1;
                MPI_Win_unlock(mpi.rank, zwin);
              }

              /* Reset signal */
              zbuf[zcount+SENT_OFFSET] = NULL_SIGNAL;
            }
          }
          else {
            /* Do comms if internal boundary */
            if (mpi.zlo != MPI_PROC_NULL) {
              printf("%d: recv k/=0\n", mpi.rank);
              /* Send safe signal */
              zbuf[zcount+SAFE_OFFSET] = SAFE_SIGNAL;
              MPI_Put(zbuf+zcount+SAFE_OFFSET, 1, MPI_DOUBLE, mpi.zlo, zcount+SAFE_OFFSET, 1, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zlo, zwin);

              /* Poll for sent signal - lock required around access */
              int sent = 0;
              while (!sent) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, zwin);
                if (zbuf[zcount+SENT_OFFSET] == SENT_SIGNAL)
                  sent = 1;
                MPI_Win_unlock(mpi.rank, zwin);
              }

              /* Reset signal */
              zbuf[zcount+SENT_OFFSET] = NULL_SIGNAL;
            }
          }
          time.comms += MPI_Wtime() - comtime;


          /*********************************************************************
          * Compute
          *********************************************************************/
          printf("%d: compute oct %d\n", mpi.rank, i+2*j+4*k);
          #pragma omp parallel for
          for (int g = 0; g < opt.ng; g++) {

            /* Do proportional "work" */
            for (int w = 0; w < opt.nang*opt.chunklen*opt.ny*opt.nz; w++) {
              compute();
            }

          } /* End group loop */

          /*********************************************************************
          * End compute
          *********************************************************************/


          /* Put (send) payload in downwind neighbours window */
          comtime = MPI_Wtime();

          if (j == 0) {
            /* Do comms if internal boundary */
            if (mpi.ylo != MPI_PROC_NULL) {
              /* Poll for safe to send signal - lock required around access */
              printf("%d: send j=0\n", mpi.rank);
              int safe = 0;
              while (!safe) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, ywin);
                if (ybuf[ycount+SAFE_OFFSET] == SAFE_SIGNAL)
                  safe = 1;
                MPI_Win_unlock(mpi.rank, ywin);
              }

              /* Reset signal */
              ybuf[ycount+SAFE_OFFSET] = NULL_SIGNAL;

              /* Put payload */
              MPI_Put(ybuf, ycount, MPI_DOUBLE, mpi.ylo, 0, ycount, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.ylo, ywin);

              /* Send sent signal */
              ybuf[ycount+SENT_OFFSET] = SENT_SIGNAL;
              MPI_Put(ybuf+ycount+SENT_OFFSET, 1, MPI_DOUBLE, mpi.ylo, ycount+SENT_OFFSET, 1, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.ylo, ywin);
            }
          }
          else {
            /* Do comms if internal boundary */
            if (mpi.yhi != MPI_PROC_NULL) {
              printf("%d: send j/=0\n", mpi.rank);
              /* Poll for safe to send signal - lock required around access */
              int safe = 0;
              while (!safe) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, ywin);
                if (ybuf[ycount+SAFE_OFFSET] == SAFE_SIGNAL)
                  safe = 1;
                MPI_Win_unlock(mpi.rank, ywin);
              }

              /* Reset signal */
              ybuf[ycount+SAFE_OFFSET] = NULL_SIGNAL;

              /* Put payload */
              MPI_Put(ybuf, ycount, MPI_DOUBLE, mpi.yhi, 0, ycount, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.yhi, ywin);

              /* Send sent signal */
              ybuf[ycount+SENT_OFFSET] = SENT_SIGNAL;
              MPI_Put(ybuf+ycount+SENT_OFFSET, 1, MPI_DOUBLE, mpi.yhi, ycount+SENT_OFFSET, 1, MPI_DOUBLE, ywin);
              MPI_Win_flush(mpi.yhi, ywin);
            }
          }
          if (k == 0) {
            /* Do comms if internal boundary */
            if (mpi.zlo != MPI_PROC_NULL) {
              printf("%d: send k=0\n", mpi.rank);
              /* Poll for safe to send signal - lock required around access */
              int safe = 0;
              while (!safe) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, zwin);
                if (zbuf[zcount+SAFE_OFFSET] == SAFE_SIGNAL)
                  safe = 1;
                MPI_Win_unlock(mpi.rank, zwin);
              }

              /* Reset signal */
              zbuf[zcount+SAFE_OFFSET] = NULL_SIGNAL;

              /* Put payload */
              MPI_Put(zbuf, zcount, MPI_DOUBLE, mpi.zlo, 0, zcount, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zlo, zwin);

              /* Send sent signal */
              zbuf[zcount+SENT_OFFSET] = SENT_SIGNAL;
              MPI_Put(zbuf+zcount+SENT_OFFSET, 1, MPI_DOUBLE, mpi.zlo, zcount+SENT_OFFSET, 1, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zlo, zwin);
            }
          }
          else {
            /* Do comms if internal boundary */
            if (mpi.zhi != MPI_PROC_NULL) {
              printf("%d: send k/=0\n", mpi.rank);
              /* Poll for safe to send signal - lock required around access */
              int safe = 0;
              while (!safe) {
                MPI_Win_lock(MPI_LOCK_SHARED, mpi.rank, 0, zwin);
                if (zbuf[zcount+SAFE_OFFSET] == SAFE_SIGNAL)
                  safe = 1;
                MPI_Win_unlock(mpi.rank, zwin);
              }

              /* Reset signal */
              zbuf[zcount+SAFE_OFFSET] = NULL_SIGNAL;

              /* Put payload */
              MPI_Put(zbuf, zcount, MPI_DOUBLE, mpi.zhi, 0, zcount, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zhi, zwin);

              /* Send sent signal */
              zbuf[zcount+SENT_OFFSET] = SENT_SIGNAL;
              MPI_Put(zbuf+zcount+SENT_OFFSET, 1, MPI_DOUBLE, mpi.zhi, zcount+SENT_OFFSET, 1, MPI_DOUBLE, zwin);
              MPI_Win_flush(mpi.zhi, zwin);
            }
          }
          time.comms += MPI_Wtime() - comtime;

        } /* End nchunks loop */

        MPI_Barrier(MPI_COMM_WORLD);
        printf("%d: done oct %d\n", mpi.rank, i+2*j+4*k);
      } /* End i loop */
    } /* End j loop */
  } /* End k loop */

  printf("%d: done!\n", mpi.rank);

  /* End the timer */
  double tock = MPI_Wtime();

  time.sweeping = tock-tick;

  end_one_sided_sweep(&ywin, &zwin, mpi.ylo, mpi.yhi, mpi.zlo, mpi.zhi);

  time.setup += MPI_Wtime() - tock;

  return time;
}

/* Init MPI buffers and set up one-sided comms */
void init_one_sided_sweep(const int ycount, const int zcount, double **ybuf, double **zbuf, MPI_Win *ywin, MPI_Win *zwin, const int ylo, const int yhi, const int zlo, const int zhi) {
  /* Allocate MPI window buffer*/
  MPI_Info info;
  MPI_Info_create(&info);
  MPI_Info_set(info, "same_disp_unit", "true");
  /* Size of payload plus 2 values to allow for safe and done signals */
  MPI_Win_allocate(sizeof(double)*(ycount+2), sizeof(double), info, MPI_COMM_WORLD, ybuf, ywin);
  MPI_Win_allocate(sizeof(double)*(zcount+2), sizeof(double), info, MPI_COMM_WORLD, zbuf, zwin);

  /* Start passive communication epoch - expose this rank's windows to its 4 neighbours */
  MPI_Win_lock(MPI_LOCK_SHARED, ylo, 0, *ywin);
  MPI_Win_lock(MPI_LOCK_SHARED, yhi, 0, *ywin);
  MPI_Win_lock(MPI_LOCK_SHARED, zlo, 0, *zwin);
  MPI_Win_lock(MPI_LOCK_SHARED, zhi, 0, *zwin);

  /* Reset the signal */
  (*ybuf)[ycount+SAFE_OFFSET] = NULL_SIGNAL;
  (*ybuf)[ycount+SENT_OFFSET] = NULL_SIGNAL;
  (*zbuf)[ycount+SAFE_OFFSET] = NULL_SIGNAL;
  (*zbuf)[ycount+SENT_OFFSET] = NULL_SIGNAL;
}

/* Finish up any MPI things */
void end_one_sided_sweep(MPI_Win *ywin, MPI_Win *zwin, const int ylo, const int yhi, const int zlo, const int zhi) {
  /* End passive communication epoch */
  MPI_Win_unlock(ylo, *ywin);
  MPI_Win_unlock(yhi, *ywin);
  MPI_Win_unlock(zlo, *zwin);
  MPI_Win_unlock(zhi, *zwin);

  /* Free MPI message buffers */
  MPI_Win_free(ywin);
  MPI_Win_free(zwin);
}

