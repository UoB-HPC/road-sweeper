
#include "comms.h"
#include <float.h>
#include <mpi.h>
#include <stdlib.h>

void decompose(mpistate *mpi) {

  /* Try options for decomposition - minimise perimeter to area ratio */
  double best = DBL_MAX;

  for (int npey = 1; npey <= mpi->nprocs; npey++) {

    /* Skip non-divisible options */
    if (mpi->nprocs % npey) continue;

    /* Number of ranks for z-dimension */
    int npez = mpi->nprocs / npey;
    if (mpi->nprocs % npez) continue;

    double perimeter = ((mpi->nprocs/npey) + (mpi->nprocs/npez)) * 2.0;
    double area = (mpi->nprocs/npey) * (mpi->nprocs/npez);
    double ratio = perimeter / area;

    /* Save best so far */
    if (ratio < best) {
      best = ratio;
      mpi->npey = npey;
      mpi->npez = npez;
    }
  }

  /* Set ranks along dimension */
  mpi->y = mpi->rank % mpi->npey;
  mpi->z = mpi->rank / mpi->npey;

  /* Set neighbour ranks */
  mpi->ylo = (mpi->y == 0) ? MPI_PROC_NULL : (mpi->y-1) + mpi->z*mpi->npey;
  mpi->yhi = (mpi->y == mpi->npey-1) ? MPI_PROC_NULL : (mpi->y+1) + mpi->z*mpi->npey;
  mpi->zlo = (mpi->z == 0) ? MPI_PROC_NULL : mpi->y + (mpi->z-1)*mpi->npey;
  mpi->zhi = (mpi->z == mpi->npez-1) ? MPI_PROC_NULL : mpi->y + (mpi->z+1)*mpi->npey;
}

/* Allocate the MPI message buffers */
void alloc_messages(messages *message, const options opt) {
  message->ybuf = malloc(sizeof(double)*opt.pencil*opt.chunklen);
  message->zbuf = malloc(sizeof(double)*opt.pencil*opt.chunklen);
}

void free_messages(messages *message) {
  free(message->ybuf);
  free(message->zbuf);
}

