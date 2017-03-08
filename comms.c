
#include "comms.h"
#include <float.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/* Divide MPI ranks, independent of mesh size - useful for weak scaling */
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

/* Decompose the mesh itself, sharing out any extra cells - useful for strong scaling */
void decompose_mesh(mpistate *mpi, options *opt) {

  /* Try options for decomposition - minimise perimeter to area ratio */
  double best = DBL_MAX;

  for (int npey = 1; npey <= mpi->nprocs; npey++) {

    /* Skip non-divisible options */
    if (mpi->nprocs % npey) continue;

    /* Number of ranks for z-dimension */
    int npez = mpi->nprocs / npey;
    if (mpi->nprocs % npez) continue;

    double perimeter = ((opt->gny/npey) + (opt->gnz/npez)) * 2.0;
    double area = (opt->gny/npey) * (opt->gnz/npez);
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

  /* Set cells per dimension, rouding down */
  opt->ny = opt->gny / mpi->npey;
  opt->nz = opt->gnz / mpi->npez;

  /* Calculate spare cells */
  int extra_y = opt->gny % mpi->npey;
  int extra_z = opt->gnz % mpi->npez;

  /* Distribute cells evenly to ranks */
  if (extra_y && mpi->y < extra_y) {
    opt->ny += 1;
  }
  if (extra_z && mpi->z < extra_z) {
    opt->nz += 1;
  }

  /* Set neighbour ranks */
  mpi->ylo = (mpi->y == 0) ? MPI_PROC_NULL : (mpi->y-1) + mpi->z*mpi->npey;
  mpi->yhi = (mpi->y == mpi->npey-1) ? MPI_PROC_NULL : (mpi->y+1) + mpi->z*mpi->npey;
  mpi->zlo = (mpi->z == 0) ? MPI_PROC_NULL : mpi->y + (mpi->z-1)*mpi->npey;
  mpi->zhi = (mpi->z == mpi->npez-1) ? MPI_PROC_NULL : mpi->y + (mpi->z+1)*mpi->npey;
}

