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


/*
 * Global structures for MPI communication
 */

#pragma once

#include "options.h"

typedef struct mpistate {

  /* Level of thread support provided by MPI implementation */
  int thread_support;

  /* Process rank */
  int rank;

  /* Number of ranks in MPI_COMM_WORLD */
  int nprocs;

  /* Number of ranks in 2D decomposition */
  int npey;
  int npez;

  /* Process rank in 2D decomposition */
  int y;
  int z;

  /* Neighbour ranks - 2D decomposition of 3D domain */
  int yhi;
  int ylo;
  int zhi;
  int zlo;

} mpistate;

void decompose(mpistate *mpi);
void decompose_mesh(mpistate *mpi, options *opt);

