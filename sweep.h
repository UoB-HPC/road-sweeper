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
 * Sweep options
 * Header provides method declarations for various sweeper options.
 * Sweeps should return a struct consisting of times in seconds spent
 * sweeping, setup costs such as buffer allocation, etc.
 */

#pragma once

#include "comms.h"
#include "options.h"

typedef struct timings {
  /* Total time spent sweeping, excluding setup */
  double sweeping;

  /* Setup and tear down costs */
  double setup;

  /* Time in comms */
  double comms;

} timings;

/*
 * Vanilla serial KBA sweeper
 * For each octant, groups are computed serially in turn.
 * We send a message of all the angles on the face for one energy group.
 */
timings serial_sweep(mpistate mpi, options opt);

/*
 * Parallel over groups, sending all groups in comms
 */
timings par_group_sweep(mpistate mpi, options opt);

/*
 * Parallel over groups, using OpenMP threads to run
 * individual group sweeps concurrently.
 * This results in calling MPI from within the threaded
 * region
 */
timings par_mpi_sweep(mpistate mpi, options opt);

/* Same as above, but uses multiple locks - one per thread */
timings par_mpi_multi_lock_sweep(mpistate mpi, options opt);

/* One sided sweeper, with parallel groups */
timings one_sided_sweep(mpistate mpi, options opt);

