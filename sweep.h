
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

