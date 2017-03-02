
/*
 * Sweep options
 * Header provides method declarations for various sweeper options.
 * Sweeps should return a time in seconds spent sweeping, excluding
 * setup costs such as buffer allocation.
 */

#pragma once

#include "comms.h"
#include "options.h"

/*
 * Vanilla serial KBA sweeper
 * For each octant, groups are computed serially in turn.
 * We send a message of all the angles on the face for one energy group.
 */
double serial_sweep(mpistate mpi, options opt);

