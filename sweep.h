
/*
 * Sweep options
 * Header provides method declarations for various sweeper options.
 */

#pragma once

#include "comms.h"
#include "options.h"

/*
 * Vanilla serial KBA sweeper
 * For each octant, groups are computed serially in turn.
 * We send a message of all the angles on the face for one energy group.
 */
void serial_sweep(mpistate mpi, options opt, messages msg);

