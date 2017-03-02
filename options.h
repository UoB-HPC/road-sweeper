
/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Number of chunks per octant */
  int nchunks;

  /* Phony compute time in seconds */
  double work;

  /* Number of doubles to send per MPI communication */
  int msglen;

}  options;

