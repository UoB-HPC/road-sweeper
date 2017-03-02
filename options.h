
/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Mesh size */
  int N;

  /* Number of chunks per octant */
  int nchunks;

  /* Phony compute time in seconds */
  double work;

}  options;

