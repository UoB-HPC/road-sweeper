
/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Number of chunks per octant */
  int nchunks;

  /* Phony compute time in seconds */
  double work;

  /* Number of cells per chunk */
  int chunklen;

  /* Pencil size - determines Y and Z dimension */
  int ny;
  int nz;

  /* Angles per cell */
  int nang;

  /* Energy groups */
  int ng;

}  options;

