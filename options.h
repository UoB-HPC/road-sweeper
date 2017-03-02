
/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Number of chunks per octant */
  int nchunks;

  /* Number of cells per chunk */
  int chunklen;

  /* Pencil size - determines Y and Z dimension */
  int ny;
  int nz;

  /* Angles per cell */
  int nang;

  /* Energy groups */
  int ng;

  /* Sweeper version */
  int version;

}  options;

