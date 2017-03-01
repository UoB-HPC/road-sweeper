
/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Mesh size */
  int N;

  /* Chunk size */
  int chunk;

  /* Number of messages to send per octant - somewhat equivalent to number of chunks */
  int msgs;

}  options;

