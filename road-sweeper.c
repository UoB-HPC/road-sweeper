
#include "comms.h"
#include <mpi.h>
#include "options.h"
#include "sweep.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERSION "0.0"

void parse_args(mpistate mpi, int argc, char *argv[], options *opt);

int main(int argc, char *argv[]) {

  /* Structure to hold MPI status */
  mpistate mpi;

  /* Initilise MPI, requesting the highest thread support possible */
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi.thread_support);

  /* Get MPI rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.nprocs);

  /* Structure to hold runtime options - set defaults */
  options opt = {
    .nchunks = 1,
    .chunklen = 1,
    .ny = 1,
    .nz = 1,
    .nang = 10,
    .ng = 16
  };

  parse_args(mpi, argc, argv, &opt);

  /* Print MPI thread support */
  if (mpi.rank == 0) {
    printf("Road Sweeper\n");
    printf("Version: %s\n", VERSION);

    printf("MPI thread support: ");
    switch (mpi.thread_support) {
      case MPI_THREAD_SINGLE:
        printf("MPI_THREAD_SINGLE\n");
        break;
      case MPI_THREAD_FUNNELED:
        printf("MPI_THREAD_FUNNELED\n");
        break;
      case MPI_THREAD_SERIALIZED:
        printf("MPI_THREAD_SERIALIZED\n");
        break;
      case MPI_THREAD_MULTIPLE:
        printf("MPI_THREAD_MULTIPLE\n");
        break;
    }

  }

  /* Perform 2D decomposition in YZ */
  decompose(&mpi);

  /* Print runtime options */
  if (mpi.rank == 0) {
    printf("MPI processes: %d\n", mpi.nprocs);
    printf("Effective mesh: %d x %d x %d\n", opt.nchunks*opt.chunklen, mpi.npey*opt.ny, mpi.npez*opt.nz);
    printf("Decomposition: %d x %d\n", mpi.npey, mpi.npez);
    printf("Subdomain: %d x %d x %d\n", opt.nchunks*opt.chunklen, opt.ny, opt.nz);
    printf("Chunks per octant: %d\n", opt.nchunks);
    printf("Cells per chunk: %d\n", opt.chunklen);
    printf("Number of angles: %d\n", opt.nang);
    printf("Number of energy groups: %d\n", opt.ng);
    printf("\n");
  }

  double tick = MPI_Wtime();

  double time = serial_sweep(mpi, opt);

  double tock = MPI_Wtime();

  if (mpi.rank == 0) {
    printf("Serial sweeps\n");
    printf("  Time:       %11.6lf s\n", time);
    printf("  + overhead: %11.6lf s\n", tock-tick-time);
  }

  MPI_Finalize();

}

void parse_args(mpistate mpi, int argc, char *argv[], options *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--nchunks") == 0) {
      opt->nchunks = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--chunklen") == 0) {
      opt->chunklen = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--ny") == 0) {
      opt->ny = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--nz") == 0) {
      opt->nz = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--nang") == 0) {
      opt->nang = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--ng") == 0) {
      opt->ng = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--help") == 0) {
      if (mpi.rank == 0) {
        printf("Usage: %s [OPTIONS]\n", argv[0]);
        printf("\t--nchunks  N\tSet number of chunks per octant\n");
        printf("\t--chunklen N\tNumber of cells in x-dimension per chunk\n");
        printf("\t--ny       N\tNumber of cells per subdomain in y-dimension\n");
        printf("\t--nz       N\tNumber of cells per subdomain in z-dimension\n");
        printf("\t--nang     N\tNumber of angles per cell\n");
        printf("\t--ng       N\tNumber of energy groups\n");
      }
      /* Exit nicely */
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }
    else {
      if (mpi.rank == 0) {
        printf("Unknown option: %s\n", argv[i]);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
      }
    }
  }
}

