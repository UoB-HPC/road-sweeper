
#include "comms.h"
#include <mpi.h>
#include "options.h"
#include "serialsweep.h"
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
    .N = 16,
    .chunk = 1
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

    printf("MPI processes: %d\n", mpi.nprocs);
  }

  /* Perform 2D decomposition */
  decompose(&mpi, opt);
  if (mpi.rank == 0) {
    printf("Decomposition: %d x %d\n", mpi.npey, mpi.npez);
  }

  serial_sweep(mpi, opt);


  MPI_Finalize();

}

void parse_args(mpistate mpi, int argc, char *argv[], options *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--mesh") == 0) {
      opt->N = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--chunk") == 0) {
      opt->chunk = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--help") == 0) {
      if (mpi.rank == 0) {
        printf("Usage: %s [OPTIONS]\n", argv[0]);
        printf("\t--mesh  N\tSet number of cells in each mesh dimension\n");
        printf("\t--chunk C\tSet chunk size for x-dimension\n");
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

