
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
    .nchunks = 1,
    .work = 0.1,
    .chunklen = 1,
    .pencil = 1
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
    printf("Decomposition: %d x %d\n", mpi.npey, mpi.npez);
    printf("Subdomain: %d x %d x %d\n", opt.nchunks*opt.chunklen, opt.pencil, opt.pencil);
    printf("Chunks per octant: %d\n", opt.nchunks);
    printf("Cells per chunk: %d\n", opt.chunklen);
    printf("\n");
  }

  /* Allocate message buffers */
  messages message;
  alloc_messages(&message, opt);

  double tick = MPI_Wtime();

  serial_sweep(mpi, opt, message);

  double tock = MPI_Wtime();

  if (mpi.rank == 0) {
    printf("Serial sweep: %11.6lf s\n", tock-tick);
  }

  /* Free MPI buffers */
  free_messages(&message);

  MPI_Finalize();

}

void parse_args(mpistate mpi, int argc, char *argv[], options *opt) {
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--nchunks") == 0) {
      opt->nchunks = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--work") == 0) {
      opt->work = atof(argv[++i]);
    }
    else if (strcmp(argv[i], "--chunklen") == 0) {
      opt->chunklen = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--pencil") == 0) {
      opt->pencil = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--help") == 0) {
      if (mpi.rank == 0) {
        printf("Usage: %s [OPTIONS]\n", argv[0]);
        printf("\t--nchunks  N\tSet number of chunks per octant\n");
        printf("\t--work     t\tSpin lock for t seconds between receive and send\n");
        printf("\t--chunklen N\tNumber of cells in x-dimension per chunk\n");
        printf("\t--pencil   N\tPencil size in Y and Z. e.g. N=1 is a Xx1x1 pencil\n");
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

