
#include "comms.h"
#include <mpi.h>
#include <stdio.h>


int main(int argc, char *argv[]) {

  /* Structure to hold MPI status */
  mpistate mpi;

  /* Initilise MPI, requesting the highest thread support possible */
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi.thread_support);

  /* Get MPI rank */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi.nprocs);

  /* Print MPI thread support */
  if (mpi.rank == 0) {
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

  MPI_Finalize();

}

