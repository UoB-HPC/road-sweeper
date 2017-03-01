
#include <stdio.h>

#include <mpi.h>

int main(int argc, char *argv[]) {

  int mpi_thread_support;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &mpi_thread_support);

  int rank, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (rank == 0) {
    printf("MPI thread support: ");
    switch (mpi_thread_support) {
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

