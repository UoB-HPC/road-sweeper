/*
 * This file is part of road-sweeper.
 *
 * road-sweeper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * road-sweeper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with road-sweeper.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "comms.h"
#include <mpi.h>
#include "options.h"
#include "sweep.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VERSION "0.0"

enum sweep {SERIAL, PARGROUP, PARMPI, MULTILOCK};

void print_timings(options opt, timings *times);
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
    .nsweeps = 1,
    .nchunks = 1,
    .chunklen = 1,
    .ny = 1,
    .nz = 1,
    .nang = 10,
    .ng = 16,
    .strong = 0
  };

  parse_args(mpi, argc, argv, &opt);

  /* Print MPI thread support */
  if (mpi.rank == 0) {
    printf("Road Sweeper\n");
    printf("Version: %s\n", VERSION);
    printf("\n");

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
  if (opt.strong) {
    decompose_mesh(&mpi, &opt);
  }
  else {
    decompose(&mpi);
    opt.gny = mpi.npey*opt.ny;
    opt.gnz = mpi.npez*opt.nz;
  }

  /* Print runtime options */
  if (mpi.rank == 0) {
    printf("MPI processes: %d\n", mpi.nprocs);
    printf("Effective mesh: %d x %d x %d\n", opt.nchunks*opt.chunklen, opt.gny, opt.gnz);
    printf("  Cells: %ld\n", (long)opt.nchunks*opt.chunklen*opt.gny*opt.gnz);
    printf("Decomposition: %d x %d\n", mpi.npey, mpi.npez);
    printf("Subdomain: %d x %d x %d\n", opt.nchunks*opt.chunklen, opt.ny, opt.nz);
    printf("Chunks per octant: %d\n", opt.nchunks);
    printf("Cells per chunk: %d\n", opt.chunklen);
    printf("Number of angles: %d\n", opt.nang);
    printf("Number of energy groups: %d\n", opt.ng);
    printf("Numer of sweeps: %d\n", opt.nsweeps);
    printf("====================\n");
    if (opt.version == SERIAL) printf("Running serial sweeper\n");
    else if (opt.version == PARGROUP) printf("Running parallel group sweeper\n");
    else if (opt.version == PARMPI) printf("Running parallel MPI sweeper\n");
    else if (opt.version == MULTILOCK) printf("Running parallel MPI sweeper (multiple locks)\n");
    printf("\n");
  }

  timings *times = malloc(opt.nsweeps*sizeof(timings));

  /* Run the benchmark multiple times */
  for (int s = 0; s < opt.nsweeps; s++) {

    if (opt.version == SERIAL)
      times[s] = serial_sweep(mpi, opt);
    else if (opt.version == PARGROUP)
      times[s] = par_group_sweep(mpi, opt);
    else if (opt.version == PARMPI)
      times[s] = par_mpi_sweep(mpi, opt);
    else if (opt.version == MULTILOCK)
      times[s] = par_mpi_multi_lock_sweep(mpi, opt);
  }

  if (mpi.rank == 0) {
    print_timings(opt, times);
  }

  free(times);

  MPI_Finalize();

}

void print_timings(options opt, timings *times) {
  double total = 0.0;
  int min = 0;
  int max = 0;
  for (int s = 0; s < opt.nsweeps; s++) {
    total += times[s].setup+times[s].sweeping;
    if (times[s].sweeping < times[min].sweeping) {
      min = s;
    }
    if (times[s].sweeping > times[max].sweeping) {
      max = s;
    }
  }
  printf("  Total for all sweeps %11.6lf s\n", total);
  printf("  Time variance:       %11.6lf s\n", times[max].sweeping-times[min].sweeping);
  printf("  Minimum sweep time (sweep #%d)\n", min+1);
  printf("    Total:       %11.6lf s\n", times[min].setup+times[min].sweeping);
  printf("      Setup:     %11.6lf s\n", times[min].setup);
  printf("      Sweeping:  %11.6lf s\n", times[min].sweeping);
  printf("        Comms:   %11.6lf s (%.1lf%%)\n", times[min].comms, times[min].comms/times[min].sweeping*100.0);
  double compute = times[min].sweeping-times[min].comms;
  printf("        Compute: %11.6lf s (%.1lf%%)\n", compute, compute/times[min].sweeping*100.0);
  printf("====================\n");
  printf("\n");
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
    else if (strcmp(argv[i], "--meshny") == 0) {
      opt->gny = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--meshnz") == 0) {
      opt->gnz = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--strong") == 0) {
      opt->strong = 1;
    }
    else if (strcmp(argv[i], "--nang") == 0) {
      opt->nang = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--ng") == 0) {
      opt->ng = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--sweep") == 0) {
      i++;
      if (strcmp(argv[i], "serial") == 0) {
        opt->version = SERIAL;
      }
      else if (strcmp(argv[i], "pargroup") == 0) {
        opt->version = PARGROUP;
      }
      else if (strcmp(argv[i], "parmpi") == 0) {
        opt->version = PARMPI;
      }
      else if (strcmp(argv[i], "multilock") == 0) {
        opt->version = MULTILOCK;
      }
      else {
        if (mpi.rank == 0) {
        printf("Unknown sweep type: %s\n", argv[i]);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
    }
    else if (strcmp(argv[i], "--nsweeps") == 0) {
      opt->nsweeps = atoi(argv[++i]);
    }
    else if (strcmp(argv[i], "--help") == 0) {
      if (mpi.rank == 0) {
        printf("Usage: %s [OPTIONS]\n", argv[0]);
        printf("\t--nsweeps  N\tRun N sweeps\n");
        printf("\t--nchunks  N\tSet number of chunks per octant\n");
        printf("\t--chunklen N\tNumber of cells in x-dimension per chunk\n");
        printf("\t--ny       N\tNumber of cells per subdomain in y-dimension\n");
        printf("\t--nz       N\tNumber of cells per subdomain in z-dimension\n");
        printf("\t--meshny   N\tNumber of cells in y-dimension - not compatible with ny option\n");
        printf("\t--meshnz   N\tNumber of cells in z-dimension - not compatible with nz option\n");
        printf("\t--strong    \tSpecify running strong scaling\n");
        printf("\t--nang     N\tNumber of angles per cell\n");
        printf("\t--ng       N\tNumber of energy groups\n");
        printf("\t--sweep type\tSweeper to run. Options: serial, pargroup, parmpi, multilock\n");
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

  /* Validation */
  if (opt->strong && (opt->gny < 1 || opt->gnz < 1)) {
    if (mpi.rank == 0) {
      printf("Must set --meshny and --meshnz with --strong option\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }
}

