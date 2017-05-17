Road Sweeper
============

This is a small benchmark code designed to investigate MPI characteristics of a KBA sweep.

> This old brooms had 17 new heads and 14 new handles in its time.

## Building
Just type `make` to produce the `road-sweeper` binary.
The following options can be set:

| Parameter  | Description                     |
|------------|---------------------------------|
| `MPICC`    | Set the MPI C compiler          |
| `CFLAGS`   | Set C compiler flags            |
| `OMP`      | Set OpenMP library flag         |
| `OPTIONS`  | Add additional compiler options |

A preprocessor definition `WORK` can also be set to alter the effective computation delay between messages.
The default value is 50 causing 50 floating point additions to be performed.
The work routine represents the solution of a single point in the entire domain.

## Runtime options
A number of options are passed on the command line.

| Parameter      | Description                                             | Default value   |
|----------------|---------------------------------------------------------|-----------------|
| `--nsweeps N`  | Number of sweeps to perform                             | 1               |
| `--nchunks N`  | Number of spatial chunks per octant                     | 1               |
| `--chunklen N` | Number of cells in x per chunk                          | 1               |
| `--ny N`       | Number of cells in subdomain in y                       | 1               |
| `--nz N`       | Number of cells in subdomain in z                       | 1               |
| `--meshny N`   | Cells in global y-dimension                             | -               |
| `--meshnz N`   | Cells in global z-dimension                             | -               |
| `--strong`     | Perform strong scaling decomposition                    | Off (i.e. weak) |
| `--nang N`     | Number of angles per cell                               | 10              |
| `--ng N`       | Number of groups per cell                               | 16              |
| `--sweep type` | Sweep type (`serial`, `pargroup`, `parmpi`, `mutilock`) | `serial`        |

The number of MPI ranks only need be specified on `mpirun`.
The number of OpenMP threads should be set via the `OMP_NUM_THREADS` environment variable.

### Weak scaling
By default the mesh parameters will be weak scaled (replicated) across the number of MPI ranks.

### Strong scaling
The `--strong`` flag must be specified, along with the `--meshny` and `--meshnz` options.
`--ny` and `--nz` must not be set.
The YZ spatial domain is as evenly as possible across the number of MPI ranks.
Each rank contains the complete X domain, and is of size `nchunks * chunklen` cells.

## Sweep types
A number of sweep types are investigated.
Each is implemented in its own source file.
Each sweeper is allowed to allocate the required MPI message buffer sizes, and any other initilisation it requires.

### Serial
MPI only sweep with no OpenMP threading.
One message is sent per chunk consisting of all the angles for the face cells for *a single* energy group.

### Parallel group
Groups are threaded with OpenMP.
One message is sent per chunk consisting of all the angles for the face cells for *all* energy groups.
The group loop is inside the loop over space and so the OpenMP threads synchronise before the messages are sent.
This is the traditional MPI+OpenMP style implementation.

### Parallel MPI
Sweeps for each groups are threaded with OpenMP.
This means seperate spatial sweeps are running concurrently.
Each *thread* sends one message per chunk consisting of all the angles for the face cells for a single energy group.
Where the MPI library is `MPI_THREAD_SERIALIZED`, a single OpenMP lock is used to syncronise threads.
In the case of `MPI_THREAD_MULTIPLE` no locks are used.
This could be described as an OpenMP+MPI implementation.

### Parallel MPI (multiple locks)
This sweep has a similar regeime to the Parallel MPI sweep, however multiple OpenMP locks are used to synchronise the threads to reduce the lock contention on a single lock.
Each thread has its own lock, which begins locked other than the first thread.
Once MPI communication is complete each thead unlocks the neighbouring lock, resulting in an ordering of the threads.
Note that this behaviour does not conform to the OpenMP specification because locks may only be locked/unlocked by a task region, which corresponds to each thread.

