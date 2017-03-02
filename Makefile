
default: road-sweeper

MPICC = mpicc
CFLAGS = -O3 -std=c99
OMP = -fopenmp

SRC = road-sweeper.c comms.c serialsweep.c compute.c pargroupsweep.c
HEADER = options.h comms.h sweep.h compute.h

road-sweeper: $(SRC) $(HEADER)
	$(MPICC) $(CFLAGS) $(SRC) $(OMP) -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

