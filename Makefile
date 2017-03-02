
default: road-sweeper

MPICC = mpicc
CFLAGS = -O3 -std=c99

SRC = road-sweeper.c comms.c serialsweep.c compute.c
HEADER = options.h comms.h serialsweep.h compute.h

road-sweeper: $(SRC) $(HEADER)
	$(MPICC) $(CFLAGS) $(SRC) -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

