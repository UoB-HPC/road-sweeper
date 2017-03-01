
default: road-sweeper

MPICC = mpicc
CFLAGS = -O3

SRC = road-sweeper.c comms.c serialsweep.c
HEADER = options.h comms.h serialsweep.h

road-sweeper: $(SRC) $(HEADER)
	$(MPICC) $(CFLAGS) $(SRC) -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

