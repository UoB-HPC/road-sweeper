
default: road-sweeper

MPICC = mpicc
CFLAGS = -O3

road-sweeper: road-sweeper.c comms.h options.h
	$(MPICC) $(CFLAGS) road-sweeper.c -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

