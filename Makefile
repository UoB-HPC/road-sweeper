
default: road-sweeper

MPICC = mpicc
CFLAGS = -O3

road-sweeper: road-sweeper.c
	$(MPICC) $(CFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

