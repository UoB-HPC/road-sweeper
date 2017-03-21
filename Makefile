# This file is part of road-sweeper.
#
# road-sweeper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# road-sweeper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with road-sweeper.  If not, see <http://www.gnu.org/licenses/>.


default: road-sweeper

MPICC = mpicc
CFLAGS = -O3 -std=c99
OMP = -fopenmp

SRC = road-sweeper.c comms.c serialsweep.c compute.c pargroupsweep.c parmpisweep.c multilocksweep.c
HEADER = options.h comms.h sweep.h compute.h

road-sweeper: $(SRC) $(HEADER)
	$(MPICC) $(CFLAGS) $(SRC) $(OPTIONS) $(OMP) -o $@

.PHONY: clean
clean:
	rm -f road-sweeper

