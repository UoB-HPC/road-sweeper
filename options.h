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


/*
 * Structure to hold runtime options
 */

#pragma once

typedef struct options {

  /* Number of chunks per octant */
  int nchunks;

  /* Number of cells per chunk */
  int chunklen;

  /* Pencil size - determines Y and Z dimension */
  int ny;
  int nz;

  /* Global mesh size for strong scaling */
  int gny;
  int gnz;

  /* Angles per cell */
  int nang;

  /* Energy groups */
  int ng;

  /* Sweeper version */
  int version;

  /* Strong scaling run? */
  int strong;

}  options;

