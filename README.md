# rayt2d_mod - Raytracing utility from Seismic Unix modified to calculate 2-way times for defined reflectors

From seismic unix package

Copyright (c) Colorado School of Mines, 2011.

All rights reserved.                       

## Getting Started

In order to run and compile You should meet these basic criteria:

### Prerequisites

Ubuntu 18.04.1 or later

gcc 7.3.0 or later

### Compiling

You need to download a repository and edit Compile.sh with correct path to source code and linked libraries:

```
 gcc -I./path/to/folder/with/headers -L./path/to/folder/with/libraries \
 /path/to/source/code/rayt2d_mod.c -lpar -lcwp -lm \
 -o /path/for/compiled/programme/rayt2d_mod
```
### Running

Use just name of the program (rayt2d_mod) to bring internal help file

Example of running script for OBS data is in \scripts\Run.sh

## Acknowledgments
Author:  Zhenyue Liu, 10/11/94,  Colorado School of Mines

Trino Salinas, 01/01/96 included the option to handle nonflat
reference surfaces.

Subroutines from Dave Hale's modeling library were adapted in
this code to define topography using cubic splines.
 
## References:
Beydoun, W. B., and Keho, T. H., 1987, The paraxial ray method:
Geophysics, vol. 52, 1639-1653.

Cerveny, V., 1985, The application of ray tracing to the numerical
modeling of seismic wavefields in complex structures, in Dohr, G.,
ED., Seismic shear waves (part A: Theory): Geophysical Press,
Number 15 in Handbook of Geophysical Exploration, 1-124.
 

