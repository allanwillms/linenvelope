# linenvelope

> [!IMPORTANT]
> This software is obsolete: LinearBand.jl is an algorithm for the same purpose which is an improvement over LINENVELOPE.

LINENVELOPE is software for bounding time series data with a piecewise linear band.

Copyright 2008 Allan Willms.

This C software source code is distributed under the GNU General Public License, Version 3.

Bug reports and comments should be sent to Allan Willms.
## Download
There are just two files: linenvelope.c and linenvelope.h.

## Description
LINENVELOPE constructs a piecewise linear band (two piecewise linear curves differing by a constant vertical offset) which bounds a given data set {(t i, y i)}. The quality of the fit is governed by two parameters, one which defines when two adjacent step discontinuities in the slope in opposite directions are considered "close", and the other which defines the maximum number of such close points.

The algorithm is described in

- A.R. Willms, Bounding data with a piecewise linear band, SIAM J. Sci. Comput. 31 (3) (2009) 2361-2367. 

