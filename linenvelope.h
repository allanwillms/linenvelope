/* Copyright 2008 Allan R. Willms
 *
 * This file is part of linenvelope.
 *
 * Linenvelope is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef TypedefRange
#define TypedefRange
typedef double Range[2];
#endif

int linenvelope(int ndata, double *t, double *y, double tmin, int maxclose,
	double **tout, Range **yout);
