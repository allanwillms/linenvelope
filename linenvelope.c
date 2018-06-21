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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linenvelope.h"

#define MIN(a,b) ((a) <= (b) ? (a) : (b))
#define MAX(a,b) ((a) >= (b) ? (a) : (b))

/* linenvelope  version 1
 * This function computes a piece-wise linear envelope around given data (t,y).  The
 * envelope has the same height all the way along.  The points where the slope of
 * the envelope changes are called control points and at a control point either the
 * top of envelope (upper control point) or the bottom of the envelope (lower
 * control point) is exactly on the data point.  Adjacent control points are
 * considered "close" if their t values are within tmin of each other and they are
 * of opposite type (one lower, one upper).  The number of close points allowed is
 * maxclose.  The algorithm halts when the envelope height is zero or the next 
 * control point to be added would make too many close points.  The user can choose
 * to have the band output at every input t value or at just the control points.  In
 * the latter case, "false" control points are added at the beginning and end of the 
 * data before output.
 *
 *   Calling sequence:
 *   number_of_control_points = linenvelope(int number_of_data_points,
 *		double *time, double *y, double tmin, int maxclose,
 *		double **tout,  Range **yout)
 *  where tout is the address of a memory location that is a pointer to a double 
 *  and yout is the address of a memory location that is a pointer to a Range data 
 *  type (double[2]).
 *  If tout is NULL then the function allocates space for the envelope at every time
 *  point and puts a pointer to this space in *yout, and fills this space.
 *  If tout is not NULL then the function allocates space only for the control
 *  points (t values in *tout and envelope values in *yout).
 */
enum pointtype {lower=0, upper=1};

struct s_control_point {
    int index;
    enum pointtype type;
    double reduction;
    int valid;
    int newcpt;
    int newtype;
};

int initialize(struct s_control_point **cpt, int ndata, double *y, double *width);

void compute_reduction(int ncpt, struct s_control_point **cpt, double width,
	int ndata, double *t, double *y);

int reduce_cpts(int ncpt, double *width, struct s_control_point **cpt);

int insert_cpt(int ind, int *ncpt, struct s_control_point **cpt, double width,
	double tmin, double *t, int maxclose, int *nclose);

void write_band(int ndata, double *t, double *y, int *ncpt, 
	struct s_control_point **cpt, double height, double **tout, Range **yout);

int linenvelope(int ndata, double *t, double *y, double tmin, int maxclose,
	double **tout, Range **yout) {

    double height;
    int i,halt,ncpt,nclose;
    struct s_control_point **cpt;

    nclose = 0;
    cpt = (struct s_control_point **) malloc(ndata*sizeof(struct s_control_point *));
    ncpt = initialize(cpt,ndata,y,&height);
    if ((ncpt==1) || (t[cpt[1]->index] - t[cpt[0]->index] < tmin && 
		++nclose > maxclose)) halt = 1;
    else halt = 0;
    while (!halt) {
	compute_reduction(ncpt,cpt,height,ndata,t,y);
	i = reduce_cpts(ncpt,&height,cpt);
	if (height > 0)
	    halt = insert_cpt(i,&ncpt,cpt,height,tmin,t,maxclose,&nclose);
	else
	    halt = 1;
    }
    write_band(ndata,t,y,&ncpt,cpt,height,tout,yout);
    return ncpt;
}

void write_band(int ndata, double *t, double *y, int *ncpt, 
	struct s_control_point **cpt, double height, double **tout, Range **yout) {

    double z[2],slope;
    int i,j,ntruecpt;

    ntruecpt = *ncpt;
    if (tout == NULL) {
	/* In this case, the user wants ranges for all the y data. */
	*yout = (Range *) malloc(ndata*sizeof(Range));
	/* First put in the points up to the first control point. */
	if (cpt[0]->type == lower)
	    z[0] = y[cpt[0]->index];
	else
	    z[0] = y[cpt[0]->index] - height;
	z[1] = z[0] + height;
	for (i=0; i<=cpt[0]->index; i++) {
	    (*yout)[i][0] = z[0];
	    (*yout)[i][1] = z[1];
	}
	/* Now put in the points between the first and last control point. */
	j = 0;
	z[1] = z[0];
	for (i=cpt[0]->index+1; i<=cpt[*ncpt-1]->index; i++) {
	    if (i > cpt[j]->index) {
		j++;
		z[0] = z[1];
		if (cpt[j]->type == lower)
		    z[1] = y[cpt[j]->index];
		else
		    z[1] = y[cpt[j]->index] - height;
		slope = (z[1]-z[0])/(t[cpt[j]->index] - t[cpt[j-1]->index]);
	    }
	    (*yout)[i][0] = z[1] + (t[i]-t[cpt[j]->index])*slope;
	    (*yout)[i][1] = (*yout)[i][0] + height;
	}
	/* Now j = ncpt-1 (the last control point).  Put in the points from the
	 * last control point to the end. */
	if (cpt[j]->type == lower)
	    z[0] = y[cpt[j]->index];
	else
	    z[0] = y[cpt[j]->index] - height;
	z[1] = z[0] + height;
	for (i=cpt[j]->index+1; i<ndata; i++) {
	    (*yout)[i][0] = z[0];
	    (*yout)[i][1] = z[1];
	}
	/* Add in the count of fake control points. */
	if (cpt[0]->index > 0) (*ncpt)++;
	if (cpt[j-1]->index < ndata-1) (*ncpt)++;
    }
    else {
	/* In this case, the user wants just the control point times and y ranges.
	 * But we also put in fake control points at the start and finish if
	 * necessary. */
	if (cpt[0]->index > 0) (*ncpt)++;
	if (cpt[ntruecpt-1]->index < ndata-1) (*ncpt)++;
	*tout = (double *) malloc(*ncpt*sizeof(double));
	*yout = (Range *) malloc(*ncpt*sizeof(Range));
	if (cpt[0]->index > 0) {
	    (*tout)[0] = t[0];
	    if (cpt[0]->type == lower) {
		(*yout)[0][0] = y[cpt[0]->index];
		(*yout)[0][1] = y[cpt[0]->index] + height;
	    }
	    else {
		(*yout)[0][0] = y[cpt[0]->index] - height;
		(*yout)[0][1] = y[cpt[0]->index];
	    }
	    j = 1;
	}
	else j = 0;
	for (i=0; i<ntruecpt; i++) {
	    (*tout)[i+j] = t[cpt[i]->index];
	    if (cpt[i]->type == lower) {
		(*yout)[i+j][0] = y[cpt[i]->index];
		(*yout)[i+j][1] = y[cpt[i]->index] + height;
	    }
	    else {
		(*yout)[i+j][0] = y[cpt[i]->index] - height;
		(*yout)[i+j][1] = y[cpt[i]->index];
	    }
	}
	if (cpt[ntruecpt-1]->index < ndata-1) {
	    i = *ncpt - 1;
	    (*tout)[i] = t[ndata-1];
	    if (cpt[ntruecpt-1]->type == lower) {
		(*yout)[i][0] = y[cpt[ntruecpt-1]->index];
		(*yout)[i][1] = y[cpt[ntruecpt-1]->index] + height;
	    }
	    else {
		(*yout)[i][0] = y[cpt[ntruecpt-1]->index] - height;
		(*yout)[i][1] = y[cpt[ntruecpt-1]->index];
	    }
	}
    }
    for (i=0; i<ntruecpt; i++) free(cpt[i]);
    free(cpt);
    return;
}

int initialize(struct s_control_point **cpt, int ndata, double *y, double *height) {
    int i;
    double ymax,ymin;
    struct s_control_point *ptr;

    /* This function puts two control points into the array.  A point at
     * the minimum y value and a point at the maximum y value. They are ordered so
     * that the lowest t value (lowest index) is in cpt[0]. */
    for (i=0; i<2; i++) {
	cpt[i] = (struct s_control_point *) malloc(sizeof(struct s_control_point));
	cpt[i]->index = 0;
	cpt[i]->type = i;
	cpt[i]->valid = 0;
    }
    ymin = y[0];
    ymax = y[0];
    for (i=1; i<ndata; i++) {
	if (y[i] < ymin) {
	    ymin = y[i];
	    cpt[0]->index = i;
	}
	else if (y[i] > ymax) {
	    ymax = y[i];
	    cpt[1]->index = i;
	}
    }
    *height = ymax - ymin;
    if (cpt[0]->index > cpt[1]->index) {
	ptr = cpt[0];
	cpt[0] = cpt[1];
	cpt[1] = ptr;
	i = 2;
    }
    else if (cpt[0]->index == cpt[1]->index) {
	i = 1;
	free(cpt[1]);
    }
    else i = 2;
    return i;
}

void compute_reduction(int ncpt, struct s_control_point **cpt, double height,
	int ndata, double *t, double *y) {
    int i,j,k;
    int last;
    enum MODE {slide,pivot} mode;
    double temp,slope;

    static int count=0;

    for (i=0; i<ncpt; i++) {
	if (!cpt[i]->valid) {
	    cpt[i]->reduction = height;
	    cpt[i]->valid = 1;
	    if (i==0) { /* Compute reductions to left of control point 0. */
		/* Loop through the left neighbouring points and determine
		 * the smallest reduction allowed. */
		for (j=0; j<cpt[0]->index; j++) {
		    if (cpt[0]->type == lower) {
			temp = height - (y[j] - y[cpt[0]->index]);
		    }
		    else {
			temp = height - (y[cpt[0]->index] - y[j]);
		    }
		    if (temp < cpt[0]->reduction) {
			cpt[0]->reduction = temp;
			cpt[0]->newcpt = j;
			cpt[0]->newtype = 1 - cpt[0]->type;
		    }
		}
		count += cpt[0]->index;
	    }
	    /* Compute reductions to the right of this control point.
	     * If the next control point is of the same type, then the mode is 
	     * "slide" and only one boundary of the band moves as the free ends are
	     * moved.  If the next control point is of opposite type, then the mode
	     * is "pivot" and both boundaries of the band move as the free points
	     * move. */
	    if (i==ncpt-1) {
		last = ndata;
		mode = slide;
		slope = 0.0;
	    }
	    else if (cpt[i+1]->type == cpt[i]->type) {
		last = cpt[i+1]->index;
		mode = slide;
		slope = (y[cpt[i+1]->index] - y[cpt[i]->index])/
		    (t[cpt[i+1]->index] - t[cpt[i]->index]);
	    }
	    else {
		last = cpt[i+1]->index;
		mode = pivot;
	    }
	    /* Loop through the neighbouring points and determine
	     * the smallest reduction allowed. */
	    if (mode == slide) {
		for (j=cpt[i]->index+1; j<last; j++) {
		    temp = y[j] + slope*(t[cpt[i]->index] - t[j]);
		    if (cpt[i]->type == lower) 
			temp = height - (temp - y[cpt[i]->index]);
		    else 
			temp = height - (y[cpt[i]->index] - temp);
		    if (temp < cpt[i]->reduction) {
			cpt[i]->reduction = temp;
			cpt[i]->newcpt = j;
			cpt[i]->newtype = 1 - cpt[i]->type;
		    }
		}
		count += last-cpt[i]->index-1;
	    }
	    else {  /* mode==pivot */
		/* First move the free end for this control point. */
		for (j=cpt[i]->index+1; j<last; j++) {
		    slope = (y[cpt[i+1]->index] - y[j])/
			(t[cpt[i+1]->index] - t[j]);
		    temp = y[cpt[i+1]->index] + slope*
			(t[cpt[i]->index] - t[cpt[i+1]->index]);
		    if (cpt[i]->type == lower) 
			temp = height - (temp - y[cpt[i]->index]);
		    else 
			temp = height - (y[cpt[i]->index] - temp);
		    if (temp < cpt[i]->reduction) {
			cpt[i]->reduction = temp;
			cpt[i]->newcpt = j;
			cpt[i]->newtype = 1 - cpt[i]->type;
		    }
		}
		/* Now move the free end for the next control point. */
		for (j=cpt[i]->index+1; j<last; j++) {
		    slope = (y[cpt[i]->index] - y[j])/
			(t[cpt[i]->index] - t[j]);
		    temp = y[cpt[i]->index] + slope*
			(t[cpt[i+1]->index] - t[cpt[i]->index]);
		    if (cpt[i+1]->type == lower) 
			temp = height - (temp - y[cpt[i+1]->index]);
		    else 
			temp = height - (y[cpt[i+1]->index] - temp);
		    if (temp < cpt[i]->reduction) {
			cpt[i]->reduction = temp;
			cpt[i]->newcpt = j;
			cpt[i]->newtype = cpt[i]->type;
		    }
		}
		count += 2*(last-cpt[i]->index-1);
	    }
	}
    }
    printf("count = %d\n",count);
    return;
}

int reduce_cpts(int ncpt, double *height, struct s_control_point **cpt) {
    int i,ind;
    double smallest;

    /* find the smallest reduction and apply it to all points. */
    ind = 0;
    smallest = cpt[0]->reduction;
    for (i=1; i<ncpt; i++) {
	if (cpt[i]->reduction < smallest) {
	    ind = i;
	    smallest = cpt[i]->reduction;
	}
    }
    *height -= smallest;
    for (i=0; i<ncpt; i++) cpt[i]->reduction -= smallest;
    return ind;
}

int insert_cpt(int ind, int *ncpt, struct s_control_point **cpt, double height,
	double tmin, double *t, int maxclose, int *nclose) {
    int i,newind,newcpt,newtype;

    newcpt = cpt[ind]->newcpt;
    newtype = cpt[ind]->newtype;
    if (newcpt > cpt[ind]->index) newind = ind+1;
    else newind = ind;
    /* If the new point is too close in time to ind, or if there is a neighbouring
     * control point that is too close in time and has the same type as
     * that of ind, then this is a close point.  Only insert if there are not too
     * many close points already.  Otherwise halt the algorithm. */
    i = 0;
    if (newtype != cpt[ind]->type) {
	if ((newind == 0 && t[cpt[ind]->index] - t[newcpt] < tmin) ||
	    (newind > 0 && t[newcpt] - t[cpt[ind]->index] < tmin)) i++;
    }
    if (newind > 0 && newind < *ncpt) {
	if (newtype != cpt[ind+1]->type && t[cpt[ind+1]->index] - t[newcpt] < tmin) 
	    i++;
	if (cpt[ind]->type != cpt[ind+1]->type && 
		t[cpt[ind+1]->index] - t[cpt[ind]->index] < tmin) 
	    i--;
    }
    if (((*nclose) += i) > maxclose) return 1;

    /* Mark this point as invalid so that its reduction is recomputed. */
    cpt[ind]->valid = 0;
    /* Shift the control points over to make room for the inserted one. */
    for (i=*ncpt; i>newind; i--) cpt[i] = cpt[i-1];
    (*ncpt)++;
    cpt[newind] = (struct s_control_point *) malloc(sizeof(struct s_control_point));
    cpt[newind]->index = newcpt;
    cpt[newind]->type = newtype;
    cpt[newind]->valid = 0;
    return 0;
}
