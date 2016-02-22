// declare all functions here

#ifndef _FUNCTIONS_
#define _FUNCTIONS_

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#include <stdio.h>
#include "scene_io.h"

extern int pixel_x;
extern int pixel_y;

 Flt *TraceRay(int i, int j, float w, float h, Point E, Vec V, Flt f, Vec U, Flt phi);
 //void *normalize(Vec A);
 //void crossproduct(Vec A, Vec B, Vec C);
 bool solveQuadratic(Flt a, Flt b, Flt c, Flt x0, Flt x1);
 bool IntersectSphere(Point P, Point C, Flt R, Point E);
 bool IntersectTriangle(Flt polyvertices, Point P);

#endif