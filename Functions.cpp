// Define all functions here

#include <windows.h>
#include <stdio.h>
#include "Functions.h"
#include <iostream> 
#include <math.h>

using namespace std;

inline Flt *TraceRay(int pixel_x, int pixel_y, int w, int h, Point E, Vec V, Flt f, Vec U,Flt phi)
{
	Flt theta;
	Point M;
	Vec A;
	Vec B;
	Vec X;
	Vec Y;
	Flt sx;
	Flt sy;
	
	static Flt P[3];

	Flt L = sqrt(*U * (*U) + *(U+1) * (*(U+1)) + *(U+2) * (*(U+2)));
	
	*U = *U / L;
	*(U+1) = *(U+1)/ L;
	*(U+2) = *(U+2) / L;

	*A = (*(V+1)) * (*(U+2)) - (*(V+2)) * (*(U+1));
	*(A+1) = (*(V + 2)) * (*U) - (*V) * (*(U + 2));
	*(A+2) = (*V) * (*(U + 1)) - (*(V + 1)) * (*U);

	*B = (*(A + 1)) * (*(V + 2)) - (*(A + 2)) * (*(V + 1));
	*(B + 1) = (*(A + 2)) * (*V) - (*A) * (*(V + 2));
	*(B + 2) = (*A) * (*(V + 1)) - (*(A + 1)) * (*V);

	theta = phi;
	
	for (int i = 0; i < 3; i++)
	{
		Flt two = (Flt)2;
		M[i] = *(E + i) + f*( *(V + i) );
		Y[i] = f*tan(phi / two)* (*(B+i));
		X[i] = f*tan(theta /two)* (*(A+i));
	}
	
	Flt w_inv =(Flt)1 / (Flt)w;
	Flt h_inv = (Flt)1 / (Flt)h;
	
	sx = ((Flt)pixel_x + 0.5f)* w_inv;
	sy = ((Flt)pixel_y + 0.5f) *h_inv;

	for (int l = 0; l < 4; l++)
	{
		P[l] = M[l] + (2.0f * sx - 1.0f)*X[l] + (2.0f * sy - 1.0f)*Y[l];
	}

	return P;
}

inline bool IntersectSphere(Point P, Point C, Flt R, Point E, Flt t)
{
	Flt dx, dy, dz, a, b, c;
	float t0, t1; // solutions for t if the ray intersects 
	dx = *P - *E;
	dy = *(P+1) - *(E+1);
	dz = *(P+2) - *(E+2);

	a = dx*dx + dy*dy + dz*dz;
	b = 2.0f * dx*( *E - *C) + 2.0f * dy*(*(E+1) - *(C+1)) + 2.0f * dz*(*(E+2) - *(C+2));
	c = *C * *C + *(C+1) * *(C+1) + *(C+2) * *(C+2) + *E * *E + *(E+1) * *(E+1) + *(E+2) * *(E+2) - 2.0f * (*C * *E + *(C+1) * *(E+1) + *(C+2) * *(E+2)) - R*R;

	if (!solveQuadratic(a, b, c, t0, t1))
		return false;
	else
	{
		if (t0 > t1) 
		{
			Flt temp;
			temp = t0; t0 = t1; t1 = temp;
		}
		if (t0 < 0) {
			t0 = t1; // if t0 is negative, let's use t1 instead 
			if (t0 < 0) return false; // both t0 and t1 are negative 
		}

		t = t0;
		return true;
	}
}
/* intersection point
	x = E[0] + t*(P[0]-E[0]);
	y = E[1] + t*dy;
	z = E[2] + t*dz;
*/
inline bool solveQuadratic(Flt a, Flt b, Flt c, Flt x0,  Flt x1)
{
	Flt D = b * b - 4.0f * a * c;
	if (D < 0.0f) return false;
	else if (D == 0.0f) x0 = x1 = -0.5f * b / a;
	else {
		Flt q = (b > 0.0f) ?
			-0.5f * (b + sqrt(D)) :
			-0.5f * (b - sqrt(D));
		x0 = q / a;
		x1 = c / q;
	}
	if (x0 > x1) {
		std::swap(x0, x1);
	}
	return true;
}

/*
inline bool IntersectTriangle(Flt polyvertices , Point P)
{	
		P = point on image plane
		// v0, v1, v2 are vertices in counter clockwise direction

		n = crossproduct((v1 - v0), (v2 - v0));

	a1 = crossproduct((P1 - P0), (x - v0));
	a = a1[0] * n[0] + a1[1] * n[1] + a1[2] * n[2];

	b1 = crossproduct((P2 - P1), (x - P1));
	b = b1[0] * n[0] + b1[1] * n[1] + b1[2] * n[2];

	c1 = crossproduct((P0 - P2), (x - P2));
	c = c1[0] * n[0] + c1[1] * n[1] + c1[2] * n[2];

	if (a >= 0 && b >= 0 && c >= 0)
	{
	D = dotproduct n, (p1-_p0)
 R = dire of ray 
	}
	

}
*/