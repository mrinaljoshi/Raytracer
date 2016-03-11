#include <windows.h>
#include <stdio.h>
#include "Functions.h"
#include <iostream> 
#include <math.h>
#include<list>
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

inline bool IntersectSphere_shadow( Vec RayDir, Point C, Flt R, Point P_shadow, Flt *t)
//	Point P, Point C, Flt R, Point E, Flt *t)
{	
	Point P; Point E;
	E[0] = *P_shadow; E[1] = *(P_shadow + 1); E[2] = *(P_shadow + 2);
	P[0] = *P_shadow + *RayDir; P[1] = *(P_shadow+1 )+ *(RayDir+1); P[2] = *(P_shadow +2)+ *(RayDir+2);
		
	Flt dx, dy, dz, a, b, c;
	float t0 = 0.0f; Flt t1 = 0.0f; // solutions for t if the ray intersects 
	dx = *RayDir;
	dy = *(RayDir + 1);
	dz = *(RayDir + 2);

	a = dx*dx + dy*dy + dz*dz;
	b = 2.0f * dx*(E[0] - *C) + 2.0f * dy*(E[1] - *(C + 1)) + 2.0f * dz*(E[2] - *(C + 2));
	c = *C * *C + *(C + 1) * *(C + 1) + *(C + 2) * *(C + 2) + E[0] * E[0] + E[1] * E[1] + E[2] * E[2] - 2.0f * (*C * E[0] + *(C + 1) * E[1] + *(C + 2) * E[2]) - R*R;

	if (!solveQuadratic(a, b, c, &t0, &t1))
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

		*t = t0;
		return true;
	}
}


inline bool IntersectSphere(Point P, Point C, Flt R, Point E, Flt *t)
{
	Flt dx, dy, dz, a, b, c;
	float t0 = 0.0f; Flt t1=0.0f; // solutions for t if the ray intersects 
	dx = *P - *E;
	dy = *(P+1) - *(E+1);
	dz = *(P+2) - *(E+2);

	a = dx*dx + dy*dy + dz*dz;
	b = 2.0f * dx*( *E - *C) + 2.0f * dy*(*(E+1) - *(C+1)) + 2.0f * dz*(*(E+2) - *(C+2));
	c = *C * *C + *(C+1) * *(C+1) + *(C+2) * *(C+2) + *E * *E + *(E+1) * *(E+1) + *(E+2) * *(E+2) - 2.0f * (*C * *E + *(C+1) * *(E+1) + *(C+2) * *(E+2)) - R*R;

	if (!solveQuadratic(a, b, c, &t0, &t1))
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

		*t = t0;
		return true;
	}
}

inline bool solveQuadratic(Flt a, Flt b, Flt c, Flt *x0,  Flt *x1)
{
	Flt D = b * b - 4.0f * a * c;
	if (D < 0.0f) return false;
	else if (D == 0.0f) *x0 = *x1 = -0.5f * b / a;
	else {
		Flt q = (b > 0.0f) ?
			-0.5f * (b + sqrt(D)) :
			-0.5f * (b - sqrt(D));
		*x0 = q / a;
		*x1 = c / q;
	}
	if (*x0 > *x1) {
		std::swap(*x0, *x1);
	}
	return true;
}


inline bool IntersectTriangle(arr polyvertices, Point P, Point E, Flt *t)
{	//E = camera center
		//P = point on image plane
		// v0, v1, v2 are vertices in counter clockwise direction
	Vec v0; Vec v1; Vec v2;
	v0[0] = *polyvertices;
	v0[1] = *(polyvertices + 1);
	v0[2] = *(polyvertices + 2);
	v1[0] = *(polyvertices + 3);
	v1[1] = *(polyvertices + 4);
	v1[2] = *(polyvertices + 5);
	v2[0] = *(polyvertices + 6);
	v2[1] = *(polyvertices + 7);
	v2[2] = *(polyvertices + 8);

	Vec v10; Vec v20; Vec N; Flt d; Vec R; Vec v21; Vec v02;
	Vec pvec, qvec;
	Vec tvec;
	for (int i = 0; i < 3; i++)
	{
		v10[i] = v1[i] - v0[i];
		v20[i] = v2[i] - v0[i];
		v21[i] = v2[i] - v1[i];
		v02[i] = v0[i] - v2[i];
		R[i] = *(P + i) - *(E + i);   // should it be *P,*E
		tvec[i] = *(E + i) - v0[i];
	}

	// pvec = R crossproduct v20

	pvec[0] = R[1] * v20[2] - R[2] * v20[1];
	pvec[1] = R[2] * v20[0] - R[0] * v20[2];
	pvec[2] = R[0] * v20[1] - R[1] * v20[0];

	Flt det;
	det = v10[0] * pvec[0] + v10[1] * pvec[1] + v10[2] * pvec[2];

	if (det < 0.05f) return false;
	if (fabs(det) < 0.05f) return false;

	Flt invDet, u, v;
	invDet = 1.0f / det;

	u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2])* invDet;
	if (u < 0 || u>1) return false;

	qvec[0] = tvec[1] * v10[2] - tvec[2] * v10[1];
	qvec[1] = tvec[2] * v10[0] - tvec[0] * v10[2];
	qvec[2] = tvec[0] * v10[1] - tvec[1] * v10[0];

	v = (R[0] * qvec[0] + R[1] * qvec[1] + R[2] * qvec[2]) *invDet;

	if (v < 0 || u + v >1) return false;

	*t = (v20[0] * qvec[0] + v20[1] * qvec[1] + v20[2] * qvec[2]) * invDet;

	return true;
}

inline bool IntersectTriangle_shadow(arr polyvertices, Point P_shadow, Vec RayDir, Flt *t)
{
	Point P; Point E;
	E[0] = *P_shadow; E[1] = *(P_shadow + 1); E[2] = *(P_shadow + 2);
	P[0] = *P_shadow + *RayDir; P[1] = *(P_shadow + 1) + *(RayDir + 1); P[2] = *(P_shadow + 2) + *(RayDir + 2);


	Vec v0; Vec v1; Vec v2;
	v0[0] = *polyvertices;
	v0[1] = *(polyvertices + 1);
	v0[2] = *(polyvertices + 2);
	v1[0] = *(polyvertices + 3);
	v1[1] = *(polyvertices + 4);
	v1[2] = *(polyvertices + 5);
	v2[0] = *(polyvertices + 6);
	v2[1] = *(polyvertices + 7);
	v2[2] = *(polyvertices + 8);

	Vec v10; Vec v20; Vec N; Flt d; Vec R; Vec v21; Vec v02;
	Vec pvec, qvec;
	Vec tvec;
	for (int i = 0; i < 3; i++)
	{
		v10[i] = v1[i] - v0[i];
		v20[i] = v2[i] - v0[i];
		v21[i] = v2[i] - v1[i];
		v02[i] = v0[i] - v2[i];
		R[i] = P[i] - E[i];   // should it be *P,*E
		tvec[i] = E[i] - v0[i];
	}

	// pvec = R crossproduct v20

	pvec[0] = R[1] * v20[2] - R[2] * v20[1];
	pvec[1] = R[2] * v20[0] - R[0] * v20[2];
	pvec[2] = R[0] * v20[1] - R[1] * v20[0];

	Flt det;
	det = v10[0] * pvec[0] + v10[1] * pvec[1] + v10[2] * pvec[2];

	if (det < 0.05f) return false;
	if (fabs(det) < 0.05f) return false;

	Flt invDet, u, v;
	invDet = 1.0f / det;

	u = (tvec[0] * pvec[0] + tvec[1] * pvec[1] + tvec[2] * pvec[2])* invDet;
	if (u < 0 || u>1) return false;

	qvec[0] = tvec[1] * v10[2] - tvec[2] * v10[1];
	qvec[1] = tvec[2] * v10[0] - tvec[0] * v10[2];
	qvec[2] = tvec[0] * v10[1] - tvec[1] * v10[0];

	v = (R[0] * qvec[0] + R[1] * qvec[1] + R[2] * qvec[2]) *invDet;

	if (v < 0 || u + v >1) return false;

	*t = (v20[0] * qvec[0] + v20[1] * qvec[1] + v20[2] * qvec[2]) * invDet;

	return true;
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*N[0] = v10[1]*v20[2] - v10[2] * v20[1];
	N[1] = v10[2]*v20[0] - v10[0] * v20[2];
	N[2] = v10[0]*v20[1] - v10[1] * v20[0];

	d = N[0] * v0[0] + N[1] * v0[1] + N[2] * v0[2];
	
	t = (N[0] * (*E) + N[1] * (*(E+1)) + N[2] * (*(E+2)) + d)/(N[0]*R[0] +N[1]*R[1] + N[2]*R[2]);
	//cout << t << "actual";
	Vec Q; // intersection point on the triangle Q= E+tR
	for (int i = 0; i < 3; i++) {
		Q[i] = (*(E+i)) + t*(R[i]);
	}
	
	Vec vp0, vp1, vp2, C1, C2, C3; 
	for (int i = 0; i < 3; i++) {
		vp0[i] = Q[i] - v0[i];
		vp1[i] = Q[i] - v1[i];
		vp2[i] = Q[i] - v2[i];
	}

	C1[0] = v10[1]*vp0[2] - v10[2] * vp0[1];
	C1[1] = v10[2]*vp0[0] - v10[0] * vp0[2];
	C1[2] = v10[0]*vp0[1] - v10[1] * vp0[0];

	Flt c11 = N[0] * C1[0] + N[1] * C1[1] + N[2] * C1[2];

	C2[0] = v21[1] * vp1[2] - v21[2] * vp1[1];
	C2[1] = v21[2] * vp1[0] - v21[0] * vp1[2];
	C2[2] = v21[0] * vp1[1] - v21[1] * vp1[0];

	Flt c22 = N[0] * C2[0] + N[1] * C2[1] + N[2] * C2[2];

	C3[0] = v02[1] * vp2[2] - v02[2] * vp2[1];
	C3[1] = v02[2] * vp2[0] - v02[0] * vp2[2];
	C3[2] = v02[0] * vp2[1] - v02[1] * vp2[0];

	Flt c33 = N[0] * C3[0] + N[1] * C3[1] + N[2] * C3[2];

	if (c11 >= 0 && c22 >= 0 && c33 >= 0){
		return true;
	}
	*/


