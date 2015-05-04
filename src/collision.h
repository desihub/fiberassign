#ifndef COLLISION_H
#define COLLISION_H

#include	<cstdlib>
#include	<cmath>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include	<string>
#include	<vector>
#include	<algorithm>
#include	<exception>
#include	<sys/time.h>
#include        <map>
#include        "omp.h"
#include        "macros.h"
#include        "misc.h"

class PosP {
	public:
	double r1, r2, d;
	Dlist dimensions;

	PosP(double r1, double r2, double d);
};

class polygon {
	public:
	Dplist pts;
	dpair axis;
	
	void add(double a, double b);
	void transl(dpair t);
	void rotation(double t);
	void rotation_origin(double t);
};
polygon create_fh(PosP posp);

class Pols : public std::vector<polygon> {};

void rot_pt(dpair& A, dpair ax, double t);
bool intersect(dpair p1, dpair q1, dpair p2, dpair q2);
bool collision(dpair O1, dpair G1, dpair O2, dpair G2);
#endif
