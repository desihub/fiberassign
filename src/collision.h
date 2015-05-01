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
	double r1, r2, d;
	Dlist dimensions;

	PosP(double r1, double r2, double d);
	void get_Xdimensions();
	void get_Ydimensions();
}

class polygon {
	public:
	Dplist pts;
	dpair axis;
	
	void add(double a, double b);
	void transl(dpair t);
	void rotation(double t);
	void rotation_origin(double t);
};

class Pols : public std::vector<polygon> {};

#endif
