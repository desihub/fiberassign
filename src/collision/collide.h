#ifndef GLOBAL_H
#define GLOBAL_H

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
#include	"../macros.h"
#include	"../misc.h"

// Positionners parameters
class PosP {
	public:
	double r1,r2,a,b,d; // d dist between two pos

	PosP(double r1, double r2, double a, double b, double d);
};

// Positionner
class pos {
	public:
	double t, l;
	dpair c; // Theta lambda center

	pos(double t0, double l0, dpair c0);
	pos(dpair c0); // Randomly
	pos(int a, int b, PosP posp); // Randomly, with vectors a (up-right) and b (right)
	dpair armc(PosP posp) const; // Arm center
	dpair fibc(PosP posp) const; // Fib center
};

class Schema {
	public:
	Table neighbors;
	std::vector<pos> poss;

	void add(int a, int b, PosP posp);
	void compute_neigh(PosP posp);
	double dist_fib(int a, int b, PosP posp);
};

Schema hexag(PosP posp);
bool colliding(dpair c1, double r1, dpair c2, double r2);
bool colliding(pos p1, pos p2, PosP posp);
int colliding(int i, Schema sch, PosP posp);
double cartesian_weight(pos p1, pos p2);
#endif
