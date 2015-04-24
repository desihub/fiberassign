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
#include        <stdlib.h>     /* srand, rand */
#include	"../macros.h"
#include	"../misc.h"
#include	"collide.h"

PosP::PosP(double r10, double r20, double a0, double b0, double d0) {
	r1 = r10; r2 = r20; a = a0; b = b0; d = d0;
}

pos::pos(double t0, double l0, dpair c0) {
	t = t0; l = l0 ; c = c0;
}

pos::pos(dpair c0) {
	t = fRand(0,pi2);
	l = fRand(0,pi2);
	c = c0;
}

pos::pos(int a, int b, PosP posp) {
	t = fRand(0,pi2);
	l = fRand(0,pi2);
	c = dpair((a/2+b)*posp.d,a*posp.d*sqrt(3)/2);
}

dpair pos::armc(PosP posp) const {
	return sum(cartesian(posp.a,t),c);
}

dpair pos::fibc(PosP posp) const {
	return sum(armc(posp),cartesian(posp.b,l));
}

void Schema::add(int a, int b, PosP posp) {
	poss.push_back(pos(a,b,posp));
}

void Schema::compute_neigh(PosP posp) {
	Table T;
	for (int i=0; i<poss.size(); i++) {
		List l;
		for (int j=0; j<poss.size(); j++) {
			if (i!=j && dist(poss[i].c,poss[j].c)<=posp.d+0.01) l.push_back(j);
		}
		T.push_back(l);
	}
	neighbors = T;
}

Schema hexag(PosP posp) {
	Schema sch;
	sch.add(0,0,posp);
	sch.add(1,0,posp);
	sch.add(-1,0,posp);
	sch.add(0,1,posp);
	sch.add(0,-1,posp);
	sch.add(-1,1,posp);
	sch.add(1,-1,posp);
	sch.neighbors = initTable(7,0);
	int l[] = {1,2,3,4,5,6};
	sch.neighbors[0] = initList(l,6);
	return sch;
}

bool colliding(dpair c1, double r1, dpair c2, double r2) {
	return (r1+r2>dist(c1,c2));
}

bool colliding(pos p1, pos p2, PosP posp) {
	return (colliding(p1.armc(posp),posp.r1,p2.fibc(posp),posp.r2) || colliding(p1.fibc(posp),posp.r1,p2.armc(posp),posp.r2) || colliding(p1.fibc(posp),posp.r1,p2.fibc(posp),posp.r2));
}

int colliding(int i, Schema sch, PosP posp) {
	for (int j=0; j<sch.neighbors.size(); j++) {
		if (i!=j && colliding(sch.poss[i],sch.poss[j],posp)) return j;
	}
	return -1;
}
