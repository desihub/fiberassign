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

// PosP

PosP::PosP(double r10, double r20, double d0) {
	r1 = r10; r2 = r20; d = d0;
}

void PosP::get_Xdimensions() {
	dimension.push_back(-0.967);
	dimension.push_back(-0.387);
	dimension.push_back(1.844);
	dimension.push_back(1.981);
	dimension.push_back(2.235);
	dimension.push_back(2.668);
	dimension.push_back(2.688);
	dimension.push_back(2.944);
	dimension.push_back(3.0);
	dimension.push_back(3.514);
	dimension.push_back(3.682);
	dimension.push_back(4.240);
}

void PosP::get_Ydimensions() {
	dimension.push_back(1);
}

// polygon
void polygon::add(double a, double b) {
	pts.push_back(dpair(a,b));
}

void polygon::transl(dpair t) {
	for (int i=0; i<pts.size(); i++) pts[i] = pts[i]+t;
	axis = axis + t;
}

void polygon::rotation(double t) {
	for (int i=0; i<pts.size(); i++) {
		rot_pt(pts[i],axis,t);
	}
}

void polygon::rotation_origin(double t) {
	for (int i=0; i<pts.size(); i++) {
		rot_pt(pts[i],dpair(),t);
	}
	rot_pt(axis,dpair(),t);
}

// Functions

void rot_pt(dpair& A, dpair ax, double t) {
	dpair po = polar(A-ax);
	po.s += t;
	A = cartesian(po);		
}

dpair angles(dpair A, PosP posp) { // t,p for a galaxy that is in A
	double tmp = acos(norm(A)/(2*posp.r1));
	return dpair(atan(A.s/A.f)-tmp,2*tmp);
}

polygon repos_fiber(polygon& fh, dpair O, dpair G, PosP posp) { // G loc of galaxy, O origin
	dpair Gp = G-O;
	if (posp.r1+posp.r2<norm(Gp)) {
		printf("Error galaxy out of range of fiber in repos_fiber\n");
		return fh;
	}
	dpair tp = angles(Gp,posp);
	fh.rotation_origin(tp.f);
	fh.rotation(tp.s);
	fh.transl(O);
}

bool find_collision(const polygon& p1, const polygon& p2) {
}

polygon create_polygon(double a) {
	polygon pol;
	pol.add(0,0);
	pol.add(0,a);
	pol.add(a,a);
	pol.add(a,0);
	pol.axis = dpair(0,0);
}

polygon create_fh(PosP posp) { // fh = ferrule holder
	polygon fh;
	fh.axis = dpair(-posp.r1,0);

	fh.add(-4.240,0.514);
	fh.add(-3.514,1.240);
	fh.add(-2.668,1.240);
	fh.add(-2.235,0.990);
	fh.add(0.387,0.990);
	fh.add(0.967,0.410);
	fh.add(0.967,-0.410);
	fh.add(0.387,-0.990);
	fh.add(-1.844,-0.990);
	fh.add(-1.981,-1.757);
	fh.add(-2.688,-2.015);
	fh.add(-2.944,-1.922);
	fh.add(-2.944,-1.339);
	fh.add(-3.682,-1.072);
	fh.add(-4.240,-0.514);

	fh.transl(dpair(posp.r2+posp.r1,0));
	return fh;
}
