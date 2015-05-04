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
#include        "collision.h"

// PosP

PosP::PosP(double r10, double r20, double d0) {
	r1 = r10; r2 = r20; d = d0;
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

void repos_fh(polygon& fh, dpair O, dpair G, PosP posp) { // G loc of galaxy, O origin
	dpair Gp = G-O;
	if (posp.r1+posp.r2<norm(Gp)) {
		printf("Error galaxy out of range of fiber in repos_fiber\n");
	}
	dpair tp = angles(Gp,posp);
	fh.rotation_origin(tp.f);
	fh.rotation(tp.s);
	fh.transl(O);
}

bool collision(polygon& p1, polygon& p2) {
	for (int i=0; i<p1.pts.size()-1; i++) {
		for (int j=0; j<p2.pts.size()-1; j++) {
			bool b = intersect(p1.pts[i],p1.pts[i+1],p2.pts[j],p2.pts[j+1]);
			if (b) return true;
		}
	}
	return false;
}

bool collision(dpair O1, dpair G1, dpair O2, dpair G2) {
	PosP posp(3,3,10);
	polygon p1 = create_fh(posp);
	polygon p2 = create_fh(posp);
	repos_fh(p1,O1,G1,posp);
	repos_fh(p2,O2,G2,posp);
	return collision(p1,p2);
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
	fh.add(-4.240,0.514);

	fh.transl(dpair(posp.r2+posp.r1,0));
	return fh;
}

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(dpair p, dpair q, dpair r)
{
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    double val = (q.s- p.s) * (r.f - q.f) - (q.f - p.f) * (r.s - q.s);
    if (val == 0) return 0;  // colinear
    return (val > 0)? 1: 2; // clock or counterclock wise
}
 
bool intersect(dpair p1, dpair q1, dpair p2, dpair q2)
{
    // Find the four orientations needed for general case
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    if (o1 == o2) return false;
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    if (o3 == o4) return false;
    else return true;
}
