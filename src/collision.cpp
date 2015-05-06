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

// Angles are dpair (cos t, sin t)

// PosP
PosP::PosP(double r10, double r20, double d0) {
	r1 = r10; r2 = r20; d = d0;
}

// Intersection of segments

int orientation(dpair p, dpair q, dpair r) {
    // See 10th slides from following link for derivation of the formula
    // http://www.dcs.gla.ac.uk/~pat/52233/slides/Geometry1x1.pdf
    double val = (q.s- p.s) * (r.f - q.f) - (q.f - p.f) * (r.s - q.s);
    /*
    if (val == 0) {
	    printf("CAUTIOUS in orientation, colinear case ! \n");
	    return 0;  // colinear
    }
    */
    return (val > 0) ? 1 : 2; // clock or counterclock wise
}
 
bool intersect(dpair p1, dpair q1, dpair p2, dpair q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    if (o1 == o2) return false;
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    if (o3 == o4) return false;
    else return true;
}

// Intersection of segment and circle

bool intersect_seg_circ(dpair A, dpair B, dpair O, double rad) {
	// Make a drawing, with C the projection of O on (AB), use scalar products
	double rad_sq = sq(rad);

	double AO_sq = sq(A,O);
	double AB_AO = scalar_prod(A,B,O);
	if (AB_AO<=0) return AO_sq<rad_sq;

	double BO_sq = sq(B,O);
	double BA_BO = scalar_prod(B,A,O);
	if (BA_BO<=0) return BO_sq<rad_sq;

	double AB_sq = sq(A,B);
	return AO_sq*(1-sq(AB_AO)/(AO_sq*AB_sq))<rad_sq;
}

// element
element::element() {}

element::element(bool b) { is_seg = b; }

element::element(dpair Ocenter, double rad_) {
	is_seg = false;
	O = Ocenter;
	rad = rad_;
}

element::element(dpair A, dpair B) {
	is_seg = true;
	add(A);
	add(B);
}

void element::add(double a, double b) {
	segs.push_back(dpair(a,b));
}

void element::add(dpair p) {
	segs.push_back(p);
}

void element::transl(dpair t) {
	if (is_seg) for (int i=0; i<segs.size(); i++) segs[i] = segs[i]+t;
	else O = O + t;
}

void element::rotation(dpair t, dpair axis) {
	if (is_seg) for (int i=0; i<segs.size(); i++) rot_pt(segs[i],axis,t);
	else rot_pt(O,axis,t);
}

bool intersect(element e1, element e2) {
	if (e1.is_seg && e2.is_seg) {
		for (int i=0; i<e1.segs.size()-1; i++) {
			for (int j=0; j<e2.segs.size()-1; j++) {
				bool b = intersect(e1.segs[i],e1.segs[i+1],e2.segs[j],e2.segs[j+1]);
				if (b) return true;
			}
		}
		return false;
	}
	if (!e1.is_seg && !e2.is_seg) return (sq(e1.O,e2.O)<sq(e1.rad+e2.rad));
	if (e1.is_seg && !e2.is_seg) {
		for (int i=0; i<e1.segs.size()-1; i++) {
			return intersect_seg_circ(e1.segs[i],e1.segs[i+1],e2.O,e2.rad);
		}
	}
	if (!e1.is_seg && e2.is_seg) {
		for (int i=0; i<e2.segs.size()-1; i++) {
			return intersect_seg_circ(e2.segs[i],e2.segs[i+1],e1.O,e1.rad);
		}
	}
}

// polygon
void polygon::add(element el) {
	elmts.push_back(el);
}

void polygon::transl(dpair t) {
	for (int i=0; i<elmts.size(); i++) elmts[i].transl(t);
	axis = axis + t;
}

void polygon::rotation(dpair t) {
	for (int i=0; i<elmts.size(); i++) {
		elmts[i].rotation(t,axis);
	}
}

void polygon::rotation_origin(dpair t) {
	dpair origin = dpair();
	for (int i=0; i<elmts.size(); i++) {
		elmts[i].rotation(t,origin);
	}
	rot_pt(axis,origin,t);
}

polygon create_fh(PosP posp) {
	polygon fh;
	fh.axis = dpair(-posp.r1,0);

	// Only segments
	//element el;
	//el.add(-4.240,0.514); el.add(-3.514,1.240); el.add(-2.668,1.240); el.add(-2.235,0.990); el.add(0.387,0.990); el.add(0.967,0.410); el.add(0.967,-0.410); el.add(0.387,-0.990); el.add(-1.844,-0.990); el.add(-1.981,-1.757); el.add(-2.688,-2.015); el.add(-2.944,-1.922); el.add(-2.944,-1.339); el.add(-3.682,-1.072); el.add(-4.240,-0.514); el.add(-4.240,0.514);
	//fh.add(el);

	// Segments and disks
	fh.add(element(dpair(-3.0,0),0.990)); // Big disk
	fh.add(element(dpair(),0.967)); // Small disk
	fh.add(element(dpair(-3.0,0.990),dpair(0,0.990))); // Upper segment
	element set(true); // Lower segment
	set.add(-2.944,-1.339); set.add(-2.944,-2.015); set.add(-1.981,-1.757); set.add(-1.844,-0.990);
	fh.add(set);

	fh.transl(dpair(posp.r2+posp.r1,0));
	return fh;
}

polygon create_cb(PosP posp) {
	polygon fh;
	fh.axis = dpair(3.0,0);

	// Segments and disks
	fh.add(element(dpair(3.0,0),3.990));
	element set(true);
	set.add(5.095,-0.474); set.add(4.358,-2.5); set.add(2.771,-2.5); set.add(1.759,-2.792); set.add(0.905,0.356);
	fh.add(set);
	return fh;
}

// Functions
void rot_pt(dpair& A, dpair ax, dpair angle) {
	dpair P = A-ax;
	double r = norm(P);
	dpair cos_sin = cos_sin_angle(P);
	dpair sum = sum_angles(cos_sin,angle);

	A = dpair(r*sum.f,r*sum.s);		
}

Dlist angles(dpair A, PosP posp) { // cos sin theta and phi for a galaxy which is in A (ref to the origin)
	double phi = 2*(acos(norm(A)/(2*posp.r1)));
	double theta = atan(A.s/A.f)-phi/2;
	Dlist ang; ang.push_back(cos(theta)); ang.push_back(sin(theta)); ang.push_back(cos(phi)); ang.push_back(sin(phi));
	return ang;
}

void repos_cb_fh(polygon& cb, polygon& fh, dpair O, dpair G, PosP posp) { // G loc of galaxy, O origin. Repositions fh and cb for 0, G
	dpair Gp = G-O;
	if (sq(posp.r1+posp.r2)<sq(Gp)) printf("Error galaxy out of range of fiber in repos_fiber\n");
	Dlist anglesT = angles(Gp,posp);
	dpair theta_ = dpair(anglesT[0],anglesT[1]);
	dpair phi_ = dpair(anglesT[2],anglesT[3]);
	cb.rotation_origin(theta_);
	fh.rotation(sum_angles(theta_,phi_));
	fh.transl(O);
	cb.transl(O);
}

bool collision(polygon& p1, polygon& p2) {
	for (int i=0; i<p1.elmts.size(); i++) {
		for (int j=0; j<p2.elmts.size(); j++) {
			if (intersect(p1.elmts[i],p2.elmts[i])) return true;
		}
	}
	return false;
}

bool collision(dpair O1, dpair G1, dpair O2, dpair G2) {
	double dist_sq = sq(G1,G2);
	if (dist_sq < sq(Collide)) return true;
	if (dist_sq > sq(NoCollide)) return false;
	Count++;
	PosP posp(3,3,10);
	polygon fh1 = create_fh(posp);
	polygon fh2 = create_fh(posp);
	polygon cb1 = create_cb(posp);
	polygon cb2 = create_cb(posp);
	repos_cb_fh(cb1,fh1,O1,G1,posp);
	repos_cb_fh(cb2,fh2,O2,G2,posp);
	if (collision(fh1,fh2)) return true;
	if (collision(cb1,fh2)) return true;
	if (collision(cb2,fh1)) return true;
	return false;
}
