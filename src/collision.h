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
#include        <map>
#include        "misc.h"
#include        "feat.h"

class PosP {
	public:
	double r1, r2;
	Dlist dimensions;

	PosP(double r1, double r2);
};

class element { // List segments, or circles
	public:
	bool is_seg; // True if segments, false if circle
	Dplist segs;
	dpair O; // Circle center
	double rad; // Circle radius

	element();
	element(bool b); // b : is_seg
	element(dpair O, double rad); // Creates circle
	element(dpair A, dpair B); // Creates segment list with only the segment AB
	void add(double a, double b); // Add points to segments list
	void add(dpair p); // Add points to segments list
	void transl(dpair t); // Translation
	void rotation(dpair t, dpair axis); // Rotation around axis
};
class Elements : public std::vector<element> {};
bool intersect(element e1, element e2); // Intersection between 2 elements

class polygon {
	public:
	Elements elmts;
	dpair axis;
	
	void add(element el);
	void transl(dpair t); // Translation
	void rotation(dpair t); // Rotation aroud polygon's axis
	void rotation_origin(dpair t); // Rotation around 0
};
polygon create_fh(PosP posp); // fh = ferrule holder
polygon create_cb(PosP posp); // cb = central body

class Pols : public std::vector<polygon> {};

void rot_pt(dpair& A, dpair axis, dpair t); // Rotation of A (angle t) aroud axis
bool intersect(dpair p1, dpair q1, dpair p2, dpair q2); // Intersection between segments p1-q1 and p2-q2
bool collision(dpair O1, dpair G1, dpair O2, dpair G2, const Feat& F); // Intersection of fh looking at galaxy G1 with fiber positionner centered in 01 and ...
#endif
