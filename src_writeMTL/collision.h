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

// Can contain information on geometry of fiber one day (hardcoded for now) not used yet
class PosP {
	public:
	double r1, r2;

	PosP(double r1, double r2);
};

void rot_pt(dpair& A, const dpair& ax, const dpair& angle); // Rotation of A (angle t) aroud axis axis
bool intersect(const dpair& p1, const dpair& q1, const dpair& p2, const dpair& q2); // Intersection between segments p1-q1 and p2-q2

class element { // List segments, or circles
	public:
	bool is_seg; // True if segments, false if circle
	Dplist segs; // Segments
	dpair O; // Circle center
	double rad; // Circle radius
	char color; // For python plot
	double transparency; // For plotting
	double radplot; // For plotting

	element();
	element(bool b); // b : is_seg
	element(const dpair& O, const double& rad); // Creates circle
	element(const dpair& A, const dpair& B); // Creates segment list with only the segment AB
	element(const dpair& A, char c, double transp, double rad0); // Point with color
	void add(const double& a, const double& b); // Add point to segments list
	void add(const dpair& p); // Add point to segments list
	void transl(const dpair& t); // Translation by t
	void rotation(const dpair& t, const dpair& axis); // Rotation of the angle def by t (cos,sin) around axis axis
	void print() const;
	void limits(Dlist& lims) const; // Xmin, Xmax, Ymin, Ymax of the element
};
class Elements : public std::vector<element> {};
bool intersect(const element& e1, const element& e2); // Intersection between 2 elements

class polygon {
	public:
	Elements elmts; // Segments and circles
	dpair axis;
	
	void add(const element& el);
	void add(const polygon& p); // Add all elements of p
	void transl(const dpair& t); // Translation
	void rotation(const dpair& t); // Rotation aroud polygon's axis
	void rotation_origin(const dpair& t); // Rotation around O
	void print() const;
	void set_color(char c); // Change color of the polygon
	Dlist limits() const; // Xmin, Xmax, Ymin, Ymax of the polygon
};
polygon create_fh(); // fh = ferrule holder
polygon create_cb(); // cb = central body

bool collision(const polygon& p1, const polygon& p2); // Check collision
Dlist angles(const dpair& A, const PosP& posp); // cos sin theta and phi for a galaxy which is in A (ref to the origin)
void repos_cb_fh(polygon& cb, polygon& fh, const dpair& O, const dpair& G, const PosP& posp); // G loc of galaxy, O origin. Repositions fh and cb in the good position according to O, G

#endif
