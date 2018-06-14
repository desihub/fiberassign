#ifndef COLLISION_H
#define COLLISION_H

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <map>
#include "misc.h"

// Can contain information on geometry of fiber one day (hardcoded for now) not
// used yet
class PosP {
    public:
        double r1, r2;

        PosP (double r1, double r2);
};

// Rotation of A (angle t) aroud axis axis
void rot_pt (dpair & A, const dpair & ax, const dpair & angle);

// Intersection between segments p1-q1 and p2-q2
bool intersect (const dpair & p1, const dpair & q1, const dpair & p2, const
                dpair & q2);

// List segments, or circles
class element {
    public:

        // True if segments, false if circle
        bool is_seg;

        // Segments
        Dplist segs;

        // Circle center
        dpair O;

        // Circle radius
        double rad;

        // For python plot
        char color;

        // For plotting
        double transparency;

        // For plotting
        double radplot;

        element ();

        // b : is_seg
        element (bool b);

        // Creates circle
        element (const dpair & O, const double & rad);

        // Creates segment list with only the segment AB
        element (const dpair & A, const dpair & B);

        // Point with color
        element (const dpair & A, char c, double transp, double rad0);

        // Add point to segments list
        void add (const double & a, const double & b);

        // Add point to segments list
        void add (const dpair & p);

        // Translation by t
        void transl (const dpair & t);

        // Rotation of the angle def by t (cos,sin) around axis axis
        void rotation (const dpair & t, const dpair & axis);

        void print () const;

        // Xmin, Xmax, Ymin, Ymax of the element
        void limits (Dlist & lims) const;
};

class Elements : public std::vector <element> {};

// Intersection between 2 elements
bool intersect (const element & e1, const element & e2);

class polygon {
    public:
        // Segments and circles
        Elements elmts;
        dpair axis;

        void add (const element & el);

        // Add all elements of p
        void add (const polygon & p);

        // Translation
        void transl (const dpair & t);

        // Rotation aroud polygon's axis
        void rotation (const dpair & t);

        // Rotation around O
        void rotation_origin (const dpair & t);

        void print () const;

        // Change color of the polygon
        void set_color (char c);

        // Xmin, Xmax, Ymin, Ymax of the polygon
        Dlist limits () const;
};

// fh = ferrule holder
polygon create_fh ();

// cb = central body
polygon create_cb ();

// Check collision
bool collision (const polygon & p1, const polygon & p2);

// cos sin theta and phi for a galaxy which is in A (ref to the origin)
Dlist angles (const dpair & A, const PosP & posp);

// G loc of galaxy, O origin. Repositions fh and cb in the good position
// according to O, G
void repos_cb_fh (polygon & cb, polygon & fh, const dpair & O, const dpair & G,
                  const PosP & posp);

#endif
