#include    <cstdlib>
#include    <cmath>
#include    <fstream>
#include    <sstream>
#include    <iostream>
#include    <iomanip>
#include    <string>
#include    <vector>
#include    <algorithm>
#include    <exception>
#include        <map>
#include        "misc.h"
#include        "collision.h"

// Angles are dpair (cos t, sin t)

// PosP
PosP::PosP(double r10, double r20) {
    r1 = r10; r2 = r20;
}

// Intersection of segments

int orientation(const dpair& p, const dpair& q, const dpair& r) {
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
 
bool intersect(const dpair& p1, const dpair& q1, const dpair& p2, const dpair& q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    if (o1 == o2) return false; // (when error) terminate called after throwing an instance of 'std::bad_alloc'
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    if (o3 == o4) return false;
    else return true;
}

// Intersection of segment and circle

bool intersect_seg_circ(const dpair& A, const dpair& B, const dpair& O, const double& rad) {
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
element::element() { color = 'k'; }

element::element(bool b) { is_seg = b; color = 'k'; }

element::element(const dpair& Ocenter, const double& rad0) {
    is_seg = false;
    O = Ocenter;
    rad = rad0;
    color = 'k';
}

element::element(const dpair& A, const dpair& B) {
    is_seg = true;
    add(A);
    add(B);
    color = 'k';
}

element::element(const dpair& A, char c, double transp, double rad0) {
    is_seg = true;
    add(A);
    color = c;
    transparency = transp;
    radplot = rad0;
}

void element::add(const double& a, const double& b) {
    segs.push_back(dpair(a,b));
}

void element::add(const dpair& p) {
    segs.push_back(p);
}

void element::transl(const dpair& t) {
    if (is_seg) for (int i=0; i<segs.size(); i++) segs[i] = segs[i]+t;
    else O = O + t;
}

void element::rotation(const dpair& t, const dpair& axis) {
    if (is_seg) for (int i=0; i<segs.size(); i++) rot_pt(segs[i],axis,t);
    else rot_pt(O,axis,t);
}

void element::print() const {
    printf("seg ? %d",is_seg); fl();
    if (!O.isnull() || rad!=0) { printf(" - center (%f,%f) rad %f",O.f,O.s,rad); fl(); }
    printf(" sizesegs= %lu ",segs.size()); fl();
    if (0<segs.size()) {
        debl(" - segs : ");
        for (int i=0; i<segs.size(); i++) segs[i].print();
    }
    else printf(" no seg ");
    printf("\n"); fl();
}

void element::limits(Dlist& lims) const {
    if (is_seg) for (int i=0; i<segs.size(); i++) {
        if (segs[i].f < lims[0]) lims[0] = segs[i].f;
        if (segs[i].s < lims[2]) lims[2] = segs[i].s;
        if (segs[i].f > lims[1]) lims[1] = segs[i].f;
        if (segs[i].s > lims[3]) lims[3] = segs[i].s;
    }
    else {
        if (O.f - rad < lims[0]) lims[0] = O.f - rad;
        if (O.s - rad < lims[2]) lims[2] = O.s - rad;
        if (O.f + rad > lims[1]) lims[1] = O.f + rad;
        if (O.s + rad > lims[3]) lims[3] = O.s + rad;
    }
}

bool intersect(const element& e1, const element& e2) {
    if (e1.is_seg && e2.is_seg) {
        for (int i=0; i<e1.segs.size()-1; i++) {
            for (int j=0; j<e2.segs.size()-1; j++) {
                bool b = intersect(e1.segs[i],e1.segs[i+1],e2.segs[j],e2.segs[j+1]);
                if (b) return true;
            }
        }
        return false;
    }
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
    if (!e1.is_seg && !e2.is_seg) return (sq(e1.O,e2.O)<sq(e1.rad+e2.rad));
}

// polygon
void polygon::add(const element& el) {
    elmts.push_back(el);
}

void polygon::add(const polygon& p) {
    for (int i=0; i<p.elmts.size(); i++) elmts.push_back(p.elmts[i]);
}

void polygon::transl(const dpair& t) {
    for (int i=0; i<elmts.size(); i++) elmts[i].transl(t);
    axis = axis + t;
}

void polygon::rotation(const dpair& t) {
    for (int i=0; i<elmts.size(); i++) elmts[i].rotation(t,axis);
}

void polygon::rotation_origin(const dpair& t) {
    dpair origin = dpair();
    for (int i=0; i<elmts.size(); i++) elmts[i].rotation(t,origin);
    rot_pt(axis,origin,t);
}

void polygon::print() const {
    for (int i=0; i<elmts.size(); i++) elmts[i].print();
    printf(" - axis : "); axis.print();
}

void polygon::set_color(char c) {
    for (int i=0; i<elmts.size(); i++) elmts[i].color = c;
}

Dlist polygon::limits() const {
    Dlist lims = initDlist(4);
    lims[0] = 1e4;
    lims[2] = 1e4;
    lims[1] = -1e4;
    lims[3] = -1e4;
    for (int i=0; i<elmts.size(); i++) elmts[i].limits(lims);
    return lims;
}

polygon create_fh() {
    polygon fh;
    fh.axis = dpair(-3.0,0);

    // Only segments
    //element el;
    //el.add(-4.240,0.514); el.add(-3.514,1.240); el.add(-2.668,1.240); el.add(-2.235,0.990); el.add(0.387,0.990); el.add(0.967,0.410); el.add(0.967,-0.410); el.add(0.387,-0.990); el.add(-1.844,-0.990); el.add(-1.981,-1.757); el.add(-2.688,-2.015); el.add(-2.944,-1.922); el.add(-2.944,-1.339); el.add(-3.682,-1.072); el.add(-4.240,-0.514); el.add(-4.240,0.514);
    //fh.add(el);

    // Segments and disks
    fh.add(element(dpair(),0.967)); // Head disk
    fh.add(element(dpair(-3.0,0.990),dpair(0,0.990))); // Upper segment
    element set(true); // Lower segment
    set.add(-2.944,-1.339); set.add(-2.944,-2.015); set.add(-1.981,-1.757); set.add(-1.844,-0.990); set.add(0,-0.990);
    fh.add(set);

    fh.transl(dpair(6.0,0));
    return fh;
}

polygon create_cb() {
    polygon cb;
    cb.axis = dpair(3.0,0);

    // Segments and disks
    cb.add(element(dpair(3.0,0),2.095));
    element set(true);
    set.add(5.095,-0.474); set.add(4.358,-2.5); set.add(2.771,-2.5); set.add(1.759,-2.792); set.add(0.905,-0.356);
    cb.add(set);
    return cb;
}

// Functions
void rot_pt(dpair& A, const dpair& ax, const dpair& angle) {
    dpair P = A-ax;
    if (P.isnull()) return;
    double r = norm(P);
    dpair cos_sin = cos_sin_angle(P);
    dpair sum = sum_angles(cos_sin,angle);
    A = ax+dpair(r*sum.f,r*sum.s);
}

Dlist angles(const dpair& A, const PosP& posp) {
    double phi = 2*acos(norm(A)/(2*posp.r1));
    double theta = atan(A.s/A.f)-phi/2 + (A.f<0 ? M_PI : 0);
    Dlist ang; ang.push_back(cos(theta)); ang.push_back(sin(theta)); ang.push_back(cos(phi)); ang.push_back(sin(phi));
    return ang;
}

void repos_cb_fh(polygon& cb, polygon& fh, const dpair& O, const dpair& G, const PosP& posp) {
    //repositions positioner (central body, ferule holder)
    dpair Gp = G-O;
    if (sq(posp.r1+posp.r2)<sq(Gp)) printf("Error galaxy out of range of fiber in repos_fiber O.f %f O.s %f G.f %f  G.s %f\n",O.f,O.s,G.f,G.s );
    Dlist anglesT = angles(Gp,posp);
    dpair theta_ = dpair(anglesT[0],anglesT[1]);
    dpair phi_ = dpair(anglesT[2],anglesT[3]);
    cb.rotation_origin(theta_);
    fh.rotation_origin(theta_);
    fh.rotation(phi_);
    fh.transl(O);
    cb.transl(O);
}

bool collision(const polygon& p1, const polygon& p2) {
    for (int i=0; i<p1.elmts.size(); i++) {
        for (int j=0; j<p2.elmts.size(); j++) {
            if (intersect(p1.elmts[i],p2.elmts[j])) return true;
        }
    }
    return false;
}
