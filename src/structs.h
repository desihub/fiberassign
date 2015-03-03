#ifndef STRUCTS_H
#define STRUCTS_H

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
#include        "omp.h"
#include        "macros.h"
#include        "misc.h"

// PP ---------------------------------------------------
class PP { // PP for plate parameters
	public:
	std::vector<double> fp; // All fiber positions (x,y) in mm
	List spectrom; // All spectrometer assignments of fibers
	Table N; // Identify neighboring positionners : neighb of fiber k are N[k]
	
	PP();
	void read_fiber_positions(const char pos_name[]);
	void get_neighbors();
};

// Plate -------------------------------------------------
struct onplate { // The position of a galaxy in plate coordinates
	double pos[2];
	int id;
};
class Onplates : public std::vector<struct onplate> {};

class plate {
	public:
	int idp;
	double nhat[3]; // Unit vector pointing to plate
	int ipass; // Pass
	Table av_gals; // av_gals[k] : available galaxies of fiber k

	void print_plate() const;
};

class Plates : public std::vector<struct plate> {};
Plates read_plate_centers(const char center_name[], int m = 1);

// galaxy -------------------------------------------------
class galaxy {
	public:
	double nhat[3], priority;
	int id, nobs;
	double ra, dec, z;
	std::vector<pair> av_tfs; // available tile/fibers

	void print_av_tfs();
};

class Gals : public std::vector<struct galaxy> {};
Gals read_galaxies(const char fname[], int n = 1);
void collect_available_tilefibers(Gals& G, const Plates& P);

// Assignment ---------------------------------------------
// 2 mappings of assignments : (j,k) -> id(gal) ; id(gal)[5] -> (j,k)
class Assignment {
	public:
	Table TF; // TF for tile fiber, #tiles X #fibers
	std::vector<std::vector<pair> > PG; // PG for pass galaxy, #passes X #galaxies
	int na; // Number of assignments

	Assignment();
	~Assignment();
	void assign(int j, int k, int g, const Plates& P);
	void unassign(int j, int k, int g, const Plates& P);
	bool is_assigned_pg(int ip, int g) const;
	bool is_assigned_tf(int j, int k) const; 
	std::vector<pair> chosen_tfs(int g) const; // Pairs (j,k) chosen by g in the Npass passes
	List unused_fibers() const;
	Table unused_fibers_by_petal(const PP& pp) const;
	List hist_petal(int binsize_petal, const PP& pp) const;
};

// Available functions for main ----------------------------
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, Time &t);

void assign_fibers(const Gals& G, const Plates& P, const PP& pp, Assignment& A);

Table conflicts(const Gals& G, const Plates& P, const PP& pp, const Assignment& A);

void display_results(const Gals& G, const List& goal, const Plates &P, const Assignment& A);

void improve(const Gals& G, const Plates&P, const PP& pp, Assignment& A);

void redistribute(const Gals& G, const Plates&P, const PP& pp, Assignment& A);

void print_free_fibers(const PP& pp, const Assignment& A);

void plot_freefibers(std::string s, const Plates&P, const Assignment& A);

void time_line(const Gals& G, const List& goal, const Assignment& A);

#endif
