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

// Feat --------------------------------------------------
class Feat { // F for features
	public:
	List prio; // Priorities 
	List goal;

	Feat();
};

// PP ---------------------------------------------------
class PP { // PP for plate parameters
	public:
	std::vector<double> fp; // All fiber positions (x,y) in mm
	List spectrom; // All spectrometer assignments of fibers
	Table fibers_of_sp; // Inv map of spectrom, fibers_of_sp[sp] are fibers of spectrom sp
	Table N; // Identify neighboring positionners : neighb of fiber k are N[k]
	
	PP();
	void read_fiber_positions(const char pos_name[]);
	void get_neighbors();
	void compute_fibsofsp(); // Computes fibers_of_sp
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
	List av_gals_plate() const; // Av gals of the plate
};

class Plates : public std::vector<struct plate> {};
Plates read_plate_centers(const char center_name[], int m = 1);
List gals_range_fibers(const Plates& P); // How many galaxies in range of a fiber

// galaxy -------------------------------------------------
class galaxy {
	public:
	double nhat[3], priority;
	int id, nobs;
	double ra, dec, z;
	std::vector<pair> av_tfs; // available tile/fibers

	void print_av_tfs();
	int prio(const Feat& F) const;
};

class Gals : public std::vector<struct galaxy> {};
Gals read_galaxies(const char fname[], int n = 1);

// Assignment ---------------------------------------------
// 2 mappings of assignments : (j,k) -> id(gal) ; id(gal)[5] -> (j,k)
class Assignment {
	public:
	Table TF; // TF for tile fiber, #tiles X #fibers
	std::vector<std::vector<pair> > PG; // PG for pass galaxy, #passes X #galaxies
	int na; // Number of assignments
	/*List Nobs; // Nobs[g] : number of times left g has to be observed*/

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
	int unused_f(int j);
	int unused_fbp(int j, int k, const PP& pp); // Number of unassigned fibers of the petal corresponding to (j,k)
	bool once_obs(int g); // Already observed ?
};
// Collect -----------------------------------------------------
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, Time &t);

void collect_available_tilefibers(Gals& G, const Plates& P);

// Assignment "global" ----------------------------------------
void assign_fibers(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A);

void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute(const Gals& G, const Plates&P, const PP& pp, Assignment& A);
// Assignment "for one" -----------------------------------------
void assign_fibers_for_one(int j, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A);

void improve_for_one(int j, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute_for_one(int j, const Gals& G, const Plates&P, const PP& pp, Assignment& A);

// Other available functions ------------------------------------
Table conflicts(const Gals& G, const Plates& P, const PP& pp, const Assignment& A);
void display_results(const Gals& G, const Plates &P, const Feat& F, const Assignment& A);
void plot_freefibers(std::string s, const Plates&P, const Assignment& A);
void print_free_fibers(const PP& pp, const Assignment& A);
void time_line(const Gals& G, const Feat& F, const Assignment& A);
#endif
