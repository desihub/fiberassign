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
	std::vector<str> kind;

	Feat();
	int id(str s) const;
	int maxgoal(int kind) const; // Gives max goal for a galaxy of this kind (goal(Ly-a) for all QSO for example)
	List maxgoal() const;
};

// PP ---------------------------------------------------
class PP { // PP for plate parameters
	public:
	std::vector<double> fp; // All fiber positions (x,y) in mm
	List spectrom; // All spectrometer assignments of fibers
	Table fibers_of_sp; // Inv map of spectrom, fibers_of_sp[sp] are fibers of spectrom sp, redundant but optimizes
	Table N; // Identify neighboring positionners : neighb of fiber k are N[k]
	
	PP();
	void read_fiber_positions(const char pos_name[]);
	void get_neighbors();
	void compute_fibsofsp(); // Computes fibers_of_sp
	List fibs_of_same_pet(int k) const;
};

// galaxy -------------------------------------------------
class galaxy {
	public:
	int id;
	double nhat[3];
	double ra, dec, z;
	std::vector<pair> av_tfs; // available tile/fibers

	void print_av_tfs();
	int prio(const Feat& F) const;
	str kind(const Feat& F) const;
};

class Gals : public std::vector<struct galaxy> {};
Gals read_galaxies(const char fname[], int n = 1);

// Plate -------------------------------------------------
struct onplate { // The position of a galaxy in plate coordinates
	int id;
	double pos[2];
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
List gals_range_fibers(const Plates& P);
List av_gals_of_kind(int kind, int j, int k, const Gals& G, const Plates& P, const Feat& F);

// Assignment ---------------------------------------------
// 2 mappings of assignments : (j,k) -> id(gal) ; id(gal)[5] -> (j,k)
class Assignment {
	public:
	Table TF; // TF for tile fiber, #tiles X #fibers
	std::vector<std::vector<pair> > PG; // PG for pass galaxy, #passes X #galaxies
	Cube kinds; // Cube[j][sp][id] : number of fibers of spectrometer sp and plate j that have the kind id. Redundtant information, but optimizes computations
	List probas; // Number of galaxies of this kind
	List plates_done; // List of done plates
	List nobsv; // List of nobs, redundant but optimizes
	List nobsv_tmp; // List of nobs, redundant but optimizes

	Assignment(const Gals& G, const Feat& F);
	~Assignment();
	void assign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp);
	void unassign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp);
	bool is_assigned_pg(int ip, int g) const;
	bool is_assigned_tf(int j, int k) const; 
	int na(); // Number of assignments
	int nobs(int g, const Gals& G, const Feat& F, bool b=false) const; // Counts how many time ob should be observed else more. If b=true, return Maxobs for this kind
	std::vector<pair> chosen_tfs(int g) const; // Pairs (j,k) chosen by g in the Npass passes
	List unused_f() const;
	int unused_f(int j);
	Table unused_fbp(const PP& pp) const; // Unused fibers  by petal
	int unused_fbp(int j, int k, const PP& pp) const; // Number of unassigned fibers of the petal corresponding to (j,k)
	Table used_by_kind(str kind, const Gals& G, const PP& pp, const Feat& F) const;
	bool once_obs(int g); // Already observed ?
	int nkind(int j, int k, str kind, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const; // Number of fibers assigned to the kind "kind" on the petal of (j,k)
	double get_proba(int i, const Gals& G, const Feat& F); // p(fake QSO | QSO) for example
	Table infos_petal(int j, int pet, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const;
	List fibs_of_kind(str kind, int j, int pet, const Gals& G, const PP& pp, const Feat& F) const; // Sublist of fibers assigned to a galaxy of type kind for this petal
	List fibs_unassigned(int j, int pet, const Gals& G, const PP& pp, const Feat& F) const;
	void update_nobsv_tmp();
};

// Global functions -----------------------------------------

void writeTFfile(const Plates& P, std::ofstream TFfile);
void writeGfile(Gals& G, std::ofstream Gfile);
void readGfile(Gals& G, galaxy Gtemp, std::ifstream GXfile);

#endif
