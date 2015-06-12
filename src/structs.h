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
#include        <map>
#include        "misc.h"
#include        "feat.h"
#include        "collision.h"

// PP ---------------------------------------------------
class PP { // PP for plate parameters
	public:
	Dlist fp; // All fiber positions (x,y) in mm
	List spectrom; // All spectrometer assignments of fibers
	Table fibers_of_sp; // Inv map of spectrom, fibers_of_sp[sp] are fibers of spectrom sp, redundant but optimizes
	Table N; // Identify neighboring positionners : neighb of fiber k are N[k]
	
	PP();
	void read_fiber_positions(const Feat& F);
	void get_neighbors(const Feat& F);
	void compute_fibsofsp(const Feat& F); // Computes fibers_of_sp
	List fibs_of_same_pet(int k) const;
	dpair coords(int k) const; // Coords of fiber k
};

// galaxy -------------------------------------------------
class galaxy {
	public:
	int id;
	double nhat[3];
	double ra, dec, z;
	Plist av_tfs; // available tile/fibers

	void print_av_tfs();
	str kind(const Feat& F) const;
};
class Gals : public std::vector<struct galaxy> {};

Gals read_galaxies(const Feat& F);

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
	List density; // density[k] is the ponderated number of objects available to (j,k)


	List av_gals_plate(const Feat& F) const; // Av gals of the plate
};
class Plates : public std::vector<struct plate> {};

Plates read_plate_centers(const Feat& F);
List av_gals_of_kind(int kind, int j, int k, const Gals& G, const Plates& P, const Feat& F);

// Assignment ---------------------------------------------
// 2 mappings of assignments : (j,k) -> id(gal) ; id(gal)[5] -> (j,k)
class Assignment {
	public:
	//// ----- Members
	Table TF; // TF for tile fiber, #tiles X #fibers TF[j][k] is the chosen galaxy, -1 if not yet chosen
	List order; // Order of tiles we want to assign, only 1-n in simple increasing order for the moment
	int next_plate; // Next plate in the order

	// Redundant information (optimizes computation time)
	Ptable GL; // GL for galaxy - list : #galaxies X (variable) #chosen TF: gives chosen tf's for galaxy g
	Cube kinds; // Cube[j][sp][id] : number of fibers of spectrometer sp and plate j that have the kind id
	Table unused; // Table [j][p] giving number of unused fibers on this petal
	List nobsv; // List of nobs, redundant but optimizes, originally true goal
	List nobsv_tmp; // List of nobs, redundant but optimizes, apparent goal, i.e. goal of category of this type, gets updated
	List once_obs; // 0 if not observed, 1 if observed  [list of all galaxies]


	//// ----- Methods
	Assignment(const Gals& G, const Feat& F);
	~Assignment();
	void assign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp);
	void unassign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp);
	int find_collision(int j, int k, int g, const PP& pp, const Gals& G, const Plates& P, const Feat& F, int col=-1) const;
	bool find_collision(int j, int k, int kn, int g, int gn, const PP& pp, const Gals& G, const Plates& P, const Feat& F, int col=-1) const;
	int is_collision(int j, int k, const PP& pp, const Gals& G, const Plates& P, const Feat& F) const;
	void verif(const Plates& P, const Gals& G, const PP& pp, const Feat& F) const; // Verif mappings are right
	int is_assigned_jg(int j, int g) const;
	int is_assigned_jg(int j, int g, const Gals& G, const Feat& F) const;
	bool is_assigned_tf(int j, int k) const; 
	int na(const Feat& F, int begin=0, int size=-1) const; // Number of assignments (changes) within plates begin to begin+size
	int nobs(int g, const Gals& G, const Feat& F, bool tmp=true) const; // Counts how many more times object should be observed. If tmp=true, return maximum for this kind (temporary information)
	//if tmp=false we actually know the true type from the start
	Plist chosen_tfs(int g, const Feat& F, int begin=0, int size=-1) const; // Pairs (j,k) chosen by g, amongst size plates from begin
	int nkind(int j, int k, int kind, const Gals& G, const Plates& P, const PP& pp, const Feat& F, bool pet=false) const; // Number of fibers assigned to the kind "kind" on the petal of (j,k). If pet=true, we don't take k but the petal p directly instead
	List fibs_of_kind(int kind, int j, int pet, const Gals& G, const PP& pp, const Feat& F) const; // Sublist of fibers assigned to a galaxy of type kind for (j,p)
	List sort_fibs_dens(int j, const List& fibs, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const; // Sort this list of fibers by decreasing density
	List fibs_unassigned(int j, int pet, const Gals& G, const PP& pp, const Feat& F) const; // Subist of unassigned fibers for (j,p)

	// Update information
	void update_nobsv_tmp_for_one(int j, const Feat& F);
	void update_once_obs(int j, const Feat& F);

	// Used to compute results at the end
	Table infos_petal(int j, int pet, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const;
	List unused_f(const Feat& F) const;
	Table unused_fbp(const PP& pp, const Feat& F) const; // Unused fibers by petal
		float colrate(const PP& pp, const Gals& G, const Plates& P, const Feat& F, int j=-1) const; // Get collision rate, j = plate number
	int nobs_time(int g, int j, const Gals& G, const Feat& F) const; // Know the number of remaining observations of g when the program is at the tile j, for pyplotTile

	// Not used (but could be useful)
	int unused_f(int j, const Feat& F) const; // Number of unused fiber on the j'th plate
	int unused_fbp(int j, int k, const PP& pp, const Feat& F) const; // Number of unassigned fibers of the petal corresponding to (j,k)

	void update_nobsv_tmp(const Feat& F);
};

bool collision(dpair O1, dpair G1, dpair O2, dpair G2, const Feat& F); // collisions from  looking at galaxy G1 with fiber positionner centered at 01 and etc calculated in mm on plate

int fprio(int g, const Gals& G, const Feat& F, const Assignment& A);//priority of galaxy g

double plate_dist(const double theta);//plate scale conversion
struct onplate change_coords(const struct galaxy& O, const struct plate& P);
dpair projection(int g, int j, const Gals& G, const Plates& P); // Projection of g on j
int num_av_gals(int j, int k, const Gals& G, const Plates& P, const Feat& F, const Assignment& A); // weighted (and only with remaining observation according to the moment in the survey), and doesn't take into account other kinds than QSO LRG ELG not used

// Pyplot -----------------------------------------------
class pyplot {
	public:
	polygon pol;
	Slist text;
	Dplist textpos;

	pyplot(polygon p);
	void addtext(dpair p, str s);
	void plot_tile(str directory, int j, const Feat& F) const;
};


#endif
