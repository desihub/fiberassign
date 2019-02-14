#ifndef STRUCTS_H
#define STRUCTS_H

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <sys/time.h>
#include <map>
#include "feat.h"

// fiber-positioner  -----------------------------------

class fpos {
    public:
        // number of fiber(!)  not positioner, read from fibpos file
        int fib_num;
        // x position in mm of positioner
        double fp_x;
        // y position in mm of positioner
        double fp_y;
        // location of fiber
        int location;
        // which spectrometer 0 - 9
        int spectrom;
        // 0=OK
        int stuck = 0;
        // 0=OK
        int broken = 0;
        // Identify neighboring positioners : neighbors of fiber k are N[k]
        std::vector <int> N;
        dpair coords;
};

class FP  : public std::vector <struct fpos> {
    public:
        // Inv map of spectrom, fibers_of_sp[sp] are fibers of spectrom sp,
        // redundant
        Table fibers_of_sp;
};

FP  read_fiber_positions (const Feat & F);

void read_fiber_status (FP & FibPos, const Feat & F);

bool int_pairCompare (const std::pair <int, int> & firstElem,
                      const std::pair <int, int> & secondElem);

// galaxy truth -------------------------------------------

class galaxy {
    public:
        // the unique identifier
        long long targetid;

        // the true type when used with a secret file
        long category;
        double z;
        uint16_t obsconditions;
};

class Gals : public std::vector <struct galaxy> {};

Gals read_Secretfile (str filename, const Feat & F);

Gals read_Secretfile_ascii (str filename, const Feat & F);

std::vector <int> count_galaxies (const Gals & G);

// target -----------------------------------------------------

class target {
    public:
        // note that this is how we read it in! not long
        long long id;
        int nobs_remain, nobs_done;
        double nhat[3];
        double ra, dec, subpriority;
        long desi_target, mws_target, bgs_target;
        int SS, SF, lastpass, priority_class, t_priority, once_obs;
        char brickname [9];
        Plist av_tfs;
        // 16bit mask indicating under what conditions this target can be
        // observed.
        uint16_t obsconditions;
};

class MTL : public std::vector <struct target> {
    public:
        std::vector <int> priority_list;
};

MTL read_MTLfile (str filename, const Feat & F, long SS, long SF);

void assign_priority_class (MTL & M);

// Plate -------------------------------------------------

struct onplate {
    // The position of a galaxy in plate coordinates
    int id;
    double pos[2];
};

class Onplates : public std::vector <struct onplate> {};

class plate {
    public:
        int tileid;
        double tilera;
        double tiledec;
        // Unit vector pointing to plate
        double nhat[3];
        // Pass
        int ipass;
        // av_gals[k] : available galaxies of fiber k
        Table av_gals;
        // density[k] is the weighted number of objects available to (j,k)
        List density;
        // SS_av_gal[p] are available standard stars on petal p of this plate
        Table SS_av_gal;
        Table SF_av_gal;
        Table SS_av_gal_fiber;
        Table SF_av_gal_fiber;
        // number of SS assigned to a petal in this plate
        std::vector <int> SS_in_petal;
        std::vector <int> SF_in_petal;
        // true if tile has some galaxies within reach
        bool is_used;
        // mask defining the kind of program (DARK, BRIGHT, GRAY)
        uint16_t obsconditions;
        // Av gals of the plate
        List av_gals_plate (const Feat & F, const MTL & M, const FP & pp) const;
};

class Plates : public std::vector <struct plate> {};

Plates read_plate_centers (const Feat & F);

// Assignment ---------------------------------------------

// 2 mappings of assignments : (j,k) -> id(gal) ; id(gal)[5] -> (j,k)
class Assignment {
    public:
        //// ----- Members
        // TF for tile fiber, #tiles X #fibers TF[j][k] is the chosen galaxy,
        // -1 if not yet chosen
        Table TF;

        // Order of tiles we want to assign, only 1-n in simple increasing
        // order for the moment
        List order;

        // Order of tiles actually containing targets. size is F.Nplate
        List suborder;

        // inverse of suborder unused tiles map to -1
        List inv_order;

        // Next plate in the order, i.e suborder(next_plate) is actually next
        // plate
        int next_plate;

        // Redundant information (optimizes computation time)
        // GL for galaxy - list : #galaxies X (variable) #chosen TF: gives
        // chosen tf's for galaxy g
        Ptable GL;

        // Cube[j][sp][id] : number of fibers of spectrometer sp and plate j
        // that have the kind id
        Cube kinds;

        // Table [j][p] giving number of unused fibers on this petal
        Table unused;

        // List of nobs, redundant but optimizes, originally true goal
        // List nobsv;
        // List of nobs, redundant but optimizes, apparent goal, i.e. goal of
        // category of this type, gets updated
        // List nobsv_tmp;

        //// ----- Methods

        Assignment (const MTL & M, const Feat & F);

        ~Assignment ();

        void assign (int j, int k, int g,  MTL & M, Plates & P, const FP & pp);

        void unassign (int j, int k, int g, MTL & M, Plates & P,
                       const FP & pp);

        int find_collision (int j, int k, int g, const FP & pp, const MTL & M,
                            const Plates & P, const Feat & F,
                            int col = -1) const;

        bool find_collision (int j, int k, int kn, int g, int gn,
                             const FP & pp, const MTL & M, const Plates & P,
                             const Feat & F, int col = -1) const;

        int is_collision (int j, int k, const FP & pp, const MTL & M,
                          const Plates & P, const Feat & F) const;

        // Verif mappings are right
        void verif (const Plates & P, const MTL & M, const FP & pp,
                    const Feat & F) const;

        int is_assigned_jg (int j, int g) const;

        int is_assigned_jg (int j, int g, const MTL & M, const Feat & F) const;

        bool is_assigned_tf (int j, int k) const;

        // Number of assignments (changes) within plates begin to begin+size
        int na (const Feat & F, int begin = 0, int size = -1) const;

        // Counts how many more times object should be observed. If tmp=true,
        // return maximum for this kind (temporary information)
        // if tmp=false we actually know the true type from the start
        int nobs (int g, const MTL & M, const Feat & F, bool tmp = true) const;

        // Pairs (j,k) chosen by g, amongst size plates from begin
        Plist chosen_tfs (int g, const Feat & F, int begin = 0) const;

        // Number of fibers assigned to the kind "kind" on the petal of (j,k).
        // If pet=true, we don't take k but the petal p directly instead
        int nkind (int j, int k, int kind, const MTL & M, const Plates & P,
                   const FP & pp, const Feat & F, bool pet = false) const;

        // Sublist of fibers assigned to a galaxy of type kind for (j,p)
        List fibs_of_kind (int kind, int j, int pet, const MTL & M,
                           const FP & pp, const Feat & F) const;

        // not used
        // Sort this list of fibers by decreasing density
        List sort_fibs_dens (int j, const List & fibs, const MTL & M,
                             const Plates & P, const FP & pp,
                             const Feat & F) const;

        // Sublist of unassigned fibers for (j,p)
        List fibs_unassigned (int j, int pet, const MTL & M, const FP & pp,
                              const Feat & F) const;

        // Used to compute results at the end
        // gives total number of unused fibers
        List unused_f (const Feat & F) const;

        // Unused fibers by petal
        Table unused_fbp (const FP & pp, const Feat & F) const;

        // Get collision rate, j = plate number
        float colrate (const FP & pp, const MTL & M, const Plates & P,
                       const Feat & F, int j = -1) const;

        // Know the number of remaining observations of g when the program is
        // at the tile j, for pyplotTile
        int nobs_time (int g, int j, const Gals & Secret, const MTL & M,
                       const Feat & F) const;

        // Number of unused fiber on the j'th plate
        int unused_f (int j, const Feat & F) const;

        // not used
        // Number of unassigned fibers of the petal corresponding to (j,k)
        int unused_fbp (int j, int k, const FP & pp, const Feat & F) const;

        // not used
        void update_nobsv_tmp (const Feat & F);
};

// collisions from  looking at galaxy G1 with fiber positioner centered at 01
// and etc calculated in mm on plate
bool collision (dpair O1, dpair G1, dpair O2, dpair G2, const Feat & F);

// int fprio(int g, const Gals& G, const Feat& F, const Assignment&
// A);//priority of galaxy g
double plate_angle (double r_plate);

// plate scale conversion
double plate_dist (const double theta);

void xy2radec (double * ra, double * dec, double telra, double teldec,
               double x, double y);

struct onplate change_coords (const struct target & O, const struct plate & P);

struct onplate radec2xy (const struct target & O, const struct plate & P);

// Projection of g on j
dpair projection (int g, int j, const MTL & M, const Plates & P);

// weighted (and only with remaining observation according to the moment in the
// survey), and doesn't take into account other kinds than QSO LRG ELG not used
int num_av_gals (int j, int k, const MTL & M, const Plates & P, const Feat & F,
                 const Assignment & A);

int A_less_than_B (int year_A, int month_A, int day_A, int year_B, int month_B,
                   int day_B);

// Pyplot -----------------------------------------------

class pyplot {
    public:
        polygon pol;
        Slist text;
        Dplist textpos;
        pyplot (polygon p);
        void addtext (dpair p, str s);
        void plot_tile (str directory, int j, const Feat & F) const;
};


#endif
