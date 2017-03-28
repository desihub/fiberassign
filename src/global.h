#ifndef GLOBAL_H
#define GLOBAL_H

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
#include    <sys/time.h>
#include    <fitsio.h>
#include    "misc.h"
#include    "feat.h"
#include    "structs.h"


// Collect -----------------------------------------------------------
// From the HTM tree, collects for each fiber and for each plate, the available galaxies
void collect_galaxies_for_all(const MTL& M, const htmTree<struct target>& T, Plates& P, const FP& pp, const Feat& F);

// From the previous computations, computes the inverse map, that is to say, for each target, computes the tile-fibers which can reach it
void collect_available_tilefibers(MTL& M, const Plates& P, const Feat& F);

// Assignment functions ----------------------------------------------
// First simple assignment plan, executing find_best on every plate on every fiber
void simple_assign(MTL& M, Plates& P, const FP& pp, const Feat& F, Assignment& A);

// More fine first assignment plan, 

void improve(MTL& M, Plates&P, const FP& pp, const Feat& F, Assignment& A, int jstart);

void improve_from_kind(MTL& M, const Plates&P, const FP& pp, const Feat& F, Assignment& A, str kind, int next=-1);

void update_plan_from_one_obs(int j, const Gals& G, MTL& M, Plates&P, const FP& pp, const Feat& F, Assignment& A);

//void redistribute_tf(MTL& M, Plates&P, const FP& pp, const Feat& F, Assignment& A, int next=-1);
void redistribute_tf(MTL& M, Plates&P, const FP& pp, const Feat& F, Assignment& A, int jstart);
void assign_sf_ss(int j, MTL& M, Plates& P, const FP& pp, const Feat& F, Assignment& A);

void assign_unused(int j, MTL& M,  Plates& P, const FP& pp, const Feat& F, Assignment& A);

void diagnostic(const Gals& G, Feat& F, const Assignment& A);
// Results functions --------------------------------------------------

void display_results(str outdir, const Gals& G, const MTL& M, const Plates &P, const FP& pp, Feat& F, const Assignment& A, int last_tile, bool latex=false);

void diagnostic(const MTL& M, const Gals& G, Feat& F, const Assignment& A);

void write_FAtile_ascii(int j, str outdir, const MTL& M, const Plates& P, const FP& pp, const Feat& F, const Assignment& A);


void fa_write(int j, str outdir, const MTL& M, const Plates& P, const FP& pp, const Feat& F, const Assignment& A);

void write_save_av_gals(int j, str outdir, const MTL& M, const Plates& P, const FP& pp, const Feat& F);

void pyplotTile(int j, str fname, const Gals& Secret, const MTL& M, const Plates& P, const FP& pp, const Feat& F, const Assignment& A);

void overlappingTiles(str fname, const Feat& F, const Assignment& A);
#endif
bool myorder(int i, int j);


List sort_by_fiber_number(FP & pp,std::vector<int> v);
Dlist sort_doublet_by_fiber_number(FP & pp, Dlist v);
