#ifndef GLOBAL_H
#define GLOBAL_H

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
#include        "structs.h"

// Collect -----------------------------------------------------
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, Time &t);

void collect_available_tilefibers(Gals& G, const Plates& P);

// Assignment "global" ----------------------------------------
// Assign next "next" fibers
void assign_fibers(int next, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, bool tmp=false);

void improve(int next, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, bool tmp=false);
void improve2(int next, str kind, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, bool tmp=false);

void redistribute_tf(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute_g(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute_g_by_kind(str kind, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

// Assignment "for one" ---------------------------------------
void assign_fibers_for_one(int j, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, bool tmp=false);

void update_plan_from_one_obs(int end, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

// Other useful functions -------------------------------------
Table conflicts(const Gals& G, const Plates& P, const PP& pp, const Assignment& A);
void results_on_inputs(const Gals& G, const Plates& P, const Feat& F);
void display_results(const Gals& G, const Plates &P, const PP& pp, const Feat& F, const Assignment& A, bool latex=false, bool tmp=false);
void print_free_fibers(const Gals& G, const PP& pp, const Feat& F, const Assignment& A, bool latex=false);
void plot_freefibers(str s, const Plates&P, const Assignment& A);
#endif
