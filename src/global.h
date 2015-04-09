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

// Collect -----------------------------------------------------------
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp);

void collect_available_tilefibers(Gals& G, const Plates& P);

// Assignment functions ----------------------------------------------
void simple_assign(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int next=-1);

void new_assign_fibers(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int next=-1);

void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next=-1);

void improve_from_kind(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, str kind, int next=-1);

void update_plan_from_one_obs(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int end);

void redistribute_tf(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next=-1);

// Results functions --------------------------------------------------
void results_on_inputs(str outdir, const Gals& G, const Plates& P, const Feat& F, bool latex=false);

void display_results(str outdir, const Gals& G, const Plates &P, const PP& pp, const Feat& F, const Assignment& A, bool latex=false);

#endif
