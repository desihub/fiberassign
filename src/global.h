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
void assign_fibers(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A);

void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);
// Assignment "for one" -----------------------------------------
void assign_fibers_for_one(int j, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A);

void improve_for_one(int j, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

void redistribute_for_one(int j, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A);

// Other useful functions ------------------------------------
List gals_range_fibers(const Plates& P); // How many galaxies in range of a fiber
void results_on_inputs(const Gals& G, const Plates& P, const Feat& F);
void display_results(const Gals& G, const Plates &P, const Feat& F, const Assignment& A);
Table conflicts(const Gals& G, const Plates& P, const PP& pp, const Assignment& A);
void plot_freefibers(std::string s, const Plates&P, const Assignment& A);
void print_free_fibers(const PP& pp, const Assignment& A);
void time_line(const Gals& G, const Feat& F, const Assignment& A);
#endif
