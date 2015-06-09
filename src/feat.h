#ifndef FEAT_H
#define FEAT_H

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
#include        "collision.h"

// Feat --------------------------------------------------
class Feat { // F for features
	public:
	// Set by features file
	str galFile;
	str tileFile;
	str fibFile;
	str outDir;

	List prio; // Priorities 
	List priopost; // Priorities when we know the kind
	List goal;
	Slist kind;
	Slist type;

	int InterPlate; 
	bool Randomize; 
	bool Pacman; 
	int Npass;
	int MaxSS;
	int MaxSF;
	double PlateRadius;
	int Analysis;
	bool InfDens;
	double TotalArea;
	double invFibArea;
	int moduloGal;
	int moduloFiber;
	
	bool Collision; // True when we take collisions into account
	bool Exact;
	double AvCollide;
	double Collide;
	double NoCollide;
	double PatrolRad;
	double NeighborRad;

	bool PlotObsTime;
	bool PlotHistLya;
	bool PlotDistLya;
	bool PlotFreeFibHist;
	bool PlotFreeFibTime;
	bool PlotSeenDens;
	bool Verif;
	// Indirectly set by features file
	int Categories;
	Slist types; // Same as type but with only QSO LRG ELG SS SF
	Smap ids; // Redundtant but optimizes (inv of kind)
	Smap ids_types; // Same but for types (inv of types) (1 to QSO, ...)

	// Set after reading other input files
	int Nplate;
	int Ngal;
	int Nfiber;
	int Npetal;
	int Nfbp; // Number of fibers by petals

	// Memorizes geometry of central body and fiber holder
	polygon cb;
	polygon fh;

	// Permit to count some things
	int Count;

	// Practical used lists
	List no_ss_sf;
	List ss_sf;

	// Methods
	Feat();
	void readInputFile(const char fname[]);
	int id(str s) const;
	int maxgoal(int kind) const; // Gives max goal for a galaxy of this kind (goal(Ly-a) for all QSO for example)
	List maxgoal() const; // List (function of kind) of max goals according to the kind (defined by prio)
	void init_ids(); // Init ids member
	void init_ids_types(); // Same
	List init_ids_list(str l[], int size) const;
	void init_types();
	void init_no_ss_sf();
	void init_ss_sf();
	bool iftype(int kind, str type) const; // Whether kind is of type "type"
};
#endif
