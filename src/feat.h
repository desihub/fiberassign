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

	int InterPlate; 
	bool Randomize; 
	bool Pacman; 
	int Npass;
	int MaxSS;
	int MaxSF;
	double PlateRadius;
	double NeighborRad;
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
	// Indirectly set by features file
	int Categories;
	Smap ids; // Redundtant but optimizes (inv of kind)

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

	// Methods
	Feat();
	void readInputFile(const char fname[]);
	int id(str s) const;
	int maxgoal(int kind) const; // Gives max goal for a galaxy of this kind (goal(Ly-a) for all QSO for example)
	List maxgoal() const; // List (function of kind) of max goals according to the kind (defined by prio)
	void init_ids();
};
#endif
