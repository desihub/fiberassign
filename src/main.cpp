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
#include	"modules/htmTree.h"
#include	"modules/kdTree.h"
#include        "macros.h"
#include        "misc.h"
#include        "structs.h"
#include        "global.h"

int Nplate; int MaxSS; int MaxSF; double PlateRadius; double Collide; double NeighborRad; double PatrolRad; double TotalArea; int Ngal; int MaxPrio; int MaxObs; int Categories; int Npass; int Nfiber; int MinUnused; int Npetal; int Nfbp; int InterPlate; bool Randomize;


int main(int argc, char **argv) {
	//// Initializations -------------------------------------------------
	srand48(1234); // Make sure we have reproducability
	check_args(argc);
	Time t,time; // t for global, time for local
	init_time(t);

	// Parameters
	Npass = 5;
	MinUnused = 50;
	MaxSS = 10; MaxSF = 40;
	PlateRadius = 1.65;
	TotalArea = 15789.0;
	Collide = 2.1;
	NeighborRad = 11.0;
	PatrolRad = 6.0;
	InterPlate = 200;
	Randomize = true;

	str kind[] = {"QSO Ly-a","QSO Tracer","LRG","ELG","Fake QSO","Fake LRG","SS","SF"};
	int prio[] = {1,1,3,5,1,3,2,4}; // has to be >= 0
	int goal[] = {5,1,2,1,1,1,1,1};

	// Taking previous features into account
	Feat F; F.kind = initList(kind); F.prio = initList(prio); F.goal = initList(goal); MaxPrio = max(F.prio); MaxObs = max(F.goal); Categories = F.kind.size();

	// Read galaxies
	Gals G;
	G = read_galaxies(argv[1],atoi(argv[5])); // ! Reads all galaxies iff arg2=1
	Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(Ngal).c_str(),argv[1]);

	// Read fiber center positions and compute things
	PP pp;
	pp.read_fiber_positions(argv[3]); 
	Nfiber = pp.fp.size()/2;
	Npetal = max(pp.spectrom)+1; // spectrom has to be identified from 0 to Npetal-1
	Nfbp = (int) (Nfiber/Npetal);
	pp.get_neighbors(); 
	pp.compute_fibsofsp();

	// Read plates
	Plates P = read_plate_centers(argv[2],1); // ! Reads all galaxies iff arg2=1
	Nplate = P.size();
	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(Nplate).c_str(),argv[2],Nfiber,argv[3]);

	//// Collect available galaxies <-> tilefibers ----------------------
	// HTM Tree of galaxies
	const double MinTreeSize = 0.01; // <--- ?
	init_time_at(time,"# Start building HTM tree",t);
	htmTree<struct galaxy> T(G,MinTreeSize);
	print_time(time,"# ... took :");
	//T.stats();
	
	// For plates/fibers, collect available galaxies
	collect_galaxies_for_all(G,T,P,pp);

	// For each galaxy, computes available tilefibers
	collect_available_tilefibers(G,P);

	//results_on_inputs(G,P,F,true);

	//// Assignment ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	Assignment A(G,F);
	print_time(t,"# Start assignment at : ");
	// Make a plan ----------------------------------------------------
	//new_assign_fibers(G,P,pp,F,A,-1,true);
	simple_assign(G,P,pp,F,A,2000);
	improve(G,P,pp,F,A,2000);
	improve_from_kind(G,P,pp,F,A,"SF",2000);
	improve_from_kind(G,P,pp,F,A,"SS",2000);
	//improve(G,P,pp,F,A,2000); // Improves almost nothing
	//improve_from_kind(G,P,pp,F,A,"SF",2000);
	//improve_from_kind(G,P,pp,F,A,"SS",2000);
	A.verif(P);

	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=0; jj<2000; jj++) {
		int j = A.next_plate;
		printf(" - Plate %d : ",j); fl();
		//improve_from_kind(G,P,pp,F,A,"SF",1);
		//improve_from_kind(G,P,pp,F,A,"SS",1);
		// here is observation time
		printf("  %s not assigned  - ",format(5,f(Nfiber-A.na(j,1))).c_str());
		update_plan_from_one_obs(G,P,pp,F,A,1999);
		A.next_plate++;
	}
	print_time(time,"# ... took :");

	// Make a plan ----------------------------------------------------
	//new_assign_fibers(G,P,pp,F,A,-1,true);
	simple_assign(G,P,pp,F,A);
	improve(G,P,pp,F,A);
	improve_from_kind(G,P,pp,F,A,"SF");
	improve_from_kind(G,P,pp,F,A,"SS");
	//improve(G,P,pp,F,A);
	//improve_from_kind(G,P,pp,F,A,"SF");
	//improve_from_kind(G,P,pp,F,A,"SS");
	A.verif(P);

	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=2000; jj<Nplate; jj++) {
		int j = A.next_plate;
		printf(" - Plate %d : ",j); fl();
		//improve_from_kind(G,P,pp,F,A,"SF",1);
		//improve_from_kind(G,P,pp,F,A,"SS",1);
		// here is observation time
		printf("  %s not assigned  - ",format(5,f(Nfiber-A.na(j,1))).c_str());
		update_plan_from_one_obs(G,P,pp,F,A,Nplate-1);
		A.next_plate++;
	}
	print_time(time,"# ... took :");

	// Results ----------------------------------------------------
	A.verif(P);
	display_results(G,P,pp,F,A,true);
	print_free_fibers(G,pp,F,A,true);

	print_time(t,"# Finished !... in");
	return(0);
}
