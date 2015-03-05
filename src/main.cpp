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

int MaxPlate;
int Ngal;

int main(int argc, char **argv) {
	//// Initialization ---------------------------------------------------
	check_args(argc);
	srand48(1234); // Make sure we have reproducability
	Time t,time; // t for global, time for local
	init_time(t);
	// Number of observations desired, and priorities
	Feat F;
	F.prio[0]=2; F.prio[1]=2; F.prio[2]=3; F.prio[3]=4; F.prio[4]=2; F.prio[5]=3;
	F.goal.assign(6,1); F.goal[0]=5; F.goal[2]=2;
	F.kind[0]="QSO Ly-a"; F.kind[1]="QSO Tracer"; F.kind[2]="LRG"; F.kind[3]="ELG"; F.kind[4]="Fake QSO"; F.kind[5]="Fake LRG";

	// Galaxies
	Gals G;
	G = read_galaxies(argv[1],1); // ! Reads all galaxies iff arg2=1
	Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(Ngal),argv[1]);

	// Read fiber center positions and compute things
	PP pp;
	pp.read_fiber_positions(argv[3]); 
	pp.get_neighbors(); 
	pp.compute_fibsofsp(); 

	// Read plates
	Plates P = read_plate_centers(argv[2]); // ! Reads all galaxies iff arg2=1
	MaxPlate = P.size();
	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(MaxPlate),argv[2],MaxFiber,argv[3]);

	//// Collect available galaxies <-> tilefibers ----------------------
	// HTM Tree of galaxies (doing it inside collect fun doesn't optimize)
	const double MinTreeSize = 0.01; // <--- ?
	init_time_at(time,"# Start building HTM tree",t);
	htmTree<struct galaxy> T(G,MinTreeSize);
	print_time(time,"# ... took :");
	//T.stats();
	// For plates/fibers, collect available galaxies
	init_time_at(time,"# Begin collecting available galaxies",t);
	collect_galaxies_for_all(G,T,P,pp,t);
	print_time(time,"# ... took :");
	//print_list("\n How many galaxies in range of a fiber",gals_range_fibers(P));
	// For each galaxy, computes available tilefibers
	init_time_at(time,"# Begin computing available tilefibers",t);
	collect_available_tilefibers(G,P);
	print_time(time,"# ... took :");

	results_on_inputs(G,P,F);
	//P[0].print_plate();

	////// Assignment -----------------------------------------------------
	//Assignment A;
	//Naive assignment
		//init_time_at(time,"# Begin naive assignment",t);
	//assign_fibers(G,P,pp,F,A);
	//print_time(time,"# ... took");

	//Table Q = conflicts(G,P,pp,A);
	//display_results(G,P,F,A);

	//// Improve --------------------------------------------------------
	//init_time_at(time,"# Begin improve",t);
	//improve(G,P,pp,F,A);
	//print_time(time,"# ... took :");
	//display_results(G,P,F,A);
	//print_free_fibers(pp,A);

	//init_time_at(time,"# Begin redistribute",t);
	//redistribute(G,P,pp,F,A);
	//print_time(time,"# ... took :"); 
	//print_free_fibers(pp,A);

	//init_time_at(time,"# Begin redistribute",t);
	//redistribute(G,P,pp,F,A);
	//print_time(time,"# ... took :"); 
	//print_free_fibers(pp,A);


	//plot_freefibers("free_fiber_plot.txt",P,A);   
	//time_line(G,F,A);

	Assignment A0;
	// Order : verificate that passes are increasing
	init_time_at(time,"# Begin real time assignment",t);
	for (int j=0; j<MaxPlate && j<1000; j++) {
		printf(" - Plate %d : ",j);
		assign_fibers_for_one(j,G,P,pp,F,A0);
		if (A0.unused_f(j)>MinUnused*MaxPetal) {
			//assign_fibers_for_one(j,G,P,pp,F,A0);
			//improve_for_one(j,G,P,pp,F,A0);
			//redistribute_for_one(j,G,P,pp,F,A0);
			//assign_fibers_for_one(j,G,P,pp,F,A0);
		}
		printf("\n");
	}
	print_time(time,"# ... took :");

	display_results(G,P,F,A0);
	print_free_fibers(pp,A0);	

	print_time(t,"# Finished !... in");
	return(0);
}
