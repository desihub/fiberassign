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

int Nplate; int MaxSS; int MaxSF; double PlateRadius; double Collide; double NeighborRad; double PatrolRad; double TotalArea; int Ngal; int MaxPrio; int MaxObs; int Categories; int Npass; int Nfiber; int MinUnused; int Npetal; int Nfbp;

int main(int argc, char **argv) {
	//// Initialization ---------------------------------------------------
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
	init_time_at(time,"# Begin collecting available galaxies",t);
	collect_galaxies_for_all(G,T,P,pp,t);
	print_time(time,"# ... took :");

	// For each galaxy, computes available tilefibers
	init_time_at(time,"# Begin computing available tilefibers",t);
	collect_available_tilefibers(G,P);
	print_time(time,"# ... took :");

	//results_on_inputs(G,P,F);

	//// Assignment plate per plate -----------------------------------
	Assignment A(G,F);
	// Order : verificate that passes are increasing
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=0; jj<2000; jj++) {
		printf(" - Plate %d : ",jj);
		int j = A.order[jj];
		assign_fibers_for_one(j,G,P,pp,F,A,true);
		improve2(1,"SF",G,P,pp,F,A,true);
		improve2(1,"SS",G,P,pp,F,A,true);
		printf("  %s not assigned\n",f(Nfiber-A.na(j,1)).c_str());
		// here is the real observation moment
		A.update_once_obs(j);
		A.update_nobsv_tmp_for_one(j);
		A.next_plate++;
		if (A.unused_f(j)>500) {}
	}
	print_time(time,"# ... took :");
	display_results(G,P,pp,F,A,false);
	
	A.verif(P);

	init_time_at(time,"# Begin global assignment ------------",t);
	assign_fibers(3000,G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve",t);
	improve(3000,G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve SF",t);
	improve2(3000,"SF",G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve SS",t);
	improve2(3000,"SS",G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	A.verif(P);
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=2000; jj<5000; jj++) {
		printf(" - Plate %d : ",jj); fl();
		//A.verif(P);
		int j = A.order[jj];
		// here is observation time, we have then access to nobsv 
		printf("  %s not assigned\n",f(Nfiber-A.na(j,1)).c_str());
		A.update_once_obs(j);
		update_plan_from_one_obs(4999,G,P,pp,F,A);
		A.update_nobsv_tmp_for_one(j);
		A.next_plate++;
	}
	print_time(time,"# ... took :");
	display_results(G,P,pp,F,A,false);


	init_time_at(time,"# Begin global assignment ------------",t);
	assign_fibers(-1,G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve",t);
	improve(-1,G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve SF",t);
	improve2(-1,"SF",G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin improve SS",t);
	improve2(-1,"SS",G,P,pp,F,A,true);
	print_time(time,"# ... took :");

	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=5000; jj<Nplate; jj++) {
		printf(" - Plate %d : ",jj); fl();
		int j = A.order[jj];
		// here is observation time
		printf("  %s not assigned\n",f(Nfiber-A.na(j,1)).c_str());
		A.update_once_obs(j);
		update_plan_from_one_obs(Nplate-1,G,P,pp,F,A);
		A.update_nobsv_tmp_for_one(j);
		A.next_plate++;
	}
	print_time(time,"# ... took :");

	display_results(G,P,pp,F,A,false,true);
	print_free_fibers(G,pp,F,A,false);

	//for (int i=1; i<=10; i++) {
		//print_table("Petal "+i2s(i),A.infos_petal(1000*i,5,G,P,pp,F));
	//}
	
	//plot_freefibers("free_fiber_plot.txt",P,A);   
	//time_line(G,F,A);

	print_time(t,"# Finished !... in");
	return(0);
}
