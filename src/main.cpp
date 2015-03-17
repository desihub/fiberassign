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
	//std::cout << (false || true && true || false) << std::endl;
	//std::cout << (false || true && true || false || false || true) << std::endl;
	//// Initialization ---------------------------------------------------
	check_args(argc);
	srand48(1234); // Make sure we have reproducability
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
	G = read_galaxies(argv[1],1); // ! Reads all galaxies iff arg2=1
	Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(Ngal),argv[1]);

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
	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(Nplate),argv[2],Nfiber,argv[3]);

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

	////// Assignment plate per plate -----------------------------------
	//Assignment A0(G,F);
	//// Order : verificate that passes are increasing
	//init_time_at(time,"# Begin real time assignment",t);
	//for (int j=0; j<Nplate; j++) {
		//printf(" - Plate %d : ",j); std::cout.flush();
		//assign_fibers_for_one(j,G,P,pp,F,A0);
		////A0.plates_done.push_back(j);
		////A0.update_probas(G,F);
		//if (A0.unused_f(j)>500) {
			////improve_for_one(j,G,P,pp,F,A0);
			////redistribute_for_one(j,G,P,pp,F,A0);
			////assign_fibers_for_one(j,G,P,pp,F,A0);
		//}
		//std::cout << std::endl;
	//}
	//print_time(time,"# ... took :");

	//display_results(G,P,pp,F,A0,true);
	//print_free_fibers(G,pp,F,A0,true);

	////// Assignment global --------------------------------------------
	Assignment A(G,F);
	// Naive assignment
	init_time_at(time,"# Begin naive assignment",t);
	assign_fibers(G,P,pp,F,A);
	print_time(time,"# ... took");
	
	print_table("Petal",A.infos_petal(1000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(2000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(3000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(4000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(5000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(6000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(7000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(8000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(9000,5,G,P,pp,F));
	print_table("Petal",A.infos_petal(10000,5,G,P,pp,F));

	init_time_at(time,"# Begin display results and free fibers",t);
	display_results(G,P,pp,F,A,false);
	print_free_fibers(G,pp,F,A,false);
	print_time(time,"# ... took");

	// Improve --------------------------------------------------------
	
	init_time_at(time,"# Begin improve2",t);
	improve2("SF",G,P,pp,F,A);
	print_time(time,"# ... took :");
	display_results(G,P,pp,F,A,false);
	print_free_fibers(G,pp,F,A,false);

	init_time_at(time,"# Begin improve",t);
	improve(G,P,pp,F,A);
	print_time(time,"# ... took :");
	display_results(G,P,pp,F,A,false);
	print_free_fibers(G,pp,F,A,false);



	//init_time_at(time,"# Begin redistribute galaxies",t);
	//redistribute_g(G,P,pp,F,A);
	//print_time(time,"# ... took :"); 
	//display_results(G,P,pp,F,A,false);
	//print_free_fibers(G,pp,F,A,false);

	init_time_at(time,"# Begin redistribute SF",t);
	redistribute_g_by_kind("SF",G,P,pp,F,A);
	print_time(time,"# ... took :"); 
	display_results(G,P,pp,F,A,false);
	print_free_fibers(G,pp,F,A,false);

	//init_time_at(time,"# Begin improve",t);
	//improve(G,P,pp,F,A);
	//print_time(time,"# ... took :");
	//display_results(G,P,pp,F,A,false);
	//print_free_fibers(G,pp,F,A,false);

	//init_time_at(time,"# Begin redistribute TF",t);
	//redistribute_tf(G,P,pp,F,A);
	//print_time(time,"# ... took :"); 
	//display_results(G,P,pp,F,A,false);
	//print_free_fibers(G,pp,F,A,false);

	//init_time_at(time,"# Begin improve",t);
	//improve(G,P,pp,F,A);
	//print_time(time,"# ... took :");
	//display_results(G,P,pp,F,A,false);
	//print_free_fibers(G,pp,F,A,false);


	//plot_freefibers("free_fiber_plot.txt",P,A);   
	//time_line(G,F,A);

	print_time(t,"# Finished !... in");
	return(0);
}
