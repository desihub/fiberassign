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
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "collision.h"
#include        "global.h"

int main(int argc, char **argv) {
	//// Initializations ---------------------------------------------
	srand48(1234); // Make sure we have reproducability
	check_args(argc);
	Time t,time; // t for global, time for local
	init_time(t);
	Feat F;

	// Read parameters file
	F.readInputFile(argv[1]);
	printFile(argv[1]);

	// Read galaxies
	Gals G;
	G = read_galaxies(F); // ! Reads all galaxies iff arg2=1
	F.Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(F.Ngal).c_str(),F.galFile.c_str());

	// Read fiber center positions and compute related things
	PP pp;
	pp.read_fiber_positions(F); 
	F.Nfiber = pp.fp.size()/2; 
	F.Npetal = max(pp.spectrom)+1; // spectrom has to be identified from 0 to F.Npetal-1
	F.Nfbp = (int) (F.Nfiber/F.Npetal);
	pp.get_neighbors(F); pp.compute_fibsofsp(F);

	// Read plates
	Plates P = read_plate_centers(F); // ! Reads all galaxies iff arg2=1
	F.Nplate = P.size();
	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(F.Nplate).c_str(),F.tileFile.c_str(),F.Nfiber,F.fibFile.c_str());

	// Computes geometry of cb and fh
	F.cb = create_cb();
	F.fh = create_fh();

	//// Collect available galaxies <-> tilefibers --------------------
	// HTM Tree of galaxies
	const double MinTreeSize = 0.01; // <--- ?
	init_time_at(time,"# Start building HTM tree",t);
	htmTree<struct galaxy> T(G,MinTreeSize);
	print_time(time,"# ... took :");
	//T.stats();
	
	// For plates/fibers, collect available galaxies
	collect_galaxies_for_all(G,T,P,pp,F);

	// For each galaxy, computes available tilefibers
	collect_available_tilefibers(G,P,F);

	//results_on_inputs("doc/figs/",G,P,F,true);

	//// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
	Assignment A(G,F);
	print_time(t,"# Start assignment at : ");
	int limit = 2000;
	// Make a plan ----------------------------------------------------
	//new_assign_fibers(G,P,pp,F,A,2000);
	simple_assign(G,P,pp,F,A,limit);
	improve_from_kind(G,P,pp,F,A,"SF",limit);

	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=0; jj<limit; jj++) {
		int j = A.next_plate;
		printf(" - Plate %d : ",j);
		//improve_from_kind(G,P,pp,F,A,"SF",1);
		//improve_from_kind(G,P,pp,F,A,"SS",1);
		if (j%500==0) pyplotTile(j,"doc/figs",G,P,pp,F,A);
		// <-- here is the real observation time
		printf(" %s not as - ",format(5,f(F.Nfiber-A.na(F,j,1))).c_str());
		if (0<=j-F.Analysis) update_plan_from_one_obs(G,P,pp,F,A,limit-1);
		else printf("\n");
		A.next_plate++;
	}
	print_time(time,"# ... took :");

	// Make a plan ----------------------------------------------------
	//new_assign_fibers(G,P,pp,F,A);
	simple_assign(G,P,pp,F,A);
	improve(G,P,pp,F,A);
	improve_from_kind(G,P,pp,F,A,"SF");
	improve_from_kind(G,P,pp,F,A,"SS");

	redistribute_tf(G,P,pp,F,A);
	improve(G,P,pp,F,A);
	redistribute_tf(G,P,pp,F,A);
	improve(G,P,pp,F,A);

	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=limit; jj<F.Nplate; jj++) {
		int j = A.next_plate;
		printf(" - Plate %d :",j);
		improve_from_kind(G,P,pp,F,A,"SF",1);
		improve_from_kind(G,P,pp,F,A,"SS",1);
		if (j%500==0) pyplotTile(j,"doc/figs",G,P,pp,F,A);
		// <-- here is the real observation time
		printf(" %s not as - ",format(5,f(F.Nfiber-A.na(F,j,1))).c_str());
		update_plan_from_one_obs(G,P,pp,F,A,F.Nplate-1);
		A.next_plate++;
	}
	print_time(time,"# ... took :");

	// Results -------------------------------------------------------
	//for (int j=0; j<F.Nplate; j++) write_FAtile(j,F.outDir,G,P,pp,F,A); // Write output
	//overlappingTiles("overlaps.txt",F,A); // Write some overlapping tiles (for S.Bailey)
	init_time_at(time,"# Display results",t);
	display_results("doc/figs/",G,P,pp,F,A,true);
	print_time(time,"# ... took :");
	if (F.Verif) A.verif(P,G,pp,F); // Verification that the assignment is sane
	print_time(t,"# Finished !... in");
	return(0);
}
