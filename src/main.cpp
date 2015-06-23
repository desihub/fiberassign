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
	Time t, time; // t for global, time for local
	init_time(t);
	Feat F;

	// Read parameters file
	F.readInputFile(argv[1]);
	printFile(argv[1]);

	// Read galaxies
	Gals G;
	G = read_galaxies(F);
	F.Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(F.Ngal).c_str(),F.galFile.c_str());

	// Read fiber center positions and compute related things
	PP pp;
	pp.read_fiber_positions(F); 
	F.Nfiber = pp.fp.size()/2; 
	F.Npetal = max(pp.spectrom)+1; // spectrometer has to be identified from 0 to F.Npetal-1
	F.Nfbp = (int) (F.Nfiber/F.Npetal);
	pp.get_neighbors(F); pp.compute_fibsofsp(F);

	// Read plates
	Plates P = read_plate_centers(F); 
	F.Nplate = P.size();
	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(F.Nplate).c_str(),F.tileFile.c_str(),F.Nfiber,F.fibFile.c_str());

	// Computes geometries of cb and fh
	F.cb = create_cb();  // cb=central body
	F.fh = create_fh();  // fh=fiber holder

	//// Collect available galaxies <-> tilefibers --------------------
	// HTM Tree of galaxies
	const double MinTreeSize = 0.01; // 
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

	// Make a plan ----------------------------------------------------
	Table unused; Table unused_time; int interv = 5; // For ploting histograms

	new_assign_fibers(G,P,pp,F,A); // Plans whole survey

	// Smooth out distribution of free fibers
	for (int i=0; i<9; i++) {
		print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
		unused.push_back(histogram(A.unused_fbp(pp,F),1));
		redistribute_tf(G,P,pp,F,A);
	}
	// Still not updated, so all QSO targets have multiple observations etc
	for (int i=0; i<5; i++) {
		print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
		unused.push_back(histogram(A.unused_fbp(pp,F),1));
		improve(G,P,pp,F,A);
		print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
		unused.push_back(histogram(A.unused_fbp(pp,F),1));
		redistribute_tf(G,P,pp,F,A);
		print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
		unused.push_back(histogram(A.unused_fbp(pp,F),1));
		redistribute_tf(G,P,pp,F,A);
	}
	for (int i=0; i<4; i++) {
		print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
		unused.push_back(histogram(A.unused_fbp(pp,F),1));
		redistribute_tf(G,P,pp,F,A);
	}

	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=0; jj<F.Nplate; jj++) {
		int j = A.next_plate;
		printf(" - Plate %d :",j);
		assign_sf_ss(j,G,P,pp,F,A);
		assign_left(j,G,P,pp,F,A); //if unused fibers, assign them
		if (j%500==0) pyplotTile(j,"doc/figs",G,P,pp,F,A);//beautiful picture of positioners, galaxies
		// <-- here is the real observation time
		printf(" %s not as - ",format(5,f(A.unused_f(j,F))).c_str()); fl();
		if (0<=j-F.Analysis) update_plan_from_one_obs(G,P,pp,F,A,F.Nplate-1); else printf("\n");
		// Update corrects all future occurrences of wrong QSOs etc and tries to observe something else
		A.next_plate++;

		if (j==100 || j==200 || j==700 || j==1500 || j==3000 || j==5000 || j==7000) {  //redistribute and improve on various occasions
			redistribute_tf(G,P,pp,F,A);
			redistribute_tf(G,P,pp,F,A);
			improve(G,P,pp,F,A);
			redistribute_tf(G,P,pp,F,A);
		}
	}
	print_time(time,"# ... took :");

	// Results -------------------------------------------------------
	if (F.Output) for (int j=0; j<F.Nplate; j++) write_FAtile_ascii(j,F.outDir,G,P,pp,F,A); // Write output
	//overlappingTiles("overlaps.txt",F,A); // Write some overlapping tiles (for S.Bailey)
	print_hist("Unused fibers",interv,histogram(A.unused_fbp(pp,F),interv),false);
	unused.push_back(histogram(A.unused_fbp(pp,F),1));
	print_mult_table_latex("Histograms of unused fibers","doc/figs/unusedfibs.dat",unused,1);

	init_time_at(time,"# Display results",t);
	display_results("doc/figs/",G,P,pp,F,A,true);
	print_time(time,"# ... took :");
	//if (F.Verif) A.verif(P,G,pp,F); // Verification that the assignment is sane
	print_time(t,"# Finished !... in");
	return(0);
}
