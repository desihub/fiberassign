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
#include        "structs.h"
#include        "misc.h"

int MaxPlate;
int Ngal;

int main(int argc, char **argv) {
	//// Initialization ---------------------------------------------------
	check_args(argc);
	srand48(1234); // Make sure we have reproducability
	Time t,time; // t for global, time for local
	init_time(t);
	// Number of observations desired for a particular galaxy type
	List goal;
	goal.assign(7,1); goal[1]=5;  goal[3]=2;
	// Galaxies
	Gals G;
	G = read_galaxies(argv[1],1); // ! Reads all galaxies iff arg2=1
	Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(Ngal),argv[1]);
	// Read fiber center positions and compute neighbors (pp for plate params)
	PP pp;
	pp.read_fiber_positions(argv[3]); 
	pp.get_neighbors(); 
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
	// For each galaxy, computes available tilefibers
	init_time_at(time,"# Begin computing available tilefibers",t);
	collect_available_tilefibers(G,P);
	print_time(time,"# ... took :");

	//P[0].print_plate();
	// av_tfs[g].size ~ 7 +/- 3

	//// Assignment -----------------------------------------------------
	Assignment A;
	// Naive assignment
	init_time_at(time,"# Begin naive assignment",t);
	assign_fibers(G,P,pp,A);
	print_time(time,"# ... took");

	//Table Q = conflicts(G,P,pp,A);
	//display_results(G,goal,P,A);

	//// Improve --------------------------------------------------------
	init_time_at(time,"# Begin improve",t);
	improve(G,P,pp,A);
	print_time(time,"# ... took :");
	display_results(G,goal,P,A);
	print_free_fibers(pp,A);

	init_time_at(time,"# Begin redistribute",t);
	redistribute(G,P,pp,A);
	print_time(time,"# ... took :"); 
	print_free_fibers(pp,A);

	init_time_at(time,"# Begin redistribute",t);
	redistribute(G,P,pp,A);
	print_time(time,"# ... took :"); 
	print_free_fibers(pp,A);


	//plot_freefibers("free_fiber_plot.txt",P,A);   
	//time_line(G,goal,A);

	// write the results -- we should only need an unsigned char to
	// store the number of passes left.
	//
	// This section will output ra, dec, redshift and type of all targets
	// along with nobs; nobs == 0 means that the target has been observed
	// sufficient number of times and its redshift is very likely measured [L Samushia].
	// For now I'm outputing the file in ascii for convenience [L Samushia].
	// (erased for now)
	print_time(t,"# Finished !... in");
	return(0);
}
