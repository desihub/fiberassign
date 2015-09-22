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
//reduce redistributes, updates  07/02/15 rnc
int main(int argc, char **argv) {
	//// Initializations ---------------------------------------------
	srand48(1234); // Make sure we have reproducability
	check_args(argc);
	Time t, time; // t for global, time for local
	init_time(t);
	Feat F;

	// Read parameters file //
	F.readInputFile(argv[1]);
	printFile(argv[1]);
	// Read galaxies
	Gals G;
    if(F.Ascii){
        G=read_galaxies_ascii(F);}
    else{
        G = read_galaxies(F);
    }
	F.Ngal = G.size();
	printf("# Read %s galaxies from %s \n",f(F.Ngal).c_str(),F.galFile.c_str());
    std::vector<int> count;
    count=count_galaxies(G);
    for(int i=0;i<8;i++){printf ("   %d  \n",count[i]);}
    // make MTL
    MTL M=make_MTL(G,F);
    for(int i=0;i<M.priority_list.size();++i){
        printf(" priority   %f",M.priority_list[i]);
    }
    printf(" \n");
    
    std::vector <int> count_class(M.priority_list.size(),0);
    for(int i;i<M.size();++i){
        count_class[M[i].priority_class]+=1;
    }
    for(int i;i<M.priority_list.size();++i){
        printf("  class  %d  number  %d\n",i,count_class[i]);
    }
    
    printf(" number of MTL galaxies  %d\n",M.size());
    printf("Read fiber center positions and compute related things\n");
	PP pp;
	pp.read_fiber_positions(F); 
	F.Nfiber = pp.fp.size()/2; 
	F.Npetal = max(pp.spectrom)+1;
    printf("spectrometer has to be identified from 0 to F.Npetal-1\n");
	F.Nfbp = (int) (F.Nfiber/F.Npetal);// fibers per petal = 500
	pp.get_neighbors(F); pp.compute_fibsofsp(F);
    printf("get neighbors of each fiber;\n");
                                                //for each spectrometer, get list of fibers

    printf("Read plates in order they are to be observed\n ");
    
	Plates P_original = read_plate_centers(F);
    F.Nplate=P_original.size();
    printf("This takes the place of Opsim or NextFieldSelector; will be replaced by appropriate code\n");
    Plates P;
    P.resize(F.Nplate);
    List permut = random_permut(F.Nplate);
    for (int jj=0; jj<F.Nplate; jj++){
        P[jj]=P_original[permut[jj]];
    }

	printf("# Read %s plate centers from %s and %d fibers from %s\n",f(F.Nplate).c_str(),F.tileFile.c_str(),F.Nfiber,F.fibFile.c_str());

	// Computes geometries of cb and fh: pieces of positioner - used to determine possible collisions
	F.cb = create_cb(); // cb=central body
	F.fh = create_fh(); // fh=fiber holder

	//// Collect available galaxies <-> tilefibers --------------------
	// HTM Tree of galaxies
	const double MinTreeSize = 0.01;
	init_time_at(time,"# Start building HTM tree",t);
	htmTree<struct target> T(M,MinTreeSize);
	print_time(time,"# ... took :");//T.stats();
	
	// For plates/fibers, collect available galaxies; done in parallel  P[plate j].av_gal[k]=[g1,g2,..]
	collect_galaxies_for_all(M,T,P,pp,F);
	
	// For each galaxy, computes available tilefibers  G[i].av_tfs = [(j1,k1),(j2,k2),..]
	collect_available_tilefibers(M,P,F);

	//results_on_inputs("doc/figs/",G,P,F,true);

	//// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
    printf("before assignment\n");
	Assignment A(M,F);
	print_time(t,"# Start assignment at : ");

	// Make a plan ----------------------------------------------------
	//new_assign_fibers(G,P,pp,F,A); // Plans whole survey without sky fibers, standard stars
                                   // assumes maximum number of observations needed for QSOs, LRGs
    printf(" Nplate %d  Ngal %d   Nfiber %d \n", F.Nplate, F.Ngal, F.Nfiber);
 
    simple_assign(M,P,pp,F,A);
    
    
	print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false); // Hist of unused fibs
                                    // Want to have even distribution of unused fibers
                                    // so we can put in sky fibers and standard stars

	// Smooth out distribution of free fibers, and increase the number of assignments
	for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A);// more iterations will improve performance slightly
	for (int i=0; i<1; i++) {                           // more iterations will improve performance slightly
		improve(M,P,pp,F,A);
		redistribute_tf(M,P,pp,F,A);
		redistribute_tf(M,P,pp,F,A);
	}
	for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A);

	print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false);

	// Still not updated, so all QSO targets have multiple observations etc
	// Apply and update the plan --------------------------------------
	init_time_at(time,"# Begin real time assignment",t);
	for (int jj=0; jj<F.Nplate; jj++) {
		int j = A.next_plate;
		//printf(" - Plate %d :",j);
		//assign_sf_ss(j,M,P,pp,F,A); // Assign SS and SF just before an observation
		assign_unused(j,M,P,pp,F,A);
		//if (j%2000==0) pyplotTile(j,"doc/figs",G,P,pp,F,A); // Picture of positioners, galaxies
		
		//printf(" %s not as - ",format(5,f(A.unused_f(j,F))).c_str()); fl();
		// Update corrects all future occurrences of wrong QSOs etc and tries to observe something else
		if (0<=j-F.Analysis) update_plan_from_one_obs(G,M,P,pp,F,A,F.Nplate-1); else printf("\n");
		A.next_plate++;
		// Redistribute and improve on various occasions  add more times if desired

		if ( j==2000 || j==4000) {
			redistribute_tf(M,P,pp,F,A);
			redistribute_tf(M,P,pp,F,A);
			improve(M,P,pp,F,A);
			redistribute_tf(M,P,pp,F,A);
		}
	}
	print_time(time,"# ... took :");

	// Results -------------------------------------------------------
	if (F.Output) for (int j=0; j<F.Nplate; j++) write_FAtile_ascii(j,F.outDir,M,P,pp,F,A); // Write output
	display_results("doc/figs/",M,P,pp,F,A,true);
	if (F.Verif) A.verif(P,M,pp,F); // Verification that the assignment is sane


	print_time(t,"# Finished !... in");
 
	return(0);
  
}
