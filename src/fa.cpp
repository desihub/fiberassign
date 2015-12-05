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
    MTL M;
   

	// Read parameters file //
	F.readInputFile(argv[1]);
	printFile(argv[1]);
	// Read Secretfile
    // Secret contains the identity of each target: QSO-Ly-a, QSO-tracers, LRG, ELG, fake QSO, fake LRG, SS, SF
    Gals Secret;
    printf("before read secretfile \n");
    init_time_at(time,"# read Secret file",t);

    Secret=read_Secretfile(F.Secretfile,F);
    printf("# Read %d galaxies from %s \n",Secret.size(),F.Secretfile.c_str());
	print_time(time,"# ... took :");
    std::vector<int> count;
    count=count_galaxies(Secret);
    printf(" Number of galaxies by type, QSO-Ly-a, QSO-tracers, LRG, ELG, fake QSO, fake LRG, SS, SF\n");
    for(int i=0;i<8;i++){printf (" type %d number  %d  \n",i, count[i]);}
    //read the three input files
    init_time_at(time,"# read target, SS, SF files",t);
    MTL Targ=read_MTLfile(F.Targfile,F,0,0);
    MTL SStars=read_MTLfile(F.SStarsfile,F,1,0);
    MTL SkyF=read_MTLfile(F.SkyFfile,F,0,1);
    print_time(time,"# ... took :");
    //combine the three input files

    M=Targ;
   
    printf(" M size %d \n",M.size());
    M.insert(M.end(),SStars.begin(),SStars.end());
    printf(" M size %d \n",M.size());
    M.insert(M.end(),SkyF.begin(),SkyF.end());
    printf(" M size %d \n",M.size());
    F.Ngal=M.size();
    assign_priority_class(M);
    
    //establish priority classes
    init_time_at(time,"# establish priority clasess",t);
    std::vector <int> count_class(M.priority_list.size(),0);
    for(int i;i<M.size();++i){
        if(!M[i].SS&&!M[i].SF){
        count_class[M[i].priority_class]+=1;
        }
    }
	print_time(time,"# ... took :");
   
    for(int i;i<M.priority_list.size();++i){
        printf("  class  %d  number  %d\n",i,count_class[i]);
    }
    
    PP pp;
	pp.read_fiber_positions(F); 
	F.Nfiber = pp.fp.size()/2; 
	F.Npetal = max(pp.spectrom)+1;
    F.Nfbp = (int) (F.Nfiber/F.Npetal);// fibers per petal = 500
	pp.get_neighbors(F); pp.compute_fibsofsp(F);
    
    //P is original list of plates
	Plates P = read_plate_centers(F);

    printf(" future plates %d\n",P.size());
    F.Nplate=P.size();
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
    init_time_at(time,"# collect galaxies at ",t);
	
	// For plates/fibers, collect available galaxies; done in parallel
    collect_galaxies_for_all(M,T,P,pp,F);
    print_time(time,"# ... took :");//T.stats();
    init_time_at(time,"# collect available tile-fibers at",t);
	// For each galaxy, computes available tilefibers  G[i].av_tfs = [(j1,k1),(j2,k2),..]
	collect_available_tilefibers(M,P,F);
	
	//results_on_inputs("doc/figs/",G,P,F,true);

	//// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
    printf(" Nplate %d  Ngal %d   Nfiber %d \n", F.Nplate, F.Ngal, F.Nfiber);
    Assignment A(M,F);
    
	print_time(t,"# Start assignment at : ");

	// Make a plan ----------------------------------------------------
    // Plans whole survey without sky fibers, standard stars
    // assumes maximum number of observations needed for QSOs, LRGs

    simple_assign(M,P,pp,F,A);

    //check to see if there are tiles with no galaxies
    //need to keep mapping of old tile list to new tile list
    //and inverse map
    A.inv_order=initList(F.Nplate,-1);
    int inv_count=0;
    for (int j=0;j<F.Nplate ;++j){

        bool not_done=true;
        for(int k=0;k<F.Nfiber && not_done;++k){
            if(A.TF[j][k]!=-1){
                A.suborder.push_back(j);
                not_done=false;
                inv_count++;
                A.inv_order[j]=inv_count;
                
            }
        }
    }
    F.NUsedplate=A.suborder.size();
    printf(" Plates after screening %d \n",F.NUsedplate);
    
    //if(F.diagnose)diagnostic(M,G,F,A);

    print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false); // Hist of unused fibs
    
	// Smooth out distribution of free fibers, and increase the number of assignments
    
	for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A,0);// more iterations will improve performance slightly
	for (int i=0; i<1; i++) {
        improve(M,P,pp,F,A,0);
		redistribute_tf(M,P,pp,F,A,0);
	}
	print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false);
    //try assigning SF and SS before real time assignment
    for (int j=0;j<F.NUsedplate;++j){

        int js=A.suborder[j];
        assign_sf_ss(js,M,P,pp,F,A); // Assign SS and SF for each tile
        assign_unused(js,M,P,pp,F,A);
    }
    if(F.diagnose)diagnostic(M,Secret,F,A);
    init_time_at(time,"# Begin real time assignment",t);

	//Execute plan, updating targets at intervals
    for(int i=0;i<F.pass_intervals.size();i++)printf(" i=%d interval %d \n",i,F.pass_intervals[i]);
    std::vector <int> update_intervals=F.pass_intervals;
    update_intervals.push_back(F.NUsedplate);//to end intervals at last plate
    for(int i=0;i<update_intervals.size()-1;++i){//go plate by used plate
        int starter=update_intervals[i];

        //display_results("doc/figs/",G,P,pp,F,A,true);
        //plan whole survey from this point out
        /*
        for (int jj=starter; jj<F.NUsedplate; jj++) {
            int js = A.suborder[jj];
            assign_sf_ss(js,M,P,pp,F,A); // Assign SS and SF
            assign_unused(js,M,P,pp,F,A);
        }
        */
        //update target information for interval i

        for (int jj=starter; jj<update_intervals[i+1]; jj++) {
            if (0<=jj-F.Analysis) update_plan_from_one_obs(jj,Secret,M,P,pp,F,A); else printf("\n no update\n");
            // Update corrects all future occurrences of wrong QSOs etc and tries to observe something else
           
        }
    
    
 

        redistribute_tf(M,P,pp,F,A,starter);
        improve(M,P,pp,F,A,starter);
        redistribute_tf(M,P,pp,F,A,starter);
        //}

        if(F.diagnose)diagnostic(M,Secret,F,A);
    }
    // check on SS and SF
/*
    for(int j=0;j<F.NUsedplate;++j){
        int js=A.suborder[j];
        printf("\n js = %d\n",js);
        for (int p=0;p<F.Npetal;++p){
            int count_SS=0;
            int count_SF=0;
            for (int k=0;k<F.Nfbp;++k){
                int kk=pp.fibers_of_sp[p][k];
                int g=A.TF[js][kk];
                if(g!=-1 && M[g].SS)count_SS++;
                if(g!=-1 && M[g].SF)count_SF++;
                
            }
            printf("  %d  %d   ",count_SS,count_SF);
        }
        printf("\n");
    }
     */
    

 
	// Results -------------------------------------------------------
    if (F.PrintAscii) for (int j=0; j<F.NUsedplate; j++){
        write_FAtile_ascii(A.suborder[j],F.outDir,M,P,pp,F,A);
    }
    
    if (F.PrintFits) for (int j=0; j<F.NUsedplate; j++){
        fa_write(A.suborder[j],F.outDir,M,P,pp,F,A); // Write output
    }
    

	display_results("doc/figs/",Secret,M,P,pp,F,A,true);
	if (F.Verif) A.verif(P,M,pp,F); // Verification that the assignment is sane


	print_time(t,"# Finished !... in");

	return(0);
  
}
