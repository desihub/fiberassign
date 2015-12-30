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
    
    init_time_at(time,"# read target, SS, SF files",t);
    MTL Targ=read_MTLfile(F.Targfile,F,0,0);
    MTL SStars=read_MTLfile(F.SStarsfile,F,1,0);
    MTL SkyF=read_MTLfile(F.SkyFfile,F,0,1);
    print_time(time,"# ... took :");
    //combine the three input files
    M=Targ;
    printf(" Target size %d \n",M.size());
    M.insert(M.end(),SStars.begin(),SStars.end());
    printf(" Standard Star size %d \n",M.size());
    M.insert(M.end(),SkyF.begin(),SkyF.end());
    printf(" Sky Fiber size %d \n",M.size());
    
    F.Ngal = M.size();
    assign_priority_class(M);
    std::vector <int> count_class(M.priority_list.size(),0);
    for(int i;i<M.size();++i){
        if(!M[i].SS&&!M[i].SF){
            count_class[M[i].priority_class]+=1;
        }
    }
    for(int i;i<M.priority_list.size();++i){
        printf("  class  %d  number  %d\n",i,count_class[i]);
    }
    print_time(time,"# ... took :");
    
    // fiber positioners
	PP pp;
	pp.read_fiber_positions(F); 
	F.Nfiber = pp.fp.size()/2; 
	F.Npetal = max(pp.spectrom)+1;
	F.Nfbp = (int) (F.Nfiber/F.Npetal);// fibers per petal = 500
	pp.get_neighbors(F);
    pp.compute_fibsofsp(F);
    
    //P is original list of plates
	Plates P = read_plate_centers(F);
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
    
	// For plates/fibers, collect available galaxies; done in parallel  P[plate j].av_gal[k]=[g1,g2,..]
	collect_galaxies_for_all(M,T,P,pp,F);
    print_time(time,"# ... took :");//T.stats();
    init_time_at(time,"# collect available tile-fibers at",t);
    
	// For each galaxy, computes available tilefibers  G[i].av_tfs = [(j1,k1),(j2,k2),..]
	collect_available_tilefibers(M,P,F);

	//results_on_inputs("doc/figs/",G,P,F,true);

	//// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
    printf(" Nplate %d  Ngal %d   Nfiber %d \n", F.Nplate, F.Ngal, F.Nfiber);
	Assignment A(M,F);
    // Make a plan ----------------------------------------------------
	print_time(t,"# Start assignment at : ");
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
                A.suborder.push_back(j);//suborder[jused] is jused-th used plate
                not_done=false;
                A.inv_order[j]=inv_count;//inv_order[j] is -1 unless used
                inv_count++;
            }
        }
    }
    F.NUsedplate=A.suborder.size();
    printf(" Plates actually used %d \n",F.NUsedplate);
    for(int i=0;i<F.NUsedplate;i++)printf(" jused  %d  j  %d\n",i,A.suborder[i]);
    

    
    print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false); // Hist of unused fibs
    
    // Smooth out distribution of free fibers, and increase the number of assignments
    
    for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A,0);// more iterations will improve performance slightly
    for (int i=0; i<3; i++) {
        improve(M,P,pp,F,A,0);
        redistribute_tf(M,P,pp,F,A,0);
    }
    print_hist("Unused fibers",5,histogram(A.unused_fbp(pp,F),5),false);
    //try assigning SF and SS before real time assignment
    for (int jused=0;jused<F.NUsedplate;++jused){
        
        int j=A.suborder[jused];
        assign_sf_ss(j,M,P,pp,F,A); // Assign SS and SF for each tile
        assign_unused(j,M,P,pp,F,A);
    }

    

	// Results -------------------------------------------------------*/
    std::vector <int> total_used_by_class(M.priority_list.size(),0);
    int total_used_SS=0;
    int total_used_SF=0;
    for (int jused=0;jused<F.NUsedplate;++jused){
        std::vector <int> used_by_class(M.priority_list.size(),0);
        int used_SS=0;
        int used_SF=0;
        int j=A.suborder[jused];
        for(int k=0;k<F.Nfiber;++k){
            int g=A.TF[j][k];
            if(g!=-1){
                if(M[g].SS){
                    total_used_SS++;
                    used_SS++;
                    }
                    else if(M[g].SF){
                        used_SF++;
                        total_used_SF++;
                    }
                    else{
                        used_by_class[M[g].priority_class]++;
                        total_used_by_class[M[g].priority_class]++;
                    }
            }
        }
        printf(" plate jused %5d j %5d  SS   %4d    SF   %4d",jused,j,used_SS,used_SF);
        for (int pr=0;pr<M.priority_list.size();++pr){
            printf(" class %2d   %5d",pr,used_by_class[pr]);
        }
        printf("\n");
    printf(" plate jused %5d j %5d  SS   %4d    SF   %4d",jused,j,total_used_SS,total_used_SF);
    for (int pr=0;pr<M.priority_list.size();++pr){
        printf(" class %2d   %5d",pr,total_used_by_class[pr]);
    }
    printf("\n");
        
        
        

    
                        
                                        
                             
                
                                            
                                    
    
    
    if (F.PrintAscii) for (int j=0; j<F.Nplate; j++){
            write_FAtile_ascii(j,F.outDir,M,P,pp,F,A);
        }

    if (F.PrintFits) for (int j=0; j<F.Nplate; j++){
        fa_write(j,F.outDir,M,P,pp,F,A); // Write output
    }
    /*
	display_results("doc/figs/",G,M,P,pp,F,A,true);
	if (F.Verif) A.verif(P,M,pp,F); // Verification that the assignment is sane
     */

    print_time(t,"# Finished !... in");
    
    return(0);    
}
