#include    <cstdlib>
#include    <cmath>
#include    <fstream>
#include    <sstream>
#include    <iostream>
#include    <iomanip>
#include    <string>
#include    <vector>
#include    <algorithm>
#include    <exception>
#include    <sys/time.h>
#include    "modules/htmTree.h"
#include    "modules/kdTree.h"
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "collision.h"
#include        "global.h"
#include     <map>
#include  <stdexcept>
//reduce redistributes, updates  07/02/15 rnc
int main(int argc, char **argv) {
    // argv[1] is the features file
    //// Initializations ---------------------------------------------
    //
    //  M is collection of targets, sky fibers, and standard stars
    //  each target has a priority provided by the mtl file
    //  all targets with the same priority are collected into a class

    check_args(argc);
    Time t, time; // t for global, time for local
    init_time(t);
    Feat F;
    MTL M;

    // Read parameters file //
    F.readInputFile(argv[1]);
    printFile(argv[1]);
    // Read input files for standards, skys and targets.
    // Try to read SS and SF before targets to avoid wasting time if these
    // smaller files can't be read.
    init_time_at(time,"# Read target, SS, SF files",t);
    MTL SStars = read_MTLfile(F.SStarsfile,F,1,0);
    MTL SkyF   = read_MTLfile(F.SkyFfile,  F,0,1);
    MTL Targ   = read_MTLfile(F.Targfile,  F,0,0);

    print_time(time,"# ...read targets  took :");
    //combine the three input files
    M=Targ;
    printf(" Target size %d \n",M.size());
    std::cout.flush();
    //need to be able to match immutable target id to position in list
    //check for duplicates on mtl only to allow duplication with SS
    std::map<long long,int> invert_target;
    std::map<long long,int>::iterator targetid_to_idx;
    std::pair<std::map<long long,int>::iterator,bool> ret;
    for(unsigned i=0;i<M.size();++i)
      {
	ret = invert_target.insert(std::make_pair(M[i].id,i));
	//check for duplicates (std::map.insert only created keys, fails on duplicate keys)
	if(ret.second == false ){
	  std::ostringstream o;
	  o<<"Duplicate targetid "<<M[i].id<<" in MTL";
	  throw std::logic_error(o.str().c_str());
	}
      }
    M.insert(M.end(),SStars.begin(),SStars.end());
    printf(" Standard Star size %d \n",M.size());
    M.insert(M.end(),SkyF.begin(),SkyF.end());
    printf(" Sky Fiber size %d \n",M.size());


    init_time_at(time,"# map position in target list to immutable targetid",t);

    init_time_at(time,"# assign priority classes",t);
    F.Ngal = M.size();
    assign_priority_class(M);
    std::vector <int> count_class(M.priority_list.size(),0);
    for(int i=0;i<M.size();++i){
        if(!M[i].SS&&!M[i].SF){
            count_class[M[i].priority_class]+=1;
        }
    }
    for(int i=0;i<M.priority_list.size();++i){
        printf("  class  %d  number  %d\n",i,count_class[i]);
    }
    print_time(time,"# ...priority list took :");

    init_time_at(time,"# Start positioners",t);
    // fiber positioners
    
    F.Npetal = 10;//spectrometers run 0 to 9 unless pacman
    FP pp =read_fiber_positions(F); 
    read_fiber_status(pp, F);
    //order the fibers by their fiber number (fib_num) not simply order in list
    //need to fix spectrom (List) and fp
    
    F.Nfiber = pp.size(); //each fiber has two co-ordinates so divide by two
    F.Nfbp = F.Nfiber/F.Npetal;// fibers per petal = 500

    print_time(time,"# ..posiioners  took :");
    //
    init_time_at(time,"# Start plates",t);
    //P is original list of plates

    Plates P = read_plate_centers(F);

    F.Nplate=P.size();

    printf("# Read %s plate centers from %s and %d fibers from %s\n",f(F.Nplate).c_str(),F.tileFile.c_str(),F.Nfiber,F.fibFile.c_str());
    //    
    print_time(time,"# ..plates   took :");
    // Computes geometries of cb and fh: pieces of positioner - used to determine possible collisions
    F.cb = create_cb(); // cb=central body
    F.fh = create_fh(); // fh=fiber holder

    //// Collect available galaxies <-> tilefibers --------------------
    // HTM Tree of galaxies
    const double MinTreeSize = 0.01;
    init_time_at(time,"# Start building HTM tree",t);
    htmTree<struct target> T(M,MinTreeSize);
    print_time(time,"# Doing kd-tree... took :");//T.stats();
    init_time_at(time,"# collect galaxies at ",t);

    // For plates/fibers, collect available galaxies; done in parallel  P[plate j].av_gal[k]=[g1,g2,..]
    collect_galaxies_for_all(M,T,P,pp,F);

    print_time(time,"# ... took :");//T.stats();
    init_time_at(time,"# collect available tile-fibers at",t);

    // For each galaxy, computes available tilefibers  G[i].av_tfs = [(j1,k1),(j2,k2),..]
    collect_available_tilefibers(M,P,F);

    //// Assignment |||||||||||||||||||||||||||||||||||||||||||||||||||
    printf(" Nplate %d  Ngal %d   Nfiber %d \n", F.Nplate, F.Ngal, F.Nfiber);
    Assignment A(M,F);
    // Make a plan ----------------------------------------------------
    print_time(t,"# Start assignment at : ");
    simple_assign(M,P,pp,F,A);

    //check to see if there are tiles with no galaxies
    //need to keep mapping of old tile list to new tile list and inverse map
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

    // Smooth out distribution of free fibers, and increase the number of assignments
    // probably should not hard wire the limits i<1, i<3 in redistribute and improve
       
    for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A,0);// more iterations will improve performance slightly
    for (int i=0; i<1; i++) {
        improve(M,P,pp,F,A,0);
        redistribute_tf(M,P,pp,F,A,0);
	}

    init_time_at(time,"# assign SS and SF ",t);

    //try assigning SF and SS before real time assignment
    for (int jused=0;jused<F.NUsedplate;++jused){
        int j=A.suborder[jused];
        assign_sf_ss(j,M,P,pp,F,A); // Assign SS and SF for each tile
        assign_unused(j,M,P,pp,F,A);
    }
    
    //fill all unused fibers with sky fibers
    for (int jused=0;jused<F.NUsedplate;++jused){
        int j=A.suborder[jused];
        fill_unused_with_sf(j,M,P,pp,F,A);
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

    }
    init_time_at(time,"# count SS and SF ",t);
    printf(" Totals SS   %4d    SF   %4d",total_used_SS,total_used_SF);

    for (int pr=0;pr<M.priority_list.size();++pr){
        printf(" class %2d   %5d",pr,total_used_by_class[pr]);
    }
    printf("\n");

    init_time_at(time,"# print fits files ",t);
    if (F.PrintFits) for (int jused=0; jused<F.NUsedplate; jused++){
        int j=A.suborder[jused];
        fa_write(j,F.outDir,M,P,pp,F,A); // Write output
    }
  
  print_time(t,"# Finished !... in");
  
  return(0);    
}

