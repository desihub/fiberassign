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
//reduce redistributes, updates  07/02/15 rnc
int main(int argc, char **argv) {
    //// Initializations ---------------------------------------------
    srand48(1234); // Make sure we have reproducability
    check_args(argc);
    Time t, time; // t for global, time for local
    init_time(t);
    Feat F;
    MTL M;
    bool diagnose=true;
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
    print_time(time,"# ... took :");
    //combine the three input files
    M=Targ;
    printf(" Target size %d \n",M.size());
    M.insert(M.end(),SStars.begin(),SStars.end());
    printf(" Standard Star size %d \n",M.size());
    M.insert(M.end(),SkyF.begin(),SkyF.end());
    printf(" Sky Fiber size %d \n",M.size());

    //need map between position in M and potentialtargetid
    std::map<long long,int> invert_target;
    std::map<long long,int>::iterator targetid_to_idx;
    std::pair<std::map<long long,int>::iterator,bool> ret;
    for(unsigned i=0;i<M.size();++i)
    {
        ret = invert_target.insert(std::make_pair(M[i].id,i));
        // Check for duplicates (std::map.insert only creates keys, fails on duplicate keys)
	/*
        if ( ret.second == false ) {
            std::ostringstream o;
            o << "Duplicate tileid " << M[i].id << " in tileFile!";
            throw std::logic_error(o.str().c_str());
        }
	*/
    }



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
    bool savetime=true; 
    Plates P = read_plate_centers(F);
    F.Nplate=P.size();
    printf("# Read %s plate centers from %s and %d fibers from %s\n",f(F.Nplate).c_str(),F.tileFile.c_str(),F.Nfiber,F.fibFile.c_str());    
    // Computes geometries of cb and fh: pieces of positioner - used to determine possible collisions
    F.cb = create_cb(); // cb=central body
    F.fh = create_fh(); // fh=fiber holder
    const char *filecheck="/project/projectdirs/desi/users/rncahn/large_mock_test/tmp/fiberassign/save_av_gals_00004.fits"; 
    //// Collect available galaxies <-> tilefibers --------------------
    // HTM Tree of galaxies

    std::ifstream infile(filecheck); 
    if(!infile.good()||!savetime){//av_gals has not been saved yet or not saving time *************

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


    } else { //  av_gals have been saved
          // each new epoch has shorter list of tiles
          size_t cfilesize = 512;
	  char filename[cfilesize];
          printf(" Will read saved files F.Nplate = %d \n",F.Nplate);
	  std::cout.flush();
	  for(int j=0;j<F.Nplate;j++){
	    Table av_gals,ss_av_gals,sf_av_gals;//here we have int not long long so Table is ok
	    //printf(" j %d to read back P\n",j);
	    //std::cout.flush();
	    std::vector<std::vector<long long> >av_gals_id,ss_av_gals_id,sf_av_gals_id;
	    int ret = snprintf(filename, cfilesize, "%s/save_av_gals_%05d.fits", F.outDir.c_str(), P[j].tileid);
	    read_save_av_gals(filename,  F,av_gals_id,ss_av_gals_id,sf_av_gals_id,diagnose);
	    //need to change potentialtargetid  place in list of mtl
	    for(int k=0;k<F.Nfiber;++k){
	      List collect;
	      for(int m=0;m<av_gals_id[k].size();++m){
		targetid_to_idx=invert_target.find(av_gals_id[k][m]);
		collect.push_back(targetid_to_idx->second);
	      }
	      av_gals.push_back(collect);
	      List collect_ss;
	      for(int m=0;m<ss_av_gals_id[k].size();++m){
		targetid_to_idx=invert_target.find(ss_av_gals_id[k][m]);
		collect_ss.push_back(targetid_to_idx->second);
	      }
	      ss_av_gals.push_back(collect_ss);
	      List collect_sf;
	      for(int m=0;m<sf_av_gals_id[k].size();++m){
		targetid_to_idx=invert_target.find(sf_av_gals_id[k][m]);
		collect_sf.push_back(targetid_to_idx->second);
	      }
	      sf_av_gals.push_back(collect_sf);
	     
	    }
	    P[j].av_gals=av_gals;
	    P[j].SS_av_gal_fiber=ss_av_gals;
	    P[j].SF_av_gal_fiber=sf_av_gals;
            //generate lsit by petal
            for(int k=0;k<F.Nfiber;++k){
	      int p=pp.spectrom[k];
	      for(int m=0;m<P[j].SS_av_gal_fiber[k].size();++m){
		P[j].SS_av_gal[p].push_back(P[j].SS_av_gal_fiber[k][m]);
	      }
	      for(int m=0;m<P[j].SF_av_gal_fiber[k].size();++m){
		P[j].SF_av_gal[p].push_back(P[j].SF_av_gal_fiber[k][m]);
	      }
	    }
	  }
	  

    //results_on_inputs("doc/figs/",G,P,F,true);
    }
    if(!infile.good() && savetime){//first epoch only, write files
        printf("no files there\n");
	std::cout.flush();
        for (int j=0; j<F.Nplate; j++){	  
	  write_save_av_gals(j,F.outDir,M,P,pp,F);
        }
    }else{
      printf("save files already there\n");
      std::cout.flush();
    }
    //Diagnostic comparison
    for(int j=0;j<F.Nplate;++j){
      int sum_av_gal=0;
      int sum_ss_av_gal=0;
      int sum_sf_av_gal=0;
      for (int k=0;k<F.Nfiber;++k){
	   sum_av_gal+=P[j].av_gals[k].size();
	   sum_ss_av_gal+=P[j].SS_av_gal_fiber[k].size();
	   sum_sf_av_gal+=P[j].SF_av_gal_fiber[k].size();
      }
      if (sum_av_gal>0)printf(" *** tile %d  sum_av_gal  %d  sum_ss_av_gal %d  sum_sf_av_gal %d\n",sum_av_gal,sum_ss_av_gal,sum_sf_av_gal);
	}
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


    // Smooth out distribution of free fibers, and increase the number of assignments
    
    for (int i=0; i<1; i++) redistribute_tf(M,P,pp,F,A,0);// more iterations will improve performance slightly
    for (int i=0; i<3; i++) {
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
       /* printf(" plate jused %5d j %5d  SS   %4d    SF   %4d",jused,j,used_SS,used_SF);
        for (int pr=0;pr<M.priority_list.size();++pr){
            printf(" class %2d   %5d",pr,used_by_class[pr]);
        }
        printf("\n");
        */
    }
    init_time_at(time,"# count SS and SF ",t);
    printf(" Totals SS   %4d    SF   %4d",total_used_SS,total_used_SF);
    std::cout.flush();
    for (int pr=0;pr<M.priority_list.size();++pr){
        printf(" class %2d   %5d",pr,total_used_by_class[pr]);
        std::cout.flush();
    }


    printf("\n");
    init_time_at(time,"# print txt files ",t);
    if (F.PrintAscii) for (int jused=0; jused<F.NUsedplate; jused++){
        int j=A.suborder[jused];
        write_FAtile_ascii(j,F.outDir,M,P,pp,F,A);
        }
    init_time_at(time,"# print fits files ",t);
    if (F.PrintFits) for (int jused=0; jused<F.NUsedplate; jused++){
        int j=A.suborder[jused];
        fa_write(j,F.outDir,M,P,pp,F,A); // Write output
    }

  
  print_time(t,"# Finished !... in");
  
  return(0);    
}
