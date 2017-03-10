#include    <cstdlib>
#include    <cmath>
#include    <cstdio>
#include    <fstream>
#include    <sstream>
#include    <iostream>
#include    <iomanip>
#include    <string>
#include    <string.h>
#include    <cstring>
#include    <vector>
#include    <algorithm>
#include    <exception>
#include    <stdexcept>
#include    <sys/time.h>
#include    <sys/stat.h>
#include    <stdlib.h>    

#include    "modules/htmTree.h"
#include    "modules/kdTree.h"
#include    "misc.h"
#include    "feat.h"
#include    "structs.h"
#include    "global.h"
    Time t; // t for global, time for local
    
// Collecting information from input -------------------------------------------------------------------------------------
void collect_galaxies_for_all(const MTL& M, const htmTree<struct target>& T, Plates& P, const PP& pp, const Feat& F) {
    //provides list of galaxies available to fiber k on tile j: P[j].av_gals[k]
    
    init_time(t,"# Begin collecting available galaxies");

    double rad = F.PlateRadius*M_PI/180.;

    #pragma omp parallel
    {
        #pragma omp master
        {
	    printf(" ");
        }
        // Collects for each plate
        // start at jj=0 not id
        #pragma omp for
        for (int j=0; j<F.Nplate; j++){
            plate p = P[j];
            // Takes neighboring galaxies that fall on this plate
            std::vector<int> nbr = T.near(M,p.nhat,rad);
            // Projects thoses galaxies on the focal plane
            Onplates O;
            for (int gg=0; gg<nbr.size(); gg++) {
                int g = nbr[gg];
                struct onplate op = change_coords(M[g],p);
                op.id = g;
		// Check that the target corresponds to the right program
		if(M[g].obsconditions & p.obsconditions){
		  O.push_back(op);
		}
            }
            // Build 2D KD tree of those galaxies
            KDtree<struct onplate> kdT(O,2);
            // For each fiber, finds all reachable galaxies within patrol radius, thanks to the tree
            for (int k=0; k<F.Nfiber; k++) {
                dpair X = pp.coords(k);
                std::vector<int> gals = kdT.near(&(pp.fp[2*k]),0.0,F.PatrolRad);
                for (int g=0; g<gals.size(); g++) {
                    dpair Xg = projection(gals[g],j,M,P);
                    if (sq(Xg,X)<sq(F.PatrolRad)){
                        P[j].av_gals[k].push_back(gals[g]);
                        int q=pp.spectrom[k];
                        //better to make SS & SF assigned to fibers

                        if(M[gals[g]].SS){
                            P[j].SS_av_gal[q].push_back(gals[g]);
                            P[j].SS_av_gal_fiber[k].push_back(gals[g]);
                        }
                        if(M[gals[g]].SF){
                            P[j].SF_av_gal[q].push_back(gals[g]);
                            P[j].SF_av_gal_fiber[k].push_back(gals[g]);
                        }

                    }
                }
            }
        }
    }

    print_time(t,"# ... took :");
}

bool pairCompare(const std::pair<double, int>& firstElem, const std::pair<double, int>& secondElem) {
	return firstElem.first < secondElem.first;//might want to reverse this
}

std::vector<int> sort_by_subpriority(MTL & M,std::vector<int> init){
    std::vector <int> out;
    std::vector<std::pair<double,int> > pairs;
    for(int gg=0;gg<init.size();++gg){
      int g=init[gg];
      std::pair <double,int> this_pair (M[g].subpriority,g);
      pairs.push_back(this_pair);
    }
    std::sort(pairs.begin(),pairs.end(),pairCompare);
    for(int gg=0;gg<init.size();++gg){
      out.push_back(pairs[gg].second);
    }
    return out;
}

void collect_available_tilefibers(MTL& M, const Plates& P, const Feat& F) { 
    //M[i].av_tfs is list of tile-fiber pairs available to galaxy i
    Time t;
    init_time(t,"# Begin computing available tilefibers");
    for(int j=0; j<F.Nplate; j++) {
        for(int k=0; k<F.Nfiber; k++) {
            for(int m=0; m<P[j].av_gals[k].size(); m++) {
                int i = P[j].av_gals[k][m];  //i is the id of the mth galaxy available to tile j and fiber k
                M[i].av_tfs.push_back(pair(j,k));  //list of tile-fibers available to galaxy i
            }
        }
    }
    print_time(t,"# ... took :");
}

// Assignment sub-functions -------------------------------------------------------------------------------------
// Allow (j,k) to observe g ?
inline bool ok_assign_g_to_jk(int g, int j, int k, const Plates& P, const MTL& M, const PP& pp, const Feat& F, const Assignment& A) {

    if (F.Collision) for (int i=0; i<pp.N[k].size(); i++) if (g==A.TF[j][pp.N[k][i]]) return false; 
    // Avoid 2 neighboring fibers observe the same galaxy (can happen only when Collision=true)
    if (A.find_collision(j,k,g,pp,M,P,F)!=-1){
        return false;} // No collision
    return true;
    //doesn't require that jk is unassigned//doesn't require that g isn't assigned already on this plate
    //use is_assigned_jg for this
}

// makes sure we don't exceed limit on SS and SF
inline bool ok_for_limit_SS_SF(int g, int j, int k, const MTL& M, const Plates& P, const PP& pp, const Feat& F){
     bool is_SF=M[g].SF;
    bool too_many_SF=P[j].SF_in_petal[pp.spectrom[k]]>F.MaxSF-1;
    bool is_SS=M[g].SS;
    bool too_many_SS=P[j].SS_in_petal[pp.spectrom[k]]>F.MaxSS-1;
    return !(is_SF && too_many_SF)&&!(is_SS && too_many_SS);
}

  

// Find, for (j,k), find the best galaxy it can reach among the possible ones
// Null list means you can take all possible kinds, otherwise you can only take, for the galaxy, a kind among this list
// Not allowed to take the galaxy of id no_g
inline int find_best(int j, int k, const MTL& M, const Plates& P, const PP& pp, const Feat& F, const Assignment& A) {        
    int best = -1; int mbest = -1; int pbest = 0; double subpbest = 0.;
    List av_gals = P[j].av_gals[k];
    // For all available galaxies
 
    for (int gg=0; gg<av_gals.size(); gg++) {
        int g = av_gals[gg];

        //if(ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){//don't assign SS, SF with find_best 11/20/15
        if(!M[g].SS && !M[g].SF){
            int m = M[g].nobs_remain; // Check whether it needs further observation
            if (m>=1) {
                int prio = M[g].t_priority;
		double subprio = M[g].subpriority;
                // Takes it if better priority, or if same, if it needs more observations, 
		//so shares observations if two QSOs are close
		// If still tied, use subpriority
                if (prio>pbest || (prio==pbest && m>mbest) || (prio==pbest && m==mbest && subprio>subpbest)){
                    // Check that g is not assigned yet on this plate, or on the InterPlate around, check with ok_to_assign
                    int isa=A.is_assigned_jg(j,g,M,F);
                    int ok=ok_assign_g_to_jk(g,j,k,P,M,pp,F,A);
                    if (isa==-1 && ok) {
                        best = g;
                        pbest = prio;
                        mbest = m;
			subpbest = subprio;

                    }
                }
            }
        }

    }
    return best;
}

// Tries to assign the fiber (j,k)
inline int assign_fiber(int j, int k, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    if (A.is_assigned_tf(j,k)) return -1;
    int best = find_best(j,k,M,P,pp,F,A);
    if (best!=-1) A.assign(j,k,best,M,P,pp);
    return best;
}



inline int assign_galaxy(int g,  MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int jstart) {
    // Tries to assign the galaxy g to one of the used plates after jstart
    //jstart runs possibly to F.Nplate
    int jb = -1; int kb = -1; int unusedb = -1;
    Plist av_tfs = M[g].av_tfs;
    // All the tile-fibers that can observe galaxy g
    for (int tfs=0; tfs<av_tfs.size(); tfs++) {
        int j = av_tfs[tfs].f;
        int k = av_tfs[tfs].s;
        // Check if the assignment is possible, if ok, if the tf is not used yet, and if the plate is in the list
        if (jstart<j && !A.is_assigned_tf(j,k) && ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)) {
            int unused = A.unused[j][pp.spectrom[k]];//unused fibers on this petal
            if (unusedb<unused) {
                jb = j; kb = k; unusedb = unused;//observe this galaxy by fiber on petal with most free fibefs
            }
        }
    }
    if (jb!=-1){
        A.assign(jb,kb,g,M,P,pp);
        return 1;}
    else return 0;
}

// Takes an unassigned fiber and tries to assign it with the "improve" technique described in the doc
// not used for SS or SF   
inline int improve_fiber(int jused_begin, int jused, int k, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1) {

    int j=A.suborder[jused];
    if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
        int g_try = assign_fiber(j,k,M,P,pp,F,A);//maybe doesn't allow SS or SF
        if (g_try!=-1) return g_try;
        else { // Improve
            int gb = -1; int bb = -1; int jpb = -1; int kpb = -1; int mb = -1; int pb = 1e3; int unusedb = -1;
	    std::vector<int>  av_g_init = P[j].av_gals[k];
	    std::vector<int> av_g;
            // For all available galaxies within reach that are already observed
	    av_g=sort_by_subpriority(M,av_g_init);
            for (int i=0; i<av_g.size(); i++) {
                int g = av_g[i];//a galaxy accessible to j,k

                if (g!=-1 && g!=no_g && !M[g].SS && !M[g].SF) {//not SS or SF
                    if (ok_assign_g_to_jk(g,j,k,P,M,pp,F,A) ) {
                        // Which tile-fibers have taken g ?
                        Plist tfs = A.chosen_tfs(g,F,A.suborder[jused_begin]);//all tile-fibers that observe g in tiles from begin to end
                        for (int p=0; p<tfs.size(); p++) {
                            int jp = tfs[p].f;
                            int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
                            // FIND BEST JP KP !!!
                            int best = find_best(jp,kp,M,P,pp,F,A); // best!=g because !A.assigned_pg(best)

                            if (best!=-1 && (A.is_assigned_jg(j,g,M,F)==-1 || jp==j)) {
                                int prio = M[g].t_priority;
                                int m = M[g].nobs_remain;
                                int unused = A.unused[jp][pp.spectrom[kp]]; // We take the most unused
                                if (prio>pb || (prio==pb && m>mb) || (prio==pb && m==mb && unused>unusedb)) {
                                    gb = g; bb = best; jpb = jp; kpb = kp; mb = m; pb = prio; unusedb = unused;
                            }}}}}}
            // Modify assignment
            if (gb!=-1) {
                A.unassign(jpb,kpb,gb,M,P,pp);
                A.assign(j,k,gb,M,P,pp);
                A.assign(jpb,kpb,bb,M,P,pp);
                return gb;
            }
        }
    }
    return -1;
}
// Assignment functions ------------------------------------------------------------------------------------------
// Assign fibers naively

void simple_assign(MTL &M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    Time t;
    init_time(t,"# Begin simple assignment :");
    int countme=0;
    for (int j=0; j<F.Nplate; j++) {

        int best=-1;
    for (int k=0; k<F.Nfiber; k++) { // Fiber
      best=assign_fiber(j,k,M,P,pp,F,A);
      if (best!=-1)countme++;
    }
    }
    print_time(t,"# ... took :");
    printf(" countme %d \n",countme);
}

void improve( MTL& M, Plates&P, const PP& pp, const Feat& F, Assignment& A, int jused_start) {
    //jstart is in list from 0 to F.NUsedplate
    Time t;
    init_time(t,"# Begin improve :");
    int improvements=0;
    for (int jused=jused_start; jused<F.NUsedplate; jused++){
        for (int k=0; k<F.Nfiber; k++){
            int worked=improve_fiber(jused_start,jused,k,M,P,pp,F,A);
            if (worked!=-1)improvements++;
        }
    }
    printf(" improvements  %d\n",improvements);
    print_time(t,"# ... took :");
}

// If there are galaxies discovered as fake for example, they won't be observed several times in the plan
// has access to G,not just M, because it needs to know the truth



void new_replace( int j, int p, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    // do standard stars,going through priority classes from least to most
    //keep track of reassignments

    int reassign_SS=0;
    int reassign_SF=0;
    int add_SS=0;
    int add_SF=0;
   
    //can get all available SS,SF on plate from P[j].av_gals_plate restricting to plate p

    for(int c=0;P[j].SS_in_petal[p]<F.MaxSS && c<M.priority_list.size();++c ){//try to do this for lowest priority first
        std::vector <int> gals_init=P[j].SS_av_gal[p]; //standard stars on this petal
	std::vector <int> gals;
	bool thisplate=false;
	//sort these by priority
	gals=sort_by_subpriority(M,gals_init);
        for(int gg=0;gg<gals.size() ;++gg){//all the standard stars on this petal

            int g=gals[gg];//a standard star
	    bool thisstar=false;

            if(A.is_assigned_jg(j,g)==-1){//not assigned on this tile
                Plist tfs=M[g].av_tfs;//all tiles and fibers that reach g

                int done=0;//quit after we've used this SS
		std::vector<int> could_be_replaced;
		std::vector<int> sorted_could_be_replaced;
                for(int i=0;i<tfs.size() && done==0;++i){//al the tile-fibers that can reach this standard star

                    if(tfs[i].f==j&&pp.spectrom[tfs[i].s]==p){//a tile fiber from this petal
                        int k=tfs[i].s;//we know g can be reached by this petal of plate j and fiber k
                        int g_old=A.TF[j][k];//what is now at (j,k)  g_old can't be -1 or we would have used it already in assign_sf
			//require that this galaxy is used only once to keep things simple
                        if(g_old!=-1 && !M[g_old].SS && !M[g_old].SF && A.GL[g_old].size()==1){
                            if (M[g_old].priority_class==c&&A.is_assigned_jg(j,g,M,F)==-1 && ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                                //right priority; this SS not already assigned on this plate
			      could_be_replaced.push_back(g_old);
			    }
			}
		    }
		}
		if(could_be_replaced.size()>0){
		  //now have list of possible targets to be replaced with lowest priority
		  //use subpriority to choose
		  std::vector<int> sorted_could_be_replaced=sort_by_subpriority(M,could_be_replaced);
		  int g_chosen=sorted_could_be_replaced[0];
		  int j=A.GL[g_chosen][0].f;
		  int k=A.GL[g_chosen][0].s;
		  
		  A.unassign(j,k,g_chosen,M,P,pp);
		  reassign_SS+=assign_galaxy(g_chosen,M,P,pp,F,A,j);//try to assign, increment counter 
		  A.assign(j,k,g,M,P,pp);
		  add_SS++;
		  done=1;	

		}
	    }
	}
    }

    for(int c=0;P[j].SF_in_petal[p]<F.MaxSF && c<M.priority_list.size();++c ){//try to do this for lowest priority
        std::vector <int> gals_init=P[j].SF_av_gal[p]; //sky fibers on this petak
	bool thisplate=false;
	//sort these by subpriority
	std::vector<std::pair<double,int> > galaxy_pairs;
        for(int gg=0;gg<gals_init.size() ;++gg){
	  int g=gals_init[gg];
	  std::pair <double,int> this_pair (M[g].subpriority,g);
	  galaxy_pairs.push_back(this_pair);
	}
	std::sort(galaxy_pairs.begin(),galaxy_pairs.end(),pairCompare);
	std::vector <int> gals;

	for(int gg=0;gg<gals_init.size();++gg){
	  gals.push_back(galaxy_pairs[gg].second);
	}
        for(int gg=0;gg<gals.size();++gg){
            int g=gals[gg];//a sky fiber
            if(A.is_assigned_jg(j,g)==-1){
                Plist tfs=M[g].av_tfs;
                int done=0;
		std::vector<int> could_be_replaced; 
                for(int i=0;i<tfs.size() && done==0;++i){

                    if(tfs[i].f==j&&pp.spectrom[tfs[i].s]==p){//g is accesible to j,k
                        int k=tfs[i].s;//we know g can be reached by this petal of plate j and fiber k
                        int g_old=A.TF[j][k];//what is now at (j,k)

                        if(g_old!=-1 && !M[g_old].SS && !M[g_old].SF && A.GL[g_old].size()==1){
                            if (M[g_old].priority_class==c&&A.is_assigned_jg(j,g,M,F)==-1 && ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
			      could_be_replaced.push_back(g_old);
			    }
			}
		    }
		}
		//now have list of possible targets to be replaced with lowest priority
		//use subpriority to choose
		sort_by_subpriority(M,could_be_replaced);
		if(could_be_replaced.size()>0){
		  int g_chosen=could_be_replaced[0];
		  int j=A.GL[g_chosen][0].f;
		  int k=A.GL[g_chosen][0].s;
		  A.unassign(j,k,g_chosen,M,P,pp);
		  reassign_SF+=assign_galaxy(g_chosen,M,P,pp,F,A,j);//try to assign
		  A.assign(j,k,g,M,P,pp);
		  add_SF++;
		  done=1;

		}
	    }
	}
    }
}




void assign_unused(int j, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    // Tries to assign remaining fibers in tile jth tile with galaxies on it
    //even taking objects observed later
    //js is a tile with galaxies on it
   
    for (int k=0; k<F.Nfiber; k++) {
        
        if (!A.is_assigned_tf(j,k)) {
	  int best = -1; int mbest = -1; int pbest = 0; int jpb = -1; int kpb = -1; 
	  double subpbest = 0;
	  std::vector<int> av_gals = P[j].av_gals[k];//all available galaxies for this fiber k
	  std::vector<int> sorted_av_gals=sort_by_subpriority(M,av_gals);
            for (int gg=0; gg<sorted_av_gals.size(); gg++) {
                int g = sorted_av_gals[gg];//available galaxies
                int m = M[g].nobs_remain;
                int prio = M[g].t_priority;
		double subprio = M[g].subpriority;
                if (prio>pbest || (prio==pbest && m>mbest) || (prio==pbest && m==mbest && subprio>subpbest)) {
                    if (A.is_assigned_jg(j,g,M,F)==-1 && ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                        //not assigned this plate or within excluded interval
                        for (int i=0; i<A.GL[g].size(); i++) { //GL[g].size() is number of tf that could observe g
                            int jp = A.GL[g][i].f;
                            int kp = A.GL[g][i].s;
                            if (j<jp && jpb<jp) {//take best opportunity
                                best = g;
                                pbest = prio;
                                mbest = m;
				subpbest = subprio;
                                jpb = jp;
                                kpb = kp;
                            }
                        }
                    }
                }
            }
            if (best!=-1) {
                A.unassign(jpb,kpb,best,M,P,pp);
                A.assign(j,k,best,M,P,pp);
                
            }
        }
    }
}


void assign_sf_ss(int j, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    //assigns sky fibers and standard stars to unused fibers
    bool thisgalaxy=false;
    for (int ppet=0; ppet<F.Npetal; ppet++) {
      //for each petal
        int p = ppet;
        std::vector <int> SS_av_init=P[j].SS_av_gal[p];
	std::vector <int> SS_av;
        SS_av=sort_by_subpriority(M,SS_av_init);
	//use subpriority to order SS, SF

        std::vector <int> SF_av_init=P[j].SF_av_gal[p];
        std::vector <int> SF_av;
        SF_av=sort_by_subpriority(M,SF_av_init);
	
        if(SS_av.size()>0 ||SF_av.size()>0){
            //look at fibers on this petal
            for (int kk=0; kk<F.Nfbp; kk++) {
                int k= pp.fibers_of_sp[p][kk];
		//for a particular fiber on this petal consider SS and SF by subpriority order
                std::vector <int> SS_av_k_init=P[j].SS_av_gal_fiber[k];
                std::vector <int> SF_av_k_init=P[j].SF_av_gal_fiber[k];
		std::vector <int> SS_av_k;
		std::vector <int> SF_av_k;
		SS_av_k=sort_by_subpriority(M,SS_av_k_init);
		SF_av_k=sort_by_subpriority(M,SF_av_k_init);		
                if (A.TF[j][k]==-1){//this fiber isn't assigned yet
                    int done=0;
		    //
                    for (int gg=0; gg<SS_av_k.size()&&done==0; gg++) {
                        int g = SS_av_k[gg];//SS on this with reach of this fiber	
			// if g isn't already assigned on this fiber and we still need it and it's ok to assign, then assign
                        if(A.is_assigned_jg(j,g,M,F)==-1&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)&&ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)){
                            A.assign(j,k,g,M,P,pp);
                            done=1;			   
                        }
                    }
                    for (int gg=0; gg<SF_av_k.size()&&done==0; gg++) {
                        int g = SF_av_k[gg];//SF on this petal						
                        if(A.is_assigned_jg(j,g,M,F)==-1&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)&&ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)){
                            A.assign(j,k,g,M,P,pp);
                            done=1;			   
                        }
                    }
                }
            }
        new_replace(j,p,M,P,pp,F,A);
        }
    }

}

void redistribute_tf(MTL& M, Plates&P, const PP& pp, const Feat& F, Assignment& A, int jused_start) {
    //diagnostic
    printf("start redistribute \n");
    Time t;
    init_time(t,"# Begin redistribute TF :");
    int red(0);
    Table Done = initTable(F.NUsedplate,F.Nfiber);//consider every occupied plate and every fiber
    for (int jused=jused_start; jused<F.NUsedplate; jused++) {
        int j=A.suborder[jused];
        for (int k=0; k<F.Nfiber; k++) {
            if (Done[jused][k]==0) {
                int g = A.TF[j][k];//current assignment of (j,k)  only look if assigned
                if (g!=-1 && !M[g].SS && !M[g].SF) {
                    int jpb = -1; int kpb = -1; int unusedb = A.unused[j][pp.spectrom[k]];
                    Plist av_tfs = M[g].av_tfs;  //all possible tile fibers for this galaxy
                    for (int i=0; i<av_tfs.size(); i++) {
                        int jp = av_tfs[i].f;
                        int kp = av_tfs[i].s;
                        int unused = A.unused[jp][pp.spectrom[kp]];//unused for jp, spectrom[kp]
                        if(A.inv_order[jp]!=-1)//necessary because underdense targets may leave some plates unused
                        {
                        if(A.inv_order[jp]>F.NUsedplate || A.inv_order[jp]<0)printf("**out range  %d\n",A.inv_order[jp]);
                        if (A.suborder[jused_start]<=jp){
                            if(!A.is_assigned_tf(jp,kp)){
                                if(Done[A.inv_order[jp]][kp]==0){
                                    if( ok_assign_g_to_jk(g,jp,kp,P,M,pp,F,A)){
                                        if(A.is_assigned_jg(jp,g,M,F)==-1){
                                            if( 0<unused) {
                                                if (unusedb<unused) { // Takes the most unused petal
                                                    jpb = jp;
                                                    kpb = kp;
                                                    unusedb = unused;
                                               }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        }
                    }
                    if (jpb!=-1) {
                        A.unassign(j,k,g,M,P,pp);
                        A.assign(jpb,kpb,g,M,P,pp);
                        Done[A.inv_order[j]][k] = 1;
                        Done[A.inv_order[jpb]][kpb] = 1;
                        red++; 
                    }
                }
            }
        }
    }
    printf("  %s redistributions of tile-fibers \n",f(red).c_str());
    print_time(t,"# ... took :");
}


void diagnostic(const MTL& M, const Gals& Secret, Feat& F, const Assignment& A){
    // diagnostic  allows us to peek at the actual id of each galaxy
    printf("Diagnostics using types:QSO-Ly-a, QSO-tracers, LRG, ELG, fake QSO, fake LRG, SS, SF\n");
    std::vector<int> count_by_kind(F.Categories-2,0);
    for (int j=0;j<F.NUsedplate;++j){
        int js=A.suborder[j];
        //printf(" js = %d\n",js);
        //printf(" Secret size %d\n",Secret.size());
        for(int k=0;k<F.Nfiber;++k){
            int g=A.TF[js][k];
            if(g!=-1&&!M[g].SS&&!M[g].SF){
                //printf("g = %d  k = %d  id = %d \n",g,k,Secret[g].category);
            count_by_kind[Secret[g].category]+=1;
            }
        }
    }
    for(int i=0;i<F.Categories-2;++i){
        printf(" i  %d    number  %d \n",i,count_by_kind[i]);
    }
    int MaxObs = max(F.goal);
    Table obsrv = initTable(F.Categories,MaxObs+1);
    
    for (int g=0; g<M.size(); g++) {
        if(!M[g].SS && !M[g].SF){
        int c= Secret[g].category;
        int m = min(M[g].nobs_done,MaxObs);
        obsrv[c][m]++; //
        }
    }
    for (int c=0;c<F.Categories-2;++c){
        int tot=0;
        for (int m=0;m<MaxObs+1;++m){
            tot+=obsrv[c][m];
        }
        for (int m=0;m<MaxObs+1;++m){
            double ratio=float(obsrv[c][m])/float(tot);
            printf("   %f  ",ratio);
        }
        printf("\n");
    }
       //end diagnostic
}

void display_results(str outdir, const Gals& Secret,const MTL& M, const Plates& P, const PP& pp, Feat& F, const Assignment& A,  int last_tile,bool latex) {
	printf("# Results :\n");


    // 1 Raw numbers of galaxies by id and number of remaining observations
    int MaxObs = max(F.goal);
    Table obsrv = initTable(F.Categories-2,MaxObs+1);

    for (int g=0; g<M.size(); g++) {
        if(!M[g].SS && !M[g].SF){
		int c= Secret[g].category;
            int m=0;
            for(int k=0;k<A.GL[g].size();++k){
                if(A.GL[g][k].f<last_tile)++m;
            }
        obsrv[c][m]++; //
        }
    }

    // Add the 3 columns of tot, fibs, obs
    Table with_tots = obsrv;
    for (int i=0; i<F.Categories-2; i++) {
        int fibs = 0; int obs = 0; int tot =0;
        for (int j=0; j<=MaxObs; j++) tot += obsrv[i][j];
        for (int j=0; j<=MaxObs; j++) fibs += obsrv[i][j]*j;
        for (int j=1; j<=MaxObs; j++) obs += obsrv[i][j];
        with_tots[i].push_back(tot);
        with_tots[i].push_back(fibs);
        with_tots[i].push_back(obs);
    }

	//print_table("  Remaining observations (without negative obs ones)",with_tots,latex,F.kind);
	Dtable obs_per_sqd = ddivide_floor(with_tots,F.TotalArea);

	// Add percentages of observation
	Dtable perc = initDtable(F.Categories-2,2);
	for (int id=0; id<F.Categories-2; id++) {
		int tot = sumlist(obsrv[id]);
		int goal = F.goal[id];

		perc[id][0] = percent(tot-obsrv[id][0],tot);

		// Weighted percentage
		int d = 0;
		for (int i=0; i<=goal; i++) d += obsrv[id][i]*i;
		perc[id][1] = percent(d,tot*goal);
	}
	print_table("Obs per sqd and percentages",concatenate(obs_per_sqd,perc),latex,F.kind);

	// 3 Observed galaxies in function of time
	// Lya 1,2,3,4,5, LRG 1,2
	if (F.PlotObsTime) {
	int interval = 100;
	int nk = 9;
	Table Ttim = initTable(nk,0);

	List galaxs = initList(F.Ngal);
    printf(" F.Ntarg = %d\n",F.Ntarg);

	for (int jused=0; jused<F.NUsedplate; jused++) {
        int j=A.suborder[jused];

		for (int k=0; k<F.Nfiber; k++) {
			int g = A.TF[j][k];
            if (g!=-1) galaxs[g]++;
		}
        //only display at certain interval probably should be set in features file
		if (jused%interval==0) {
			List l = initList(20);//was 9
			for (int g=0; g<F.Ntarg; g++) {//assumes targets have lowest numbers!
				int n = galaxs[g];
                if (1<=n) {
                    if (Secret[g].category == 0){
                        l[n-1]++;
                        if(n>5)printf(" QSO Ly-a observed %d times\n",n);
                    }
                    if (Secret[g].category == 2){
                        l[n-1+5]++;
                        if(n>2)printf(" LRG observed %d times\n",n);
                    }
                    if (Secret[g].category == 1) {
                        l[n+6]++;
                    }
                    if (Secret[g].category == 3){
                        l[n+7]++;
                     
                    }
				}
         
			}
			for (int id=0; id<nk; id++) Ttim[id].push_back(l[id]);
        }         
          
    }
        print_mult_table_latex("Observed galaxies complete (interval 100)",outdir+"time2.dat",Ttim,interval);
        printf("done with time table\n");
    }

	// 4 Histogram of percentages of seen Ly-a
	if (F.PlotHistLya) {
	int id = F.ids.at("QSOLy-a");
	int goal = F.goal[id];
	Table Percseen = initTable(goal+1,0);
	for (int g=0; g<F.Ntarg; g++) {
		if (M[g].id==id) {
			int n = M[g].av_tfs.size();
			int p = A.chosen_tfs(g,F).size();
			if (n>=Percseen[p].size()) Percseen[p].resize(n+1);
			Percseen[p][n]++;
		}
	}
	make_square(Percseen);
	//print_table("Number of QSO Ly-a : x - Number of available TF - y - Number of observations",Percseen);
	for (int j=0; j<Percseen[0].size(); j++) {
		for (int i=Percseen.size()-1; i!=0; i--) {
			Percseen[i-1][j] += Percseen[i][j];
		}
	}
	print_mult_table_latex("Available tile-fibers for a galaxy (by kind)",outdir+"obsly.dat",Percseen,1);
    }

	// 5 Histogram of time between 2 obs of Ly a
	if (F.PlotDistLya) {
	Table deltas;
	for (int g=0; g<F.Ntarg; g++) {
		if (M[g].id == F.ids.at("QSOLy-a")) {
			Plist tfs = A.chosen_tfs(g,F);
			if (tfs.size()>=2) {
				List unsorted;
				List del;
				for (int i=0; i<tfs.size(); i++) {
					unsorted.push_back(tfs[i].f);
				}
				List sorted = sort(unsorted);
				for (int i=0; i<sorted.size()-1; i++) {
					int p1 = sorted[i];
					int p2 = sorted[i+1];
					del.push_back(p2-p1);
				}
				deltas.push_back(del);
			}
		}
	}
	List histo0 = histogram(deltas,10);
	//print_hist("Plate interval between 2 consecutive obs of Ly-a (interval 100)",100,histogram(deltas,100));
	Table delts; delts.push_back(histo0); delts.push_back(cumulate(histo0));
	print_mult_table_latex("Plate interval between 2 consecutive obs of Ly-a (interval 10)",outdir+"dist2ly.dat",delts,10);
    }

	// 6 Free fibers histogram
	if (F.PlotFreeFibHist) {
	Table unused_fbp = A.unused_fbp(pp,F);
	make_square(unused_fbp);
	Table hist0; hist0.push_back(histogram(unused_fbp,1));
	print_mult_table_latex("Number of petals with this many free fiber (interval 1)",outdir+"freefib.dat",hist0,1);
	}

	// 7 Free fibers in function of time (plates)
	if (F.PlotFreeFibTime) {
	List freefibtime = initList(F.Nplate);
	for (int j=0; j<F.Nplate; j++) freefibtime[j] = A.unused_f(j,F);
	Table fft; fft.push_back(freefibtime);
	print_mult_table_latex("Free fibers in function of time (plates)",outdir+"fft.dat",fft);
    }
 
	// 8 Percentage of seen objects as a function of density of objects
	if (F.PlotSeenDens) {
	Dcube densities = initDcube(F.Categories+-21,0,0);
	for (int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			// For all
			int size = P[j].av_gals[k].size();
			int oc = 0;
			for (int i=0; i<size; i++) if (A.is_assigned_jg(j,P[j].av_gals[k][i])!=-1) oc++;
			if (size!=0 && 1<=oc) { 
				double d = percent(oc,size);
				//printf("%f %f %f %d %d %d \n",d,x,invFibArea,size,oc,densities.size()); fl();
				if (size>=densities[F.Categories].size()) densities[F.Categories].resize(size+1);
				densities[F.Categories][size].push_back(d);
			}

			// For kind
			for (int t=0; t<F.Categories-2; t++) {
				int nkind = 0;
				int ock = 0;
				for (int i=0; i<size; i++) {
					int g = P[j].av_gals[k][i];
					if (M[g].id == t) {
						nkind++;
						if (A.is_assigned_jg(j,g)!=-1) ock++;
					}
				}
				if (nkind!=0 && 1<=ock) { 
					double d = percent(ock,nkind);
					if (nkind>=densities[t].size()) densities[t].resize(nkind+1);
					densities[t][nkind].push_back(d);
				}
			}
		}
	}
	Dtable densit = initDtable(F.Categories-2+1,max_row(densities));
	for (int t=0; t<F.Categories-2+1; t++) for (int i=0; i<densities[t].size(); i++) densit[t][i] = sumlist(densities[t][i])/densities[t][i].size();
	print_mult_Dtable_latex("Perc of seen obj as a fun of dens of objs",outdir+"seendens.dat",densit,1);
    }
	
	// 9 Collision histogram of distances between galaxies
	if (F.Collision) {
	Dlist coldist;
	for (int j=0; j<F.Nplate; j++) {
		List done = initList(F.Nfiber);
		for (int k=0; k<F.Nfiber; k++) {
			if (done[k]==0) {
				int c = A.is_collision(j,k,pp,M,P,F);
				if (c!=-1) {
					done[c] = 1;
					dpair G1 = projection(A.TF[j][k],j,M,P);
					dpair G2 = projection(A.TF[j][c],j,M,P);
					double d = norm(G2-G1);
					coldist.push_back(d);
				}
			}
		}
	}
	double intervaldist = 0.01;
	
	Dlist histcoldist = histogram(coldist,intervaldist);
	Dlist redhistcol = percents(histcoldist,sumlist(histcoldist));
	Dtable Dtd; Dtd.push_back(redhistcol); Dtd.push_back(cumulate(redhistcol));
	print_mult_Dtable_latex("Collision histogram of distances between galaxies",outdir+"coldist.dat",Dtd,intervaldist);
    }

	// Collision rate
	if (F.Collision) printf("Collision rate : %f %% \n",A.colrate(pp,M,P,F));

	// Percentage of fibers assigned
	//printf("  %s assignments in total (%.4f %% of all fibers)\n",f(A.na(F)).c_str(),percent(A.na(F),F.Nplate*F.Nfiber));

	// Count
	if (F.Count!=0) printf("Count = %d \n",F.Count);
    // print no. of times each galaxy is observed up to max of F.PrintGalObs
    if (F.PrintGalObs>0){
        printf(" F.PrintGalObs  %d \n",F.PrintGalObs);
        for(int g=0;g<F.PrintGalObs;++g){
                int id = M[g].id;
                int m = M[g].nobs_remain;
                int n = F.goal[id]-m;
            printf(" galaxy number %d  times observed %d\n",g,n);
        }
    }
}



void fa_write (int j, str outdir, const MTL & M, const Plates & P, const PP & pp, const Feat & F, const Assignment & A) {
    
    // generate a quiet NaN to use for invalid entries.  We cannot
    // guarantee that we have C++11, so we can't use the nice functions
    // included in that standard...
    
    const unsigned maxU = ~0;
    const float qNan = *((float*)&maxU);
    
    // constants for the filename length and fixed object
    // type length
    
    size_t cfilesize = 512;
    size_t objtypelen = 8;
    size_t bricklen = 8;
    // check if the file exists, and if so, throw an exception
    
    char filename[cfilesize];
    int ret = snprintf(filename, cfilesize, "%s/tile_%05d.fits", outdir.c_str(), P[j].tileid);
    
    struct stat filestat;
    ret = ::stat(filename, &filestat );
    
    if (ret == 0) {
        std::ostringstream o;
        o << "output file " << filename << " already exists";
        throw std::runtime_error(o.str().c_str());
    }
    
    // create the file
    
    int status = 0;
    fitsfile * fptr;
    fits_create_file (&fptr, filename, &status);
    fits_report_error (stderr, status);
    
    // Set up the schema for the table.  We explicitly malloc these
    // string arrays, since the CFITSIO API requires non-const pointers
    // to them (i.e. arrays of literals won't work).
    
    size_t ncols = 13;
    
    char ** ttype;
    char ** tform;
    char ** tunit;
    
    ttype = (char**) malloc ( ncols * sizeof(char*) );
    tform = (char**) malloc ( ncols * sizeof(char*) );
    tunit = (char**) malloc ( ncols * sizeof(char*) );
    
    if ( ! ( ttype && tform && tunit ) ) {
        std::ostringstream o;
        o << "cannot allocate column info for binary table";
        throw std::runtime_error(o.str().c_str());
    }
    
    for ( size_t c = 0; c < ncols; ++c ) {
        ttype[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        tform[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        tunit[c] = (char*) malloc ( FLEN_VALUE * sizeof(char) );
        if ( ! ( ttype[c] && tform[c] && tunit[c] ) ) {
            std::ostringstream o;
            o << "cannot allocate column info for binary table";
            throw std::runtime_error(o.str().c_str());
        }
    }
    
    strcpy(ttype[0], "FIBER");
    strcpy(tform[0], "J");
    strcpy(tunit[0], "");
    
    strcpy(ttype[1], "POSITIONER");
    strcpy(tform[1], "J");
    strcpy(tunit[1], "");
    
    strcpy(ttype[2], "NUMTARGET");
    strcpy(tform[2], "I");//int not long
    strcpy(tunit[2], "");
 
    strcpy(ttype[3], "PRIORITY");
    strcpy(tform[3], "J");
    strcpy(tunit[3], "");

    strcpy(ttype[4], "TARGETID");
    strcpy(tform[4], "K");
    strcpy(tunit[4], "");
    
    strcpy(ttype[5], "DESI_TARGET");
    strcpy(tform[5], "K");
    strcpy(tunit[5], "");
    
    strcpy(ttype[6], "BGS_TARGET");
    strcpy(tform[6], "K");
    strcpy(tunit[6], "");
    
    strcpy(ttype[7], "MWS_TARGET");
    strcpy(tform[7], "K");
    strcpy(tunit[7], "");
    
    strcpy(ttype[8], "RA");
    strcpy(tform[8], "D");
    strcpy(tunit[8], "deg");
    
    strcpy(ttype[9], "DEC");
    strcpy(tform[9], "D");
    strcpy(tunit[9], "deg");
    
    strcpy(ttype[10], "XFOCAL_DESIGN");
    strcpy(tform[10], "E");
    strcpy(tunit[10], "mm");
    
    strcpy(ttype[11], "YFOCAL_DESIGN");
    strcpy(tform[11], "E");
    strcpy(tunit[11], "mm");

    strcpy(ttype[12], "BRICKNAME");
    snprintf(tform[12], FLEN_VALUE, "%dA", (int)bricklen);
    strcpy(tunit[12], "");

    
    char extname[FLEN_VALUE];
    
    strcpy(extname, "FIBER_ASSIGNMENTS");
    
    // create the table with the full size on disk.
    
    ret = fits_create_tbl(fptr, BINARY_TBL, F.Nfiber, ncols, ttype, tform, tunit, extname, &status);
    fits_report_error(stderr, status);
    
    // get the number of rows to write for each internal FITS buffer.
    
    long optimal;
    ret = fits_get_rowsize(fptr, &optimal, &status);
    fits_report_error(stderr, status);
    
    // initialize arrays to the optimal number of rows for writing.
    int tile_id[optimal];
    int fiber_id[optimal];
    int positioner_id[optimal];
    int num_target[optimal];
    char objtype[optimal][objtypelen];
    char brickname[optimal][bricklen+1];
    char * bn_tmp[optimal];
    char * ot_tmp[optimal];
    for (int i = 0; i < optimal; i++) {
        ot_tmp[i] = objtype[i];
	bn_tmp[i] = brickname[i];
    }
    long long target_id[optimal];
    long long desi_target[optimal];
    long long bgs_target[optimal];
    long long mws_target[optimal];
    float ra[optimal];
    float dec[optimal];
    float x_focal[optimal];
    float y_focal[optimal];
    //new
    int t_priority[optimal];
    
    std::vector <long long> potentialtargetid;
 
    
    long long offset = 0;
    long long n = optimal;
    
    while ( n == optimal ) {
        
        if ( offset + optimal > F.Nfiber ) {
            n = F.Nfiber - offset;
        }
        
        if ( n > 0 ) {
            
            for (int i = 0; i < n; ++i) {
                int fib = offset + i;
                int g = A.TF[j][fib];
                
                fiber_id[i] = fib;
                positioner_id[i] = fib;
                num_target[i] = P[j].av_gals[fib].size();


                //target_id[i] = g; ********
                if(g>0) target_id[i] = M[g].id;
                else target_id[i]=-1;

                if (g < 0) {
                    //strcpy(objtype[i], "NA");
                    ra[i] = qNan;
                    dec[i] = qNan;
                    x_focal[i] = qNan;
                    y_focal[i] = qNan;
                    desi_target[i] = 0;
                    bgs_target[i] = 0;
                    mws_target[i] = 0;
		    strncpy(brickname[i], "notbrick", bricklen+1);
                } else {
                    //we aren't supposed to know the kind  use priority instead
                    //strncpy(objtype[i], F.kind[G[g].id].c_str(), objtypelen);
                    ra[i] = M[g].ra;
                    dec[i] = M[g].dec;
                    dpair proj = projection(g,j,M,P);
                    x_focal[i] = proj.f;
                    y_focal[i] = proj.s;
                    t_priority[i]=M[g].t_priority;//new
                    desi_target[i] = M[g].desi_target;
                    bgs_target[i] = M[g].bgs_target;
                    mws_target[i] = M[g].mws_target;
		    strncpy(brickname[i], M[g].brickname, bricklen+1);
                }

                // Store the potential targetids accesible to this fibre (the actual targetid, not the index).
                for (int k = 0; k < P[j].av_gals[fib].size(); ++k) {
                    int gal_idx = P[j].av_gals[fib][k]; // MTL index for k'th target accessible to this fibre
                    if (gal_idx >= 0) {
                        potentialtargetid.push_back(M[gal_idx].id);
                }}
            }
        int tileid = P[j].tileid;
        float tilera = P[j].tilera;
        float tiledec = P[j].tiledec;

        fits_write_key(fptr, TINT, "TILEID", &(tileid), "Tile ID number", &status);
            fits_report_error(stderr, status);
        fits_write_key(fptr, TFLOAT, "TILERA", &(tilera), "Tile RA [deg]", &status);
        fits_report_error(stderr, status);
        fits_write_key(fptr, TFLOAT, "TILEDEC", &(tiledec), "Tile DEC [deg]", &status);
            fits_report_error(stderr, status);

            fits_write_col(fptr, TINT, 1, offset+1, 1, n, fiber_id, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TINT, 2, offset+1, 1, n, positioner_id, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TINT, 3, offset+1, 1, n, num_target, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TINT, 4, offset+1, 1, n, t_priority, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TLONGLONG, 5, offset+1, 1, n, target_id, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TLONGLONG, 6, offset+1, 1, n, desi_target, &status);
            fits_report_error(stderr, status);

            fits_write_col(fptr, TLONGLONG, 7, offset+1, 1, n, bgs_target, &status);
            fits_report_error(stderr, status);

            fits_write_col(fptr, TLONGLONG, 8, offset+1, 1, n, mws_target, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 9, offset+1, 1, n, ra, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 10, offset+1, 1, n, dec, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 11, offset+1, 1, n, x_focal, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 12, offset+1, 1, n, y_focal, &status);
            fits_report_error(stderr, status);

	    fits_write_col(fptr, TSTRING, 13, offset+1, 1, n, bn_tmp, &status);
	    fits_report_error(stderr, status);
        }
        
        offset += n;
    }
    
    // PotentialFiberMap table.  We have only one column, so it is safe
    // from a performance perspective to write the whole thing.
    //NO  12/15/16
    
    strcpy(ttype[0], "POTENTIALTARGETID");
    strcpy(tform[0], "K");
    strcpy(tunit[0], "");
    
    strcpy(extname, "POTENTIAL_ASSIGNMENTS");
    
    ret = fits_create_tbl(fptr, BINARY_TBL, potentialtargetid.size(), 1, ttype, tform, tunit, extname, &status);
    fits_report_error(stderr, status);
    
    fits_write_col(fptr, TLONGLONG, 1, 1, 1, potentialtargetid.size(), &(potentialtargetid[0]), &status);
    fits_report_error(stderr, status);
    
    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);
    
    for ( size_t c = 0; c < ncols; ++c ) {
        free ( ttype[c] );
        free ( tform[c] );
        free ( tunit[c] );
    }
    free ( ttype );
    free ( tform );
    free ( tunit );
    
    return;
}


void pyplotTile(int jused, str directory, const Gals& Secret, const MTL& M,const Plates& P, const PP& pp, const Feat& F, const Assignment& A) {
    std::vector<char> colors;
    colors.resize(F.Categories);
    colors[0] = 'k'; colors[1] = 'g'; colors[2] = 'r'; colors[3] = 'b'; colors[4] = 'm'; colors[5] = 'y'; colors[6] = 'w'; colors[7] = 'c';
    polygon pol;
    PosP posp(3,3);

    int j=A.suborder[jused];
    for (int k=0; k<F.Nfiber; k++) {
        dpair O = pp.coords(k);
        int g = A.TF[j][k];
        if (g!=-1) {
            dpair Ga = projection(g,j,M,P);
            
            polygon fh = F.fh;
            polygon cb = F.cb;
            repos_cb_fh(cb,fh,O,Ga,posp);
            //if (A.is_collision(j,k,pp,G,P,F)!=-1) {
            //cb.set_color('r');
            //fh.set_color('r');
            //}
            //fix color assignment to account for Secret and SkyF
            char this_color;
            if (g<F.Ntarg){
                this_color=colors[Secret[g].category];
            }
            else if (g<F.Ntarg+F.NSStars) this_color='w';
            else this_color='c';

            cb.set_color(this_color);
            fh.set_color(this_color);
            pol.add(cb);
            pol.add(fh);
            if(this_color!='w'){
            pol.add(element(O,this_color,0.3,5));
            }
            else             pol.add(element(O,'k',0.1,5));
        }
        else pol.add(element(O,'k',0.1,3));//unassigned fiber
        List av_gals = P[j].av_gals[k];

        for (int i=0; i<av_gals.size(); i++) {
            int gg = av_gals[i];
            char this_color;
            if (gg>F.Ntarg+F.NSStars){
                this_color='c';
                dpair Ga = projection(gg,j,M,P);
                pol.add(element(Ga,this_color,1,0.5));
             }
            else if (gg>F.Ntarg){
                this_color='w';
                dpair Ga = projection(gg,j,M,P);
                pol.add(element(Ga,this_color,1,0.5));
            }
            else if(1<=A.nobs_time(gg,j,Secret,M,F)){
                this_color=colors[Secret[gg].category];
                dpair Ga = projection(gg,j,M,P);
                if (this_color=='k') pol.add(element(Ga,'k',1,A.is_assigned_jg(j,gg)==-1?0.9:0.5));
                else pol.add(element(Ga,this_color,1,0.5));
            }
        }
    }

    pyplot pyp(pol);

    //for (int k=0; k<F.Nfiber; k++) pyp.addtext(pp.coords(k),i2s(k)); // Plot fibers identifiers
    
    pyp.plot_tile(directory,j,F); 
        
    
}

void overlappingTiles(str fname, const Feat& F, const Assignment& A) {
  FILE * file;
  file = fopen(fname.c_str(),"w");
  for (int g=0; g<F.Ngal; g++) {
    if (A.GL[g].size()==5) {
      fprintf(file,"%d ",g);
      for (int i=0; i<A.GL[g].size(); i++) fprintf(file,"(%d,%d) ",A.GL[g][i].f,A.GL[g][i].s);
      fprintf(file,"\n");
    }
  }
  fclose(file);
}
