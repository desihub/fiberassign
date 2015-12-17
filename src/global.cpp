#include	<cstdlib>
#include	<cmath>
#include        <cstdio>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include	<string>
#include        <cstring>
#include	<vector>
#include	<algorithm>
#include	<exception>
#include        <stdexcept>
#include	<sys/time.h>
#include        <stdlib.h>     /* srand, rand */
#include	"modules/htmTree.h"
#include	"modules/kdTree.h"
#include        "omp.h"
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "global.h"
#include    <sys/stat.h>

extern "C" {
	#include <sys/stat.h>
}

#include <stdexcept>


// Collecting information from input -------------------------------------------------------------------------------------
void collect_galaxies_for_all(const MTL& M, const htmTree<struct target>& T, Plates& P, const PP& pp, const Feat& F) {
    //provides list of galaxies available to fiber k on tile j: P[j].av_gals[k]
	Time t;
	init_time(t,"# Begin collecting available galaxies");
	List permut = random_permut(F.Nplate);
	double rad = F.PlateRadius*M_PI/180.;
	//int jj;
	//omp_set_num_threads(24);
    #pragma omp parallel
	{ 	int id = omp_get_thread_num(); if (id==0) printf(" ");
		// Collects for each plate
        // start at jj=0 not id
        #pragma omp for
        for (int jj=0; jj<F.Nplate; jj++){
			int j = permut[jj];
			plate p = P[j];
			// Takes neighboring galaxies that fall on this plate
			std::vector<int> nbr = T.near(M,p.nhat,rad);
			// Projects thoses galaxies on the focal plane
			Onplates O;
			for (int gg=0; gg<nbr.size(); gg++) {
				int g = nbr[gg];

				struct onplate op = change_coords(M[g],p);

				op.id = g;
				O.push_back(op);
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
                        if(M[gals[g]].t_priority==9900){
                            P[j].SS_av_gal[q].push_back(gals[g]);
                        }
                        if(M[gals[g]].t_priority==9800){
                            P[j].SF_av_gal[q].push_back(gals[g]);
                        }

                    }
                }
            }
        }
    }

	print_time(t,"# ... took :");
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
    if(M[g].t_priority==9900 || M[g].t_priority==9900) return false;
    if (P[j].ipass==4 && M[g].lastpass==0){
        return false;} // Only ELG at the last pass
	if (F.Collision) for (int i=0; i<pp.N[k].size(); i++) if (g==A.TF[j][pp.N[k][i]]) return false; // Avoid 2 neighboring fibers observe the same galaxy (can happen only when Collision=true)
    if (A.find_collision(j,k,g,pp,M,P,F)!=-1){
        return false;} // No collision
    //if (M[g].t_priority==9800 && )
	return true;
    //doesn't require that jk is unassigned//doesn't require that g isn't assigned already on this plate
    //use is_assigned_jg for this
}

// makes sure we don't exceed limit on SS and SF
inline bool ok_for_limit_SS_SF(int g, int j, int k, const MTL& M, const Plates& P, const PP& pp, const Feat& F){
    bool is_SF=M[g].t_priority==9800;
    bool too_many_SF=P[j].SF_in_petal[pp.spectrom[k]]>F.MaxSF-1;
    bool is_SS=M[g].t_priority==9900;
    bool too_many_SS=P[j].SS_in_petal[pp.spectrom[k]]>F.MaxSS-1;
    return !(is_SF&&too_many_SF)&&!(is_SS&&too_many_SS);
}

    

// Find, for (j,k), find the best galaxy it can reach among the possible ones
// Null list means you can take all possible kinds, otherwise you can only take, for the galaxy, a kind among this list
// Not allowed to take the galaxy of id no_g
inline int find_best(int j, int k, const MTL& M, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, int no_g=-1, List kind=Null()) {
	int best = -1; int mbest = -1; int pbest = 10000;
	List av_gals = P[j].av_gals[k];
	// For all available galaxies
	for (int gg=0; gg<av_gals.size(); gg++) {
		int g = av_gals[gg];
        if(ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
            int m = M[g].nobs_remain; // Check whether it needs further observation
            if (m>=1) {
                int prio = M[g].t_priority;
            
                // Takes it if better priority, or if same, if it needs more observations, so shares observations if two QSOs are close
                if (prio<pbest || (prio==pbest && m>mbest)) {

                    // Check that g is not assigned yet on this plate, or on the InterPlate around, check with ok_to_assign
                    int isa=A.is_assigned_jg(j,g,M,F);
                    int ok=ok_assign_g_to_jk(g,j,k,P,M,pp,F,A);
                    if (isa==-1 && ok && g!=no_g ) {
                        best = g;
                        pbest = prio;
                        mbest = m;
                    }
                }
            }
        }
    }
	return best;
}

// Tries to assign the fiber (j,k)
inline int assign_fiber(int j, int k, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1, List kind=Null()) {
	if (A.is_assigned_tf(j,k)) return -1;
	int best = find_best(j,k,M,P,pp,F,A,no_g,kind);
    int g=best;

    if (best!=-1) A.assign(j,k,best,M,P,pp);
	return best;
}


// Tries to assign the galaxy g to one of the plates list starting with the jstart one, and of size size
//default jstart is as A.next_plate
//default size is number of plates to go
//used only in replace
inline void assign_galaxy(int g,  MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int jstart=-1, int size=-1) {
	int j0 = (jstart==-1) ? A.next_plate : jstart;
	int n = (size==-1) ? F.Nplate-j0 : size;// number of plates to do
	int jb = -1; int kb = -1; int unusedb = -1;
	Plist av_tfs = M[g].av_tfs;
	// All the tile-fibers that can observe galaxy g
	for (int tfs=0; tfs<av_tfs.size(); tfs++) {
		int j = av_tfs[tfs].f;
		int k = av_tfs[tfs].s;
		// Check if the assignment is possible, if ok, if the tf is not used yet, and if the plate is in the list
		if (j0<j && j<j0+n && !A.is_assigned_tf(j,k) && ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)) {
			int unused = A.unused[j][pp.spectrom[k]];//unused fibers on this petal
			if (unusedb<unused) {
				jb = j; kb = k; unusedb = unused;//observe this galaxy by fiber on petal with most free fibefs
			}
		}
	}
	if (jb!=-1) A.assign(jb,kb,g,M,P,pp);
}

// Takes an unassigned fiber and tries to assign it with the "improve" technique described in the doc
// not used for SS or SF   
inline int improve_fiber(int begin, int next, int j, int k, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1) {
	if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
		// tries to assign it in the conventional way to galaxy available to it
		int g_try = assign_fiber(j,k,M,P,pp,F,A,no_g);
		if (g_try!=-1) return g_try;
		else { // Improve
			int gb = -1; int bb = -1; int jpb = -1; int kpb = -1; int mb = -1; int pb = 1e3; int unusedb = -1;
			List av_g = P[j].av_gals[k];
			// For all available galaxies within reach that are already observed
			for (int i=0; i<av_g.size(); i++) {
				int g = av_g[i];
				if (g!=-1 && g!=no_g) {
					if (ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)) {//this doesn't check to see that jk isnt assigned: it is
						// Which tile-fibers have taken g ?
						Plist tfs = A.chosen_tfs(g,F,begin,next);//all tile-fibers that observe g in tiles from begin to next
						for (int p=0; p<tfs.size(); p++) {
							int jp = tfs[p].f;
							int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
							// FIND BEST JP KP !!!
							int best = find_best(jp,kp,M,P,pp,F,A); // best!=g because !A.assigned_pg(best)

							if (best!=-1 && (A.is_assigned_jg(j,g,M,F)==-1 || jp==j)) {
								int prio = M[g].t_priority;
								int m = M[g].nobs_remain;
								int unused = A.unused[jp][pp.spectrom[kp]]; // We take the most unused
								if (prio<pb || (prio==pb && m>mb) || (prio==pb && m==mb && unused>unusedb)) {
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
//not used !
// Assignment functions ------------------------------------------------------------------------------------------
// Assign fibers naively
// Not used at present

void simple_assign(MTL &M, Plates& P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin simple assignment :");
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-j0 : next; // Not F.Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	//List randPlates = F.Randomize ? random_permut(plates) : plates;
    
	for (int jj=0; jj<n; jj++) {
		//int j = randPlates[jj];
        int j=jj;
		List randFibers = random_permut(F.Nfiber);
		for (int kk=0; kk<F.Nfiber; kk++) { // Fiber
			int k = randFibers[kk];
			assign_fiber(j,k,M,P,pp,F,A);
		}
	}
	str next_str = next==-1 ? "all left" : f(n);
	if (next!=1) print_time(t,"# ... took :");
}

void improve( MTL& M, Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin improve :");
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-j0 : next;
	int na_start = A.na(F,j0,n);//number of assigned tile-fibers from j0 to jo+n-1
	List plates = sublist(j0,n,A.order);
	//List randPlates = F.Randomize ? random_permut(plates) : plates;
	//for (int jj=0; jj<n; jj++) for (int k=0; k<F.Nfiber; k++) improve_fiber(j0,n,randPlates[jj],k,M,P,pp,F,A);
    for (int jj=0; jj<n; jj++) for (int k=0; k<F.Nfiber; k++) improve_fiber(j0,n,plates[jj],k,M,P,pp,F,A);
	int na_end = A.na(F,j0,n);
	printf("  %s more assignments (%.3f %% improvement)\n",f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));//how many new assigned tf's
	if (next!=1) print_time(t,"# ... took :");
}

// If there are galaxies discovered as fake for example, they won't be observed several times in the plan
// haas access to G,not just M, because it needs to know the truth



void new_replace( int j, int p, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    // do standard stars,going through priority classes from least to most
    // skip SS and SF, so start at size -3
    //can get all available SS,SF on plate from P[j].av_gals_plate restricting to plate p
    for(int c=M.priority_list.size()-3;P[j].SS_in_petal[p]<F.MaxSS && c>-1;--c ){//try to do this for lowest priority
        // aside from SS and SF, so size()-3
        std::vector <int> gals=P[j].SS_av_gal[p]; //standard stars on this petal
        for(int gg=0;gg<gals.size() ;++gg){
            int g=gals[gg];//a standard star
            if(A.is_assigned_jg(j,g)==-1){
            Plist tfs=M[g].av_tfs;//all tiles and fibers that reach g
            int done=0;//quit after we've used this SS
            for(int i=0;i<tfs.size() && done==0;++i){
                if(tfs[i].f==j&&pp.spectrom[tfs[i].s]==p){//a combination on this plate
                    int k=tfs[i].s;//we know g can be reached by this petal of plate j and fiber k
                    int g_old=A.TF[j][k];//what is now at (j,k)  g_old can't be -1 or we would have used it already in assign_sf
                    if (M[g_old].priority_class==c&&A.is_assigned_jg(j,g,M,F)==-1&& ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                        //right priority; this SS not already assigned on this plate
                        A.unassign(j,k,g_old,M,P,pp);
                        assign_galaxy(g_old,M,P,pp,F,A);//try to assign
                        A.assign(j,k,g,M,P,pp);
                        done=1;
                    }
                }
            }
            }
        }
    }
    for(int c=M.priority_list.size()-3;P[j].SF_in_petal[p]<F.MaxSF && c>-1;--c ){//try to do this for lowest priority
        // aside from SS and SF, so size()-3
        std::vector <int> gals=P[j].SF_av_gal[p]; //standard stars on this plate
        for(int gg=0;gg<gals.size();++gg){//what tfs for this SS?  M[g].av_tfs
            int g=gals[gg];//a sky fiber
            if(A.is_assigned_jg(j,g)==-1){
            Plist tfs=M[g].av_tfs;
            int done=0;
            for(int i=0;i<tfs.size() && done==0;++i){
                if(tfs[i].f==j&&pp.spectrom[tfs[i].s]==p){
                    int k=tfs[i].s;//we know g can be reached by this petal of plate j and fiber k
                    int g_old=A.TF[j][k];//what is now at (j,k)
                   if (M[g_old].priority_class==c&&A.is_assigned_jg(j,g,M,F)==-1 && ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                        A.unassign(j,k,g_old,M,P,pp);
                        assign_galaxy(g_old,M,P,pp,F,A);//try to assign
                        A.assign(j,k,g,M,P,pp);
                        done=1;
                    }
                }
            }
            }
        }
    }
   

}


void assign_unused(int j, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
    // Tries to assign remaining fibers in tile j
    //even taking objects observed later
	for (int k=0; k<F.Nfiber; k++) {
		if (!A.is_assigned_tf(j,k)) {
			int best = -1; int mbest = -1; int pbest = 100000; int jpb = -1; int kpb = -1;
			List av_gals = P[j].av_gals[k];//all available galaxies for this fiber k
			for (int gg=0; gg<av_gals.size(); gg++) {
				int g = av_gals[gg];//available galaxies
				int m = M[g].nobs_remain;
				int prio = M[g].t_priority;
				if (prio<pbest || (prio==pbest && m>mbest)) {
                    if (A.is_assigned_jg(j,g,M,F)==-1 && ok_assign_g_to_jk(g,j,k,P,M,pp,F,A)&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                        //not assigned this plate or within excluded interval
						for (int i=0; i<A.GL[g].size(); i++) { //GL[g].size() is number of tf that could observe g
							int jp = A.GL[g][i].f;
							int kp = A.GL[g][i].s;
							if (j<jp && jpb<jp) {//take best opportunity
								best = g;
								pbest = prio;
								mbest = m;
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

// If not enough SS and SF,

void assign_sf_ss(int j, MTL& M, Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	List randPetals = random_permut(F.Npetal);
	for (int ppet=0; ppet<F.Npetal; ppet++) {
		int p = randPetals[ppet];
		List randFibers = random_permut(pp.fibers_of_sp[p]);//fibers for this petal
        //first use any free fibers
			for (int kk=0; kk<F.Nfbp; kk++) {
				int k = randFibers[kk];
                if (A.TF[j][k]==-1){
                    //look at available galaxies for (j.k)
                    int done=0;
                    List av_gals = P[j].av_gals[k];
                    for (int gg=0; gg<av_gals.size()&&done==0; gg++) {
                        int g = av_gals[gg];//galaxy at (j,k)
                        if(M[g].t_priority==9900&&A.is_assigned_jg(j,g,M,F)==-1&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                            A.assign(j,k,g,M,P,pp);
                            done=1;
                        }
                        else{
                            if(M[g].t_priority==9800&&A.is_assigned_jg(j,g,M,F)==-1&&ok_for_limit_SS_SF(g,j,k,M,P,pp,F)){
                                A.assign(j,k,g,M,P,pp);
                                done=1;
                            }
                        }
                    }
                }
            }
			// If not enough SS and SF, replace galaxies with lowest priority
        new_replace(j,p,M,P,pp,F,A);
    }
}


void redistribute_tf(MTL& M, Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin redistribute TF :");
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-A.next_plate : next; //from next_plate on
	List plates = sublist(j0,n,A.order);
	//List randPlates = F.Randomize ? random_permut(plates) : plates;
	List randPlates = random_permut(plates);
	int red(0);
	Table Done = initTable(F.Nplate,F.Nfiber);//consider every plate and every fiber
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randFiber = random_permut(F.Nfiber);
		for (int kk=0; kk<F.Nfiber; kk++) {
			int k = randFiber[kk];
			if (Done[j][k]==0) {
				int g = A.TF[j][k];//current assignment of (j,k)  only look if assigned
				if (g!=-1) {
					int jpb = -1; int kpb = -1; int unusedb = A.unused[j][pp.spectrom[k]];//unused for j, spectrom[k]
					Plist av_tfs = M[g].av_tfs;  //all possible tile fibers for this galaxy
					for (int i=0; i<av_tfs.size(); i++) {
						int jp = av_tfs[i].f;
						int kp = av_tfs[i].s;
						int unused = A.unused[jp][pp.spectrom[kp]];//unused for jp, spectrom[kp]
						if (j0<=jp && jp<j0+n && !A.is_assigned_tf(jp,kp) && Done[jp][kp]==0 && ok_assign_g_to_jk(g,jp,kp,P,M,pp,F,A) && A.is_assigned_jg(jp,g,M,F)==-1 && 0<unused) {
							if (unusedb<unused) { // Takes the most usused petal
                                jpb = jp;
								kpb = kp;
								unusedb = unused;
							}
						}
					}
					if (jpb!=-1&&ok_for_limit_SS_SF(g,jpb,kpb,M,P,pp,F)) {
						A.unassign(j,k,g,M,P,pp);
						A.assign(jpb,kpb,g,M,P,pp);
						Done[j][k] = 1;
						Done[jpb][kpb] = 1;
						red++; 
					}
				}
			}
		}
	}
	printf("  %s redistributions of tile-fibers \n",f(red).c_str());
	if (next!=1) print_time(t,"# ... took :");
}


void write_FAtile_ascii(int j, str outdir, const MTL& M, const Plates& P, const PP& pp, const Feat& F, const Assignment& A) {
	FILE * FA;
	str s = outdir+"tile"+i2s(j)+".txt";
	FA = fopen(s.c_str(),"w");
	for (int k=0; k<F.Nfiber; k++) {
		int g = A.TF[j][k];
		// k
		fprintf(FA,"%d ",k);
		List av_gals = P[j].av_gals[k];
		// Number of potential galaxies
		fprintf(FA,"%lu ",av_gals.size());
		// IDs of potential galaxies
		for (int i=0; i<av_gals.size(); i++) fprintf(FA,"%ld ",M[av_gals[i]].id);
		// galaxy number, ra, dec, x, y
		if (g!=-1) {
			dpair Gal = projection(g,j,M,P);
            fprintf(FA,"%ld %f %f %f %f\n",M[g].id,M[g].ra,M[g].dec,Gal.f,Gal.s);
		}
		else fprintf(FA,"-1\n");
	}
	fclose(FA);
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
    
    // check if the file exists, and if so, throw an exception
    
    char filename[cfilesize];
    // int ret = snprintf(filename, cfilesize, "%s/tile_%05d.fits", outdir.c_str(), j);
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
    
    size_t ncols = 10;
    
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
    strcpy(tform[2], "J");
    strcpy(tunit[2], "");
 
    strcpy(ttype[3], "PRIORITY");
    strcpy(tform[3], "J");
    strcpy(tunit[3], "");

    //strcpy(ttype[3], "objtype");
    //snprintf(tform[3], FLEN_VALUE, "%dA", (int)objtypelen);
    //strcpy(tunit[3], "");
    
    strcpy(ttype[4], "TARGETID");
    strcpy(tform[4], "K");
    strcpy(tunit[4], "");
    
    strcpy(ttype[5], "TARGETFLAG");
    strcpy(tform[5], "K");
    strcpy(tunit[5], "");
    
    strcpy(ttype[6], "RA");
    strcpy(tform[6], "E");
    strcpy(tunit[6], "deg");
    
    strcpy(ttype[7], "DEC");
    strcpy(tform[7], "E");
    strcpy(tunit[7], "deg");
    
    strcpy(ttype[8], "XFOCAL_DESIGN");
    strcpy(tform[8], "E");
    strcpy(tunit[8], "mm");
    
    strcpy(ttype[9], "YFOCAL_DESIGN");
    strcpy(tform[9], "E");
    strcpy(tunit[9], "mm");
    
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
    
    int fiber_id[optimal];
    int positioner_id[optimal];
    int num_target[optimal];
    char objtype[optimal][objtypelen];
    char * ot_tmp[optimal];
    for (int i = 0; i < optimal; i++) {
        ot_tmp[i] = objtype[i];
    }
    long long target_id[optimal];
    long long desi_target[optimal];
    float ra[optimal];
    float dec[optimal];
    float x_focal[optimal];
    float y_focal[optimal];
    //new
    float t_priority[optimal];
    
    std::vector <long long> potentialtargetid;
    
    // write data in buffered way
    
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

                desi_target[i] = 0;
		//PRUEBA
                //target_id[i] = g;
		target_id[i] = M[g].id;

                if (g < 0) {
                    //strcpy(objtype[i], "NA");
                    ra[i] = qNan;
                    dec[i] = qNan;
                    x_focal[i] = qNan;
                    y_focal[i] = qNan;
                } else {
                    //we aren't supposed to know the kind  use priority instead
                    //strncpy(objtype[i], F.kind[G[g].id].c_str(), objtypelen);
                    ra[i] = M[g].ra;
                    dec[i] = M[g].dec;
                    dpair proj = projection(g,j,M,P);
                    x_focal[i] = proj.f;
                    y_focal[i] = proj.s;
                    t_priority[i]=M[g].t_priority;//new
                }
                
                for (int k = 0; k < P[j].av_gals[fib].size(); ++k) {
                    potentialtargetid.push_back(P[j].av_gals[fib][k]);
                }
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
            
            fits_write_col(fptr, TFLOAT, 7, offset+1, 1, n, ra, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 8, offset+1, 1, n, dec, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 9, offset+1, 1, n, x_focal, &status);
            fits_report_error(stderr, status);
            
            fits_write_col(fptr, TFLOAT, 10, offset+1, 1, n, y_focal, &status);
            fits_report_error(stderr, status);
        }
        
        offset += n;
    }
    
    // PotentialFiberMap table.  We have only one column, so it is safe
    // from a performance perspective to write the whole thing.
    
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
