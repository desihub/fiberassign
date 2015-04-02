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
#include        <stdlib.h>     /* srand, rand */
#include	"modules/htmTree.h"
#include	"modules/kdTree.h"
#include        "omp.h"
#include        "macros.h"
#include        "misc.h"
#include        "structs.h"
#include        "global.h"

// Useful sub-functions -------------------------------------------------------------------------------------------------
// Returns the radial distance on the plate (mm) given the angle,
// theta (radians).  This is simply a fit to the data provided.
inline double plate_dist(const double theta) {
	const double p[4] = {8.297e5,-1750.,1.394e4,0.0};
	double rr=0;
	for (int i=0; i<4; i++) rr = theta*rr + p[i];
	return(rr);
}

// Returns the x-y position on the plate centered at P for galaxy O.
inline struct onplate change_coords(const struct galaxy& O, const struct plate& P) {
	struct onplate obj;
	// Rotate the "galaxy" vector so that the plate center is at z-hat.
	double nhat1[3],nhat2[3];
	const double ct=P.nhat[2],st=sqrt(1-P.nhat[2]*P.nhat[2])+1e-30;
	const double cp=P.nhat[0]/st,sp=P.nhat[1]/st;
	// First rotate by -(Phi-Pi/2) about z. Note sin(Phi-Pi/2)=-cos(Phi)
	// and cos(Phi-Pi/2)=sin(Phi).
	nhat1[0] =  O.nhat[0]*sp - O.nhat[1]*cp;
	nhat1[1] =  O.nhat[0]*cp + O.nhat[1]*sp;
	nhat1[2] =  O.nhat[2];
	// then rotate by Theta about x
	nhat2[0] =  nhat1[0];
	nhat2[1] =  nhat1[1]*ct - nhat1[2]*st;
	nhat2[2] =  nhat1[1]*st + nhat1[2]*ct;
	// now work out the "radius" on the plate
	double tht=sqrt(sq(nhat2[0],nhat2[1]));
	double rad=plate_dist(tht);
	// the x-y position is given by our nhat's scaled by this
	obj.pos[0] = nhat2[0]/tht * rad;
	obj.pos[1] = nhat2[1]/tht * rad;
	return(obj);
}

// Collecting information from input -------------------------------------------------------------------------------------
// Fast because no big function is called ! Keep it in only 1 function
// For each ~10,000 plate, finds ~25,000 galaxies reachable by the plate,
// Projects them, and for each fiber, finds reachable ones
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp) {
	Time t;
	init_time(t,"# Begin collecting available galaxies");
	List permut = random_permut(Nplate);
	double rad = PlateRadius*M_PI/180.;
	int i;
	//omp_set_num_threads(24);
#pragma omp parallel
	{ 	int id = omp_get_thread_num();
		double quarter = floor(Nplate/(4*omp_get_num_threads()));
		double cnt, avg, std; int max = 0; int min = 1e7; 
		// Collects for each plate ; shuffle order of plates (useless !?)
		for (i=id; i<Nplate; i++) { // <- begins at id, otherwise all begin at 0 -> conflict. Do all plates anyway
			int j = permut[i];
			plate p = P[j];
			// Takes neighboring ~25000 galaxies that can be reached by this plate
			std::vector<int> nbr = T.near(G,p.nhat,rad); // nbr for neighbours
			// Projects thoses galaxies on the focal plane
			Onplates O;
			for (int k=0; k<nbr.size(); k++) {
				int g = nbr[k];
				struct onplate op = change_coords(G[g],p); 
				op.id = g;
				O.push_back(op);
			}
			// For each fiber, finds all reachable galaxies
			KDtree<struct onplate> kdT(O,2);
			for (int k=0; k<Nfiber; k++) {
				std::vector<int> gals = kdT.near(&(pp.fp[2*k]),0.0,PatrolRad);
				P[j].av_gals[k] = initList(gals);
			}
			// Avancement 
			//if (cnt==quarter && id==0) {printf("  Thread 0 has done 1/4 of his job"); print_time(t," at");}
		}
		//printf("  Thread %2d finished",id); print_time(t," at");
	} // End parallel
	print_time(t,"# ... took :");
}

void collect_available_tilefibers(Gals& G, const Plates& P) {
	Time t;
	init_time(t,"# Begin computing available tilefibers");
	for(int j=0; j<Nplate; j++) {
		for(int k=0; k<Nfiber; k++) {
			for(int m=0; m<P[j].av_gals[k].size(); m++) {
				int i = P[j].av_gals[k][m];
				G[i].av_tfs.push_back(pair(j,k));
			}
		}
	}
	print_time(t,"# ... took :");
}

// (On plate p) finds if there is a collision if fiber k would watch at galaxy g (collision with neighb)
inline int find_collision(int j, int k, int g, const PP& pp, const Gals& G, const Plates& P, const Assignment& A) {
	struct onplate op = change_coords(G[g],P[j]);
	double x = op.pos[0];
	double y = op.pos[1];
	for (int i=0; i<pp.N[k].size(); i++) {
		int kn = pp.N[k][i];
		int gn = A.TF[j][kn];
		if (gn!=-1) {
			struct onplate opn = change_coords(G[gn],P[j]);
			double xn = opn.pos[0];
			double yn = opn.pos[1];
			if (sq(x-xn,y-yn) < sq(Collide)) return kn;
		}
	}
	return -1;
}

inline bool ok_assign_g_to_jk(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Feat& F, const Assignment& A) {
	str kind = G[g].kind(F);
	if (P[j].ipass==Npass-1 && !(kind=="ELG" || kind=="SS" || kind=="SF")) return false; // Only ELG, SF, SS at the last pass
	if (find_collision(j,k,g,pp,G,P,A)!=-1) return false;
	if (kind=="SS" && A.nkind(j,k,"SS",G,P,pp,F)>=MaxSS) return false;
	if (kind=="SF" && A.nkind(j,k,"SF",G,P,pp,F)>=MaxSF) return false;
	return true;
}

inline bool ok_assign_tot(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Feat& F, const Assignment& A) {
	// No collision, only ELG in last pass
	if (!ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) return false;
	if (A.is_assigned_tf(j,k)) return false;
	if (A.is_assigned_jg(j,g,InterPlate)!=-1) return false;
	if (A.nobs(g,G,F)==0) return false;
	return true;
}

// Assignment sub-functions -------------------------------------------------------------------------------------
// There is not problem with the fact that we only have knowledge on previous assignment, because when we call nobs (the only moment it can raise a problem) there can't be crossings with this plate, because g can only be assigned once
// Be very careful of (j,k) calling this function in other ! (if improvement function doesn't work anymore it's likely that bad argu,ent where called
inline int find_best(int j, int k, const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, bool as_tf=false, int no_g=-1, List no_kind=Null()) {
	int best = -1; int mbest = -1; float pbest = 1e30;
	List av_gals = P[j].av_gals[k];
	for (int gg=0; gg<av_gals.size(); gg++) {
		int g = av_gals[gg];
		int prio = G[g].prio(F);
		int m = A.nobs(g,G,F);
		bool tfb = as_tf ? !A.is_assigned_tf(j,k) : true;
		if (tfb && ok_assign_g_to_jk(g,j,k,P,G,pp,F,A) && A.is_assigned_jg(j,g,InterPlate)==-1 && A.nobs(g,G,F)>=1 && g!=no_g && !isfound(G[g].id,no_kind)) {
			if (prio<pbest || A.once_obs[g] && m>mbest) { // Then g!=-1 because prio sup pbest
				best = g;
				pbest = prio;
				mbest = m;
			}
		}
	}
	return best;
}

inline int assign_fiber(int j, int k, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1, List no_kind=Null()) {
	int best = find_best(j,k,G,P,pp,F,A,true);
	if (best!=-1) A.assign(j,k,best,G,P,pp);
	return best;
}

// Try to use unused fibers by reassigning some used ones
// We could improve this function by taking best amongst all possible improvements !
inline int improve_fiber(int begin, int next, int j, int k, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1) {
	if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
		int g_try = assign_fiber(j,k,G,P,pp,F,A,no_g); // Doesn't add anything in practice, if there was an assignment first, though it can be useful in other cases
		if (g_try!=-1) return g_try;
		else { // Improve
			List av_g = P[j].av_gals[k];
			for (int i=0; i<av_g.size(); i++) { // Not shuffled
				int g = av_g[i]; // g : possible galaxy for (j,k)
				if (g!=-1 && g!=no_g) {
					// Is it allowed for jk to take g ?
					if (ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) {
						// What tfs have taken g ? Could they take someone else ?
						Plist tfs = A.chosen_tfs(g,begin,next);
						for (int p=0; p<tfs.size(); p++) {
							int jp = tfs[p].f;
							int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
							// FIND BEST JP KP !!!
							int best = find_best(jp,kp,G,P,pp,F,A,false,-1,Null()); // best!=g because !A.assigned_pg(best)
							//int best = best2(j,k,G,P,pp,F,A); // best!=g because !A.assigned_pg(best)
							if (best!=-1 && (jp==j || A.is_assigned_jg(j,g,InterPlate)==-1)) { // If j and jp are not on the same ip, g has to be unassgned

								// Modify assignment
								A.unassign(jp,kp,g,G,P,pp);
								A.assign(j,k,g,G,P,pp);
								A.assign(jp,kp,best,G,P,pp);
								return g;
							}}}}}}}
	return -1;
}

int improve_fiber_from_kind(int id, int j, int k, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
		int p = pp.spectrom[k];
		List fibskind = A.fibs_of_kind(id,j,p,G,pp,F);
		List no_kind; no_kind.push_back(id);
		bool finished(false);
		List av_gals = P[j].av_gals[k];
		for (int i=0; i<av_gals.size() && !finished; i++) {
			int g = av_gals[i];
			if (G[g].id==id && find_collision(j,k,g,pp,G,P,A)==-1 && A.is_assigned_jg(j,g,InterPlate)==-1 && A.nobs(g,G,F)>=1) {
				for (int kkp=0; kkp<fibskind.size() && !finished; kkp++) {
					int kp = fibskind[kkp];
					int gp = A.TF[j][kp];
					int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
					if (best!=-1) {
						A.unassign(j,kp,gp,G,P,pp);
						A.assign(j,k,g,G,P,pp);
						A.assign(j,kp,best,G,P,pp);
						finished = true;
						return g;
					}
				}
			}
		}
	}
	return -1;
}

// Assignment functions ------------------------------------------------------------------------------------------
// Assign fibers naively
void simple_assign(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin simple assignment :");
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-j0 : next; // Not Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randFibers = random_permut(Nfiber);
		for (int kk=0; kk<Nfiber; kk++) { // Fiber
			int k = randFibers[kk];
			assign_fiber(j,k,G,P,pp,F,A);
		}
	}
	str next_str = next==-1 ? "all left" : f(n);
	printf("  %s assignments on %s next plates\n",f(A.na(j0,n)).c_str(),next_str.c_str());
	if (next!=1) print_time(t,"# ... took :");
}

void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin improve :");
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-j0 : next;
	int na_start = A.na(j0,n);
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) for (int k=0; k<Nfiber; k++) improve_fiber(j0,n,randPlates[jj],k,G,P,pp,F,A);
	int na_end = A.na(j0,n);
	printf("  %s more assignments (%.3f %% improvement)\n",f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));
	if (next!=1) print_time(t,"# ... took :");
}

void improve_from_kind(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, str kind, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin improve "+kind+" :");
	int id = F.id(kind);
	List no_kind; no_kind.push_back(id);
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-A.next_plate : next;
	int na_start(A.na(j0,n));
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randPetals = random_permut(Npetal);
		for (int ppet=0; ppet<Npetal; ppet++) {
			int p = randPetals[ppet];
			// Take sublist of fibers assigned to kind, and unassigned ones
			List fibskind = A.fibs_of_kind(id,j,p,G,pp,F);
			List fibsunas = A.fibs_unassigned(j,p,G,pp,F);
			for (int kk=0; kk<fibsunas.size(); kk++) {
				// Take an unassigned tf and its available galaxies
				int k = fibsunas[kk];
				bool finished(false);
				List av_gals = P[j].av_gals[k];
				for (int i=0; i<av_gals.size() && !finished; i++) {
					int g = av_gals[i];
					if (G[g].id==id && find_collision(j,k,g,pp,G,P,A)==-1 && A.is_assigned_jg(j,g,InterPlate)==-1 && A.nobs(g,G,F)>=1) {
						// If the av gal is of the good kind, try to assign it, and improving an other one
						for (int kkp=0; kkp<fibskind.size() && !finished; kkp++) {
							int kp = fibskind[kkp];
							int gp = A.TF[j][kp];
							int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
							if (best!=-1) {
								A.unassign(j,kp,gp,G,P,pp);
								A.assign(j,k,g,G,P,pp);
								A.assign(j,kp,best,G,P,pp);
								finished = true;
								erase(kkp,fibskind);
	}}}}}}}
	int na_end(A.na(j0,n));
	printf("  %s more assignments (%.3f %% improvement)\n",f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));
	if (next!=1) print_time(t,"# ... took :");
}

// If there are galaxies discovered as fake for example, they won't be observed several times in the plan
void update_plan_from_one_obs(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int end) {
	int cnt(0);
	int j0 = A.next_plate;
	int n = end-j0+1;
	int na_start(A.na(j0,n));
	List to_update;
	// Declare that we've seen those galaxies
	A.update_once_obs(j0);
	// Get the list of galaxies to update in the plan
	for (int k=0; k<Nfiber; k++) {
		int g = A.TF[j0][k];
		// Only if once_obs, we delete all further assignment. obs!=obs_tmp means that the galaxy is a fake one for example (same priority but different goal)
		if (g!=-1 && A.nobsv_tmp[g]!=A.nobsv[g] && A.once_obs[g]) to_update.push_back(g);
	}
	// Update information on previously seen galaxies
	A.update_nobsv_tmp_for_one(j0);
	// Update further in the plan
	for (int gg=0; gg<to_update.size(); gg++) {
		int g = to_update[gg];
		Plist tfs = A.chosen_tfs(g,j0+1,n-1); // Begin at j0+1, because we can't change assignment at j0 (already watched)
		while (tfs.size()!=0) {
			int jp = tfs[0].f; int kp = tfs[0].s;
			//print_Plist("Before",tfs); // Debug
			A.unassign(jp,kp,g,G,P,pp);
			int ggg = -1;
			ggg = improve_fiber(j0+1,n-1,jp,kp,G,P,pp,F,A,g);
			if (ggg==-1) ggg = improve_fiber_from_kind(F.id("SF"),jp,kp,G,P,pp,F,A);
			if (ggg==-1) ggg = improve_fiber_from_kind(F.id("SS"),jp,kp,G,P,pp,F,A);
			erase(0,tfs);
			//print_Plist("After",tfs); // Debug
			cnt++;
		}
	}
	int na_end(A.na(j0,n));
	printf(" %4d unassignments and %4d replaced by update\n",cnt,na_end-na_start+cnt); fl();
}

// If no enough SS and SF, remove old_kind an replace to SS-SF (new_kind) on petal (j,p)
void replace(int old_kind, int new_kind, int j, int p, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	int m = A.nkind(j,p,F.kind[new_kind],G,P,pp,F,true);
	List fibskindd = A.fibs_of_kind(old_kind,j,p,G,pp,F);
	List fibskind = random_permut(fibskindd);
	int Max = new_kind==F.id("SS") ? MaxSS : MaxSF;
	while (m<Max && fibskind.size()!=0) {
		bool fin(false);
		int k = fibskind[0];
		List av_g = P[j].av_gals[k];
		for (int gg=0; gg<av_g.size() && !fin; gg++) {
			int g = av_g[gg];
			if (G[g].id==new_kind && find_collision(j,k,g,pp,G,P,A)==-1 && A.nobs(g,G,F)>=1 && pp.spectrom[k]==p) {
				int g0 = A.TF[j][k];
				A.unassign(j,k,g0,G,P,pp);
				A.assign(j,k,g,G,P,pp);
				fin = true;
				m++;
			}
		}
		erase(0,fibskind);
	}
}

// For each petal, assign QSOs, LRGs, ELGs, ignoring SS and SF. Then if there are free fibers, try to assign them first to SS and then SF. Now if we don't have 10 SS and 40 SF in a petal, take SS and SF at random from those that are available to the petal and if their fiber is assigned to an ELG, remove that assignment and give it instead to the SS or SF.
void new_assign_fibers(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin complexe assignment :");
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-j0 : next; // Not Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	List no_kind; no_kind.push_back(F.id("SS")); no_kind.push_back(F.id("SF")); 
	List no_reg; no_reg.push_back(F.id("QSO Ly-a")); no_reg.push_back(F.id("QSO Tracer"));  no_reg.push_back(F.id("LRG")); no_reg.push_back(F.id("ELG")); no_reg.push_back(F.id("Fake QSO")); no_reg.push_back(F.id("Fake LRG"));
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randPetals = random_permut(Npetal);
		for (int ppet=0; ppet<Npetal; ppet++) {
			// Assign without SS and SF
			int p = randPetals[ppet];
			List fibs = pp.fibers_of_sp[p];
			List randFibers = random_permut(fibs);
			for (int kk=0; kk<Nfbp; kk++) { // Fiber
				int k = randFibers[kk];
				assign_fiber(j,k,G,P,pp,F,A,-1,no_kind);
			}
			// Assign free fibers to SS-SF
			List fibsunas = A.fibs_unassigned(j,p,G,pp,F);
			for (int kk=0; kk<fibsunas.size(); kk++) { // Fiber
				int k = fibsunas[kk];
				assign_fiber(j,k,G,P,pp,F,A,-1,no_reg);
			}
			// If no enough SS and SF, remove ELG an replace to SS-SF
			replace(F.id("ELG"),F.id("SS"),j,p,G,P,pp,F,A);
			replace(F.id("ELG"),F.id("SF"),j,p,G,P,pp,F,A);
			replace(F.id("LRG"),F.id("SS"),j,p,G,P,pp,F,A);
			replace(F.id("LRG"),F.id("SF"),j,p,G,P,pp,F,A);
		}
	}
	str next_str = next==-1 ? "all left" : f(next);
	printf("  %s assignments on %s next plates\n",f(A.na(j0,next)).c_str(),next_str.c_str());
	if (next!=1) print_time(t,"# ... took :");
}

// Other useful functions --------------------------------------------------------------------------------------------
void results_on_inputs(const Gals& G, const Plates& P, const Feat& F, bool latex) {
	printf("# Results on inputs :\n");
	// Print features
	print_list("  Kinds corresponding :",F.kind);
	print_list("  Priorities :",F.prio);
	print_list("  Goals of observations :",F.goal);
	print_list("  Max goals of observations :",F.maxgoal());

	// How many galaxies in range of a fiber ?
	print_list("  How many galaxies in range of a fiber :",gals_range_fibers(P));

	// Histograms on number of av gals per plate and per fiber
	Cube T = initCube(Categories,Nplate,Nfiber);
	for (int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			List gals = P[j].av_gals[k];
			for (int g=0; g<gals.size(); g++) T[G[gals[g]].id][j][k]++;
		}
	}
	Table hist1;
	for (int id=0; id<Categories; id++) hist1.push_back(histogram(T[id],1));
	print_mult_hist_latex("  Available galaxies (by kind) for a TF",1,hist1,true);

	// Histograms on number of av tfs per galaxy
	Table Tg = initTable(Categories,0);
	for (int g=0; g<Ngal; g++) {
		int n = G[g].av_tfs.size();
		Tg[G[g].id].push_back(n); 
	}
	Table hist2;
	for (int id=0; id<Categories; id++) hist2.push_back(histogram(Tg[id],1));
	print_mult_hist_latex("   Available tile-fibers for a galaxy (by kind)",1,hist2,true);
}

// Detemine how many galaxies need to be dropped to guarantee 40 free fibers for each petal
// Count how many free fibers there are beyond 500 in each plate
void print_free_fibers(const Gals& G, const PP& pp, const Feat& F, const Assignment& A, bool latex) {
	printf("# Free fibers statistics\n");
	// Prints histogram of free fibers
	Table unused_fbp = A.unused_fbp(pp);
	Table hist; hist.push_back(histogram(unused_fbp,5));
	print_mult_hist_latex("  Number of petals with this many free fiber",5,hist,true);

	// Histogram of SS
	Table usedSS = A.used_by_kind("SS",G,pp,F);
	print_hist("  UsedSS (number of petals)",1,histogram(usedSS,1));

	// Histogram of SF
	Table usedSF = A.used_by_kind("SF",G,pp,F);
	print_hist("  UsedSF (number of petals)",1,histogram(usedSF,1));
}

void display_results(const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, bool latex) {
	printf("# Results :\n");
	// Raw numbers of galaxies by id and number of remaining observations
	Table hist2 = initTable(Categories,MaxObs+1);
	Table ob = initTable(Categories+1,2*MaxObs+1);
	List maxobs = F.maxgoal();
	for (int g=0; g<Ngal; g++) {
		int m = A.nobs(g,G,F,false);
		if (m>=0 && m<=MaxObs) hist2[G[g].id][m]++;
		ob[G[g].id][m+MaxObs]++;
	}
	for (int m=0; m<2*MaxObs+1; m++) ob[Categories][m] = m-MaxObs;
	print_table("  Remaining observations (without negative obs ones)",with_tot(hist2),latex);
	print_table("  Real remaining observations",with_tot(ob));

	// Percentages of observation
	Dtable percents = initDtable(Categories,2);
	for (int id=0; id<Categories; id++) {
		int tot = sumlist(ob[id]);
		int goal = F.goal[id];

		percents[id][0] = percent(tot-ob[id][MaxObs+goal],tot);

		double d(0.0);
		for (int i=0; i<MaxObs; i++) d += ob[id][MaxObs+i]*(goal-i);
		percents[id][1] = percent(d,tot*goal);
	}
	print_table("Percentages of once observed, and with ponderation",percents,latex);
	
	// Observed galaxies in function of time
	Table Ttime_scaled = initTable(Categories,Nplate);
	Table Ttime = initTable(Categories,Nplate);
	for (int id=0; id<Categories; id++) {
		int n(0);
		for (int j=0; j<Nplate; j++) {
			for (int p=0; p<Npetal; p++) n += A.nkind(j,p,F.kind[id],G,P,pp,F,true);
			Ttime[id][j] = n;
			Ttime_scaled[id][j] = n;
		}
		for (int j=0; j<Ttime[id].size(); j++) Ttime_scaled[id][j] /= n;
	}
	print_mult_hist_latex("Observed galaxies in function of time (scaled)",1,Ttime_scaled);
	print_mult_hist_latex("Observed galaxies in function of time",1,Ttime);

	// Histogram of time between 2 obs of Ly a
	Table deltas;
	for (int g=0; g<Ngal; g++) {
		if (G[g].id == F.id("QSO Ly-a")) {
			Plist tfs = A.chosen_tfs(g);
			if (tfs.size()>=2) {
				List del;
				for (int i=0; i<tfs.size()-1; i++) {
					int p1 = tfs[i].f;
					int p2 = tfs[i+1].f;
					del.push_back(p2-p1);
				}
				deltas.push_back(del);
			}
		}
	}
	deb(99999);
	List h1 = histogram(deltas,100);
	deb(555);
	print_mult_hist_latex("Plate interval between 2 consecutive obs of Ly-a (interval 100)",100,histogram(deltas,100),true);

	// Some petal
	//for (int i=1; i<=10; i++) {
		//print_table("Petal "+i2s(i),A.infos_petal(1000*i,5,G,P,pp,F));
	//}
	//plot_freefibers("free_fiber_plot.txt",P,A);   

	// Percentage of fiber assigned
	printf("  %s assignments in total (%.4f %% of all fibers)\n",f(A.na()).c_str(),percent(A.na(),Nplate*Nfiber));

	//// Some stats
	//Table done = initTable(Categories,MaxObs+1);
	//List fibers_used = initList(Categories);
	//List targets = initList(Categories);
	//for (int id=0; id<Categories; id++) {
		//List hist = hist2[id];
		//for (int i=0; i<=MaxObs; i++) {
			////i here is number of observations lacking  i = goal -done
			////done=goal -i for i<=goal
			////0 for done>goal, i.e. i<0
			//if (i<=F.goal[id]) {
				//done[id][i]=hist[F.goal[id]-i];
				//fibers_used[id]+=i*done[id][i];
				//targets[id] += done[id][i];
			//}
		//}
	//}

	//Dtable stats = initDtable(Categories,4);
	//for (int id=0; id<Categories; id++) {
		//stats[id][0] = (targets[id]-done[id][0])/TotalArea;
		//stats[id][1] = fibers_used[id]/TotalArea;
		//stats[id][2] = targets[id]/TotalArea;
		//stats[id][3] = percent(targets[id]-done[id][0],targets[id]);
	//}
	//print_table("  Id on lines. (Per square degrees) Total, Fibers Used, Available, Percent of observations",stats,latex);
}

void plot_freefibers(str s, const Plates& P, const Assignment& A) { // ? Not rechecked
	printf("# Plot free fibers positions\n");
	double pi = 3.1415926535;
	List unused_fibers = A.unused_f();
	FILE * fplot;
	const char * c = s.c_str();
	fplot = fopen(c,"w");
	fprintf(fplot,"j --- ra --- dec --- unused fibers \n");
	for (int j=0; j<Nplate; j++) {
		double phi = atan2(P[j].nhat[1],P[j].nhat[0]);
		double theta = acos(P[j].nhat[2]);
		double ra = 180.*phi/pi;
		double dec = 90.*(1.-2.*theta/pi);
		fprintf(fplot,"%d  %2d  %2d  %4d \n",j,ra,dec,unused_fibers[j]);
	}
	fclose(fplot);
}

// Wreck functions - could be useful later... or not ----------------------------------------------------------------------------
// Redistrubte assignments trying to get at least 500 free fibers on each plate/tile
// Redo so we look first at plates with too few free fibers
// Before : (j,k) <-> g, with Sp(k) too much used
// After : (jreassign,kreassign) <-> g & (j,k) free, such that (jreassign,kreassign) comes from most unused (ji,ki)
// Reassign if this is improvement
void redistribute(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	List Sp = pp.spectrom;
	Table unused_fbp = A.unused_fbp(pp); // We manipulate our own one to make computations quicker
	//print_table("unused",unused_fbp);
	List randPlates = random_permut(Nplate);
	// Consider all petals with too few free fibers
	for (int jj=0; jj<Nplate; jj++) {
		int j = randPlates[jj];
		List randFibers = random_permut(Nfiber);
		for (int kk=0; kk<Nfiber; kk++) {
			int k = randFibers[kk];
			int nfree = unused_fbp[j][Sp[k]];
			//int nfree = A.unused_fbp(j,k,pp);
			int g = A.TF[j][k];
			if (g!=-1 && nfree<MinUnused/* && G[g].id==4*/) { // Only ELG
				// Consider other ways to observe this galaxy
				std::vector<pair> tfs = G[g].av_tfs;
				int jreassign(-1); int kreassign(-1); int mostunused(-1);
				for (int c=0; c<tfs.size(); c++) {
					int jp = tfs[c].f;
					int kp = tfs[c].s;
					int nfreep = unused_fbp[jp][Sp[kp]];
					//int nfreep = A.unused_fbp(jp,kp,pp);
					// Use freest plate and petal
					if (nfreep>mostunused && ok_assign_g_to_jk(g,jp,kp,P,G,pp,F,A)) {
						mostunused = nfreep;
						jreassign = jp;
						kreassign = kp; 
					}
				}
				if (mostunused > nfree) {
					//printf("%d %d %d %d %d %d %d\n",mostunused, nfree, j,k,jreassign,kreassign,g);
					A.unassign(j,k,g,G,P,pp);
					A.assign(jreassign,kreassign,g,G,P,pp);
					unused_fbp[j][Sp[k]]++;
					unused_fbp[jreassign][Sp[kreassign]]--;		
					redistributions++; 
				}}}}
				printf("  %s redistributions (~%.4f %% redistributed)\n",f(redistributions).c_str(),percent(redistributions,A.na()));
}

void redistribute_tf(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	List randPlates = random_permut(Nplate);
	for (int jj=0; jj<Nplate; jj++) {
		printf("%d ",jj); std::cout.flush();
		int red(0);
		int j = randPlates[jj];
		for (int kk=0; kk<Nfiber; kk++) {
			List randFibers = random_permut(Nfiber);
			int k = randFibers[kk];
			int g = A.TF[j][k];
			if (g!=-1) {
				int best = -1; int mbest = -1; float pbest = 1e30;
				List av_gals = P[j].av_gals[k];
				List randGals = random_permut(av_gals.size());
				for (int gg=0; gg<av_gals.size(); gg++) {
					int g = av_gals[randGals[gg]];
					int prio = G[g].prio(F);
					int m = A.nobs(g,G,F);
					if (/*g!=g_no &&*/ ok_assign_g_to_jk(g,j,k,P,G,pp,F,A) && A.is_assigned_jg(j,g,InterPlate)==-1 && A.nobs(g,G,F)>=1) {
						if (prio<pbest || A.once_obs[g] && m>mbest) {
							best = g;
							pbest = prio;
							mbest = m;
						}
					}
				}
				if (best!=-1) {
					A.unassign(j,k,g,G,P,pp);
					A.assign(j,k,best,G,P,pp);
					redistributions++; red++;
				}
			}
		}
		printf("%d ",red);
	}
	printf("  %s redistributions (~%.4f %% of TF redistributed)\n",f(redistributions).c_str(),percent(redistributions,Nfiber*Nplate));
}

void redistribute_g(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	List randGals = random_permut(Ngal);
	for (int gg=0; gg<Ngal; gg++) {
		int g = randGals[gg];
		bool finished(false);
		std::vector<pair> av_tfs = G[g].av_tfs;
		std::vector<pair> tfs = A.chosen_tfs(g);
		for (int p=0; p<tfs.size() && !finished; p++){
			int j = tfs[p].f;
			int k = tfs[p].s; 
			for (int tf=0; tf<av_tfs.size() && !finished; tf++) {
				int jp = tfs[p].f;
				int kp = tfs[p].s; 
				if (find_collision(jp,kp,g,pp,G,P,A)==-1 && !A.is_assigned_tf(jp,kp) && P[j].ipass==P[jp].ipass) {
					A.unassign(j,k,g,G,P,pp);
					A.assign(jp,kp,g,G,P,pp);
					redistributions++;
					finished = true;
				}
			}
		}
	}
	printf("  %s redistributions (~%.4f %% of galaxies redistributed)\n",f(redistributions).c_str(),percent(redistributions,Ngal));
}

// Redistribute by kind
void redistribute_g_by_kind(str kind, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	int id = F.id(kind);
	List randGals = random_permut(Ngal);
	Table Done = initTable(Nplate,Nfiber); // Already redistributed TF
	for (int gg=0; gg<Ngal; gg++) {
		int g = randGals[gg];
		if (G[g].id==id) {
			bool finished(false);
			std::vector<pair> av_tfs = G[g].av_tfs;
			std::vector<pair> tfs = A.chosen_tfs(g);
			for (int p=0; p<tfs.size() && !finished; p++) {
				int j = tfs[p].f;
				int k = tfs[p].s;
				if (Done[j][k]==0) {
					for (int tf=0; tf<av_tfs.size() && !finished; tf++) {
						int jp = tfs[p].f;
						int kp = tfs[p].s; 
						// Ce qui bloque est assigned(jp,kp), déjà assigné !
						if (Done[jp][kp]==0 && find_collision(jp,kp,g,pp,G,P,A)==-1 && !A.is_assigned_tf(jp,kp) && P[j].ipass==P[jp].ipass) {
							A.unassign(j,k,g,G,P,pp);
							A.assign(jp,kp,g,G,P,pp);
							redistributions++;
							Done[j][k] = 1;
							Done[jp][kp] = 1;
							finished = true;
						}
					}
				}
			}
		}
	}
	printf("  %s redistributions (~%.4f %% of galaxies redistributed)\n",f(redistributions).c_str(),percent(redistributions,Ngal));
}
