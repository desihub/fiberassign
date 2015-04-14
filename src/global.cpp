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

inline bool ok_assign_g_to_jk(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Feat& F, const Assignment& A) {
	int kind = G[g].id;
	if (A.find_collision(j,k,g,pp,G,P)!=-1) return false;
	if (kind==F.ids.at("SF") && A.nkind(j,k,F.ids.at("SF"),G,P,pp,F)>=MaxSF) return false;
	if (kind==F.ids.at("SS") && A.nkind(j,k,F.ids.at("SS"),G,P,pp,F)>=MaxSS) return false;
	if (P[j].ipass==Npass-1 && !(kind==F.ids.at("ELG") || kind==F.ids.at("SS") || kind==F.ids.at("SF"))) return false; // Only ELG, SF, SS at the last pass
	return true;
}

// Assignment sub-functions -------------------------------------------------------------------------------------
// There is not problem with the fact that we only have knowledge on previous assignment, because when we call nobs (the only moment it can raise a problem) there can't be crossings with this plate, because g can only be assigned once
// Be very careful of (j,k) calling this function in other ! (if improvement function doesn't work anymore it's likely that bad argu,ent where called
inline int find_best(int j, int k, const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, bool as_tf=false, int no_g=-1, List no_kind=Null()) {
	int best = -1; int mbest = -1; int pbest = 1e3;
	List av_gals = P[j].av_gals[k];
	for (int gg=0; gg<av_gals.size(); gg++) {
		int g = av_gals[gg];
		int prio = G[g].prio(F);
		int m = A.nobs(g,G,F);
		bool tfb = as_tf ? !A.is_assigned_tf(j,k) : true;
		if (m>=1 && tfb && ok_assign_g_to_jk(g,j,k,P,G,pp,F,A) && A.is_assigned_jg(j,g,InterPlate)==-1 && g!=no_g && !isfound(G[g].id,no_kind)) {
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
inline int improve_fiber(int begin, int next, int j, int k, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A, int no_g=-1) {
	if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
		int g_try = assign_fiber(j,k,G,P,pp,F,A,no_g); // Doesn't add anything in practice, if there was an assignment first, though it can be useful in other cases
		if (g_try!=-1) return g_try;
		else { // Improve
			int gb = -1;
			int bb = -1;
			int jpb = -1;
			int kpb = -1;
			int mb = -1;
			int pb = 1e3;
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

							if (best!=-1 && (A.is_assigned_jg(j,g,InterPlate)==-1 || jp==j)) {
								int prio = G[best].prio(F);
								int m = A.nobs(best,G,F);
								if (prio<pb || A.once_obs[g] && m>mb) {
									gb = g;
									bb = best;
									jpb = jp;
									kpb = kp;
									mb = m;
									pb = prio;
							}}}}}}
			// Modify assignment
			if (gb!=-1) {
				A.unassign(jpb,kpb,gb,G,P,pp);
				A.assign(j,k,gb,G,P,pp);
				A.assign(jpb,kpb,bb,G,P,pp);
				return gb;
			}
		}
	}
	return -1;
}

int improve_fiber_from_kind(int id, int j, int k, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
		int p = pp.spectrom[k];
		List fibskind = A.fibs_of_kind(id,j,p,G,pp,F);
		List no_kind; no_kind.push_back(id);
		int gb = -1;
		int gpb = -1;
		int bb = -1;
		int kpb = -1;
		int mb = -1;
		int pb = 1e3;
		List av_gals = P[j].av_gals[k];
		for (int i=0; i<av_gals.size(); i++) {
			int g = av_gals[i];
			if (G[g].id==id && A.find_collision(j,k,g,pp,G,P)==-1 && A.is_assigned_jg(j,g,InterPlate)==-1 && A.nobs(g,G,F)>=1) {
				for (int kkp=0; kkp<fibskind.size(); kkp++) {
					int kp = fibskind[kkp];
					int gp = A.TF[j][kp];
					int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
					if (best!=-1) {
						int prio = G[best].prio(F);
						int m = A.nobs(best,G,F);
						if (prio<pb || A.once_obs[g] && m>mb) {
							gb = g;
							gpb = gp;
							bb = best;
							kpb = kp;
							mb = m;
							pb = prio;
		}}}}}
		// Modify assignment
		if (gb!=-1) {
			A.unassign(j,kpb,gpb,G,P,pp);
			A.assign(j,k,gb,G,P,pp);
			A.assign(j,kpb,bb,G,P,pp);
			return gb;
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
	int id = F.ids.at(kind);
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
				int gb = -1;
				int gpb = -1;
				int bb = -1;
				int kpb = -1;
				int kkpb = -1;
				int mb = -1;
				int pb = 1e3;
				List av_gals = P[j].av_gals[k];
				for (int i=0; i<av_gals.size(); i++) {
					int g = av_gals[i];
					if (G[g].id==id && A.nobs(g,G,F)>=1 && A.is_assigned_jg(j,g,InterPlate)==-1 && A.find_collision(j,k,g,pp,G,P)==-1) {
						// If the av gal is of the good kind, try to assign it, and improving an other one
						for (int kkp=0; kkp<fibskind.size(); kkp++) {
							int kp = fibskind[kkp];
							int gp = A.TF[j][kp];
							int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
							if (best!=-1) {
								int prio = G[best].prio(F);
								int m = A.nobs(best,G,F);
								if (prio<pb || A.once_obs[g] && m>mb) {
									gb = g;
									gpb = gp;
									bb = best;
									kpb = kp;
									kkpb = kkp;
									mb = m;
									pb = prio;
							}}}}}
			// Modify assignment
				if (gb!=-1) {
					A.unassign(j,kpb,gpb,G,P,pp);
					A.assign(j,k,gb,G,P,pp);
					A.assign(j,kpb,bb,G,P,pp);
					erase(kkpb,fibskind);
				}
			}
		}
	}
	int na_end(A.na(j0,n));
	if (next==1) printf(" %s-imp +%3s (+%.3f %%)",kind.c_str(),f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));
	else {
		printf("  %s more assignments (%.3f %% improvement)\n",f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));
		print_time(t,"# ... took :");
	}
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
			int gp = -1;
			gp = improve_fiber(j0+1,n-1,jp,kp,G,P,pp,F,A,g);
			if (gp==-1) gp = improve_fiber_from_kind(F.ids.at("SF"),jp,kp,G,P,pp,F,A);
			if (gp==-1) gp = improve_fiber_from_kind(F.ids.at("SS"),jp,kp,G,P,pp,F,A);
			erase(0,tfs);
			//print_Plist("After",tfs); // Debug
			cnt++;
		}
	}
	int na_end(A.na(j0,n));
	printf(" %4d unas & %4d replaced\n",cnt,na_end-na_start+cnt); fl();
}

// If no enough SS and SF, remove old_kind an replace to SS-SF (new_kind) on petal (j,p)
void replace(int old_kind, int new_kind, int j, int p, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	int m = A.nkind(j,p,new_kind,G,P,pp,F,true);
	List fibskindd = A.fibs_of_kind(old_kind,j,p,G,pp,F);
	List fibskind = random_permut(fibskindd);
	int Max = new_kind==F.ids.at("SS") ? MaxSS : MaxSF;
	while (m<Max && fibskind.size()!=0) {
		bool fin(false);
		int k = fibskind[0];
		List av_g = P[j].av_gals[k];
		for (int gg=0; gg<av_g.size() && !fin; gg++) {
			int g = av_g[gg];
			if (G[g].id==new_kind && A.find_collision(j,k,g,pp,G,P)==-1 && A.nobs(g,G,F)>=1 && pp.spectrom[k]==p) {
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
	if (next!=1) init_time(t,"# Begin new assignment :");
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-j0 : next; // Not Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	List no_kind; no_kind.push_back(F.ids.at("QSO Ly-a")); no_kind.push_back(F.ids.at("QSO Tracer")); no_kind.push_back(F.ids.at("Fake QSO")); no_kind.push_back(F.ids.at("SS")); no_kind.push_back(F.ids.at("SF")); 
	List no_reg; no_reg.push_back(F.ids.at("QSO Ly-a")); no_reg.push_back(F.ids.at("QSO Tracer"));  no_reg.push_back(F.ids.at("LRG")); no_reg.push_back(F.ids.at("ELG")); no_reg.push_back(F.ids.at("Fake QSO")); no_reg.push_back(F.ids.at("Fake LRG"));
	List only_qso; only_qso.push_back(F.ids.at("LRG")); only_qso.push_back(F.ids.at("ELG")); only_qso.push_back(F.ids.at("Fake LRG")); only_qso.push_back(F.ids.at("SS")); only_qso.push_back(F.ids.at("SF"));
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randPetals = random_permut(Npetal);
		for (int ppet=0; ppet<Npetal; ppet++) {
			// Assign only QSO
			int p = randPetals[ppet];
			List fibs = pp.fibers_of_sp[p];
			List randFibers = random_permut(fibs);
			for (int kk=0; kk<Nfbp; kk++) { // Fiber
				int k = randFibers[kk];
				assign_fiber(j,k,G,P,pp,F,A,-1,only_qso);
			}
			// Assign without SS and SF
			for (int kk=0; kk<Nfbp; kk++) { // Fiber
				int k = randFibers[kk];
				if (!A.is_assigned_tf(j,k)) assign_fiber(j,k,G,P,pp,F,A,-1,no_kind);
			}
			// Assign free fibers to SS-SF
			List fibsunas = A.fibs_unassigned(j,p,G,pp,F);
			for (int kk=0; kk<fibsunas.size(); kk++) { // Fiber
				int k = fibsunas[kk];
				assign_fiber(j,k,G,P,pp,F,A,-1,no_reg);
			}
			// If no enough SS and SF, remove ELG an replace to SS-SF
			replace(F.ids.at("ELG"),F.ids.at("SS"),j,p,G,P,pp,F,A);
			replace(F.ids.at("ELG"),F.ids.at("SF"),j,p,G,P,pp,F,A);
			replace(F.ids.at("LRG"),F.ids.at("SS"),j,p,G,P,pp,F,A);
			replace(F.ids.at("LRG"),F.ids.at("SF"),j,p,G,P,pp,F,A);
		}
	}
	str next_str = next==-1 ? "all left" : f(next);
	printf("  %s assignments on %s next plates\n",f(A.na(j0,n)).c_str(),next_str.c_str());
	if (next!=1) print_time(t,"# ... took :");
}

void redistribute_tf(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin redistribute TF :");
	int j0 = A.next_plate;
	int n = next==-1 ? Nplate-A.next_plate : next;
	List plates = sublist(j0,n,A.order);
	List randPlates = Randomize ? random_permut(plates) : plates;
	int red(0);
	Table Done = initTable(Nplate,Nfiber);
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randFiber = random_permut(Nfiber);
		for (int kk=0; kk<Nfiber; kk++) {
			int k = randFiber[kk];
			if (Done[j][k]==0) {
				int g = A.TF[j][k];
				if (g!=-1) {
					bool finished = false;
					Plist av_tfs = G[g].av_tfs;
					for (int i=0; i<av_tfs.size() && !finished; i++) {
						int jp = av_tfs[i].f;
						int kp = av_tfs[i].s;
						if (j0<=jp && jp<j0+n && !A.is_assigned_tf(jp,kp) && Done[jp][kp]==0 && ok_assign_g_to_jk(g,jp,kp,P,G,pp,F,A) && A.is_assigned_jg(jp,g,InterPlate)==-1) {
							A.unassign(j,k,g,G,P,pp);
							A.assign(jp,kp,g,G,P,pp);
							finished = true;
							Done[j][k] = 1;
							Done[jp][kp] = 1;
							red++; }}}}}}
	printf("  %s redistributions of couples of TF\n",f(red).c_str());
	if (next!=1) print_time(t,"# ... took :");
}

// Other useful functions --------------------------------------------------------------------------------------------
void results_on_inputs(str outdir, const Gals& G, const Plates& P, const Feat& F, bool latex) {
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
	print_mult_table_latex("Available galaxies (by kind) for a TF",outdir+"avgalhist.dat",hist1,1);

	// Histograms on number of av tfs per galaxy
	Table Tg = initTable(Categories,0);
	for (int g=0; g<Ngal; g++) {
		int n = G[g].av_tfs.size();
		Tg[G[g].id].push_back(n);
	}
	Table hist2;
	for (int id=0; id<Categories; id++) hist2.push_back(histogram(Tg[id],1));
	print_mult_table_latex("Available tile-fibers for a galaxy (by kind)",outdir+"avtfhist.dat",hist2,1);
}

void display_results(str outdir, const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, bool latex) {
	printf("# Results :\n");
	// 0 List of total number of galaxies
	List tots = initList(Categories);
	for (int g=0; g<Ngal; g++) tots[G[g].id]++;
	
	// 1 Raw numbers of galaxies by id and number of remaining observations
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

	// 2 Percentages of observation
	Dtable perc = initDtable(Categories,2);
	for (int id=0; id<Categories; id++) {
		int tot = sumlist(ob[id]);
		int goal = F.goal[id];

		perc[id][0] = percent(tot-ob[id][MaxObs+goal],tot);

		double d(0.0);
		for (int i=0; i<MaxObs; i++) d += ob[id][MaxObs+i]*(goal-i);
		perc[id][1] = percent(d,tot*goal);
	}
	print_table("Percentages of once observed, and with ponderation",perc,latex);

	// 3 Observed galaxies in function of time
	int interval = 10; // Interval of plates for graphic
	Dtable Ttime_scaled = initDtable(Categories,0);
	Table Ttime = initTable(Categories,0);
	for (int id=0; id<Categories; id++) {
		List l;
		int n = 0;
		for (int j=0; j<Nplate; j++) {
			for (int p=0; p<Npetal; p++) n += A.nkind(j,p,id,G,P,pp,F,true);
			if (j%interval==0) l.push_back(n);
		}
		Ttime[id] = l;
		Ttime_scaled[id] = percents(l,l[l.size()-1]);
	}
	print_mult_Dtable_latex("Observed galaxies in function of time (scaled) (interval 10)",outdir+"time.dat",Ttime_scaled,interval);
	//print_mult_table_latex("Observed galaxies in function of time",,F.kind,Ttime);
	
	// 4 Histogramme of percentages of seen Ly-a
	int id = F.ids.at("QSO Ly-a");
	int goal = F.goal[id];
	Table Percseen = initTable(goal+1,0);
	for (int g=0; g<Ngal; g++) {
		if (G[g].id==id) {
			int n = G[g].av_tfs.size();
			int p = A.chosen_tfs(g).size();
			if (n>=Percseen[p].size()) Percseen[p].resize(n+1);
			Percseen[p][n]++;
		}
	}
	make_square(Percseen);
	print_table("Number of QSO Ly-a : x - Number of available TF - y - Number of observations",Percseen);
	for (int j=0; j<Percseen[0].size(); j++) {
		for (int i=Percseen.size()-1; i!=0; i--) {
			Percseen[i-1][j] += Percseen[i][j];
		}
	}
	print_mult_table_latex("Available tile-fibers for a galaxy (by kind)",outdir+"obsly.dat",Percseen,1);

	// 5 Histogram of time between 2 obs of Ly a
	Table deltas;
	for (int g=0; g<Ngal; g++) {
		if (G[g].id == F.ids.at("QSO Ly-a")) {
			Plist tfs = A.chosen_tfs(g);
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
	print_hist("Plate interval between 2 consecutive obs of Ly-a (interval 100)",100,histogram(deltas,100));
	Table delts; delts.push_back(histogram(deltas,10));
	print_mult_table_latex("Plate interval between 2 consecutive obs of Ly-a (interval 10)",outdir+"dist2ly.dat",delts,10);

	// 6 Free fibers histogram
	Table unused_fbp = A.unused_fbp(pp);
	make_square(unused_fbp);
	Table hist0; hist0.push_back(histogram(unused_fbp,5));
	print_mult_table_latex("Number of petals with this many free fiber (interval 5)",outdir+"freefib.dat",hist0,5);

	// 7 Free fibers in function of time (plates)
	List freefibtime = initList(Nplate);
	for (int j=0; j<Nplate; j++) freefibtime[j] = A.unused_f(j);
	Table fft; fft.push_back(freefibtime);
	print_mult_table_latex("Free fibers in function of time (plates)",outdir+"fft.dat",fft);
 
	// 8 Percentage of seen objects as a function of density of objects
	


	// Misc
	// Histogram of SS
	Table usedSS = A.used_by_kind("SS",G,pp,F);
	print_hist("UsedSS (number of petals)",1,histogram(usedSS,1));

	// Histogram of SF
	Table usedSF = A.used_by_kind("SF",G,pp,F);
	print_hist("UsedSF (number of petals)",1,histogram(usedSF,1));

	// Some petal
	for (int i=1; i<=10; i++) print_table("Petal of plate "+i2s(i),A.infos_petal(1000*i,5,G,P,pp,F));

	// Percentage of fiber assigned
	printf("  %s assignments in total (%.4f %% of all fibers)\n",f(A.na()).c_str(),percent(A.na(),Nplate*Nfiber));

	// Some stats
	Table done = initTable(Categories,MaxObs+1);
	List fibers_used = initList(Categories);
	List targets = initList(Categories);
	for (int id=0; id<Categories; id++) {
		List hist = hist2[id];
		for (int i=0; i<=MaxObs; i++) {
			//i here is number of observations lacking  i = goal -done
			//done=goal -i for i<=goal
			//0 for done>goal, i.e. i<0
			if (i<=F.goal[id]) {
				done[id][i]=hist[F.goal[id]-i];
				fibers_used[id]+=i*done[id][i];
				targets[id] += done[id][i];
			}
		}
	}

	Dtable stats = initDtable(Categories,4);
	for (int id=0; id<Categories; id++) {
		stats[id][0] = (targets[id]-done[id][0])/TotalArea;
		stats[id][1] = fibers_used[id]/TotalArea;
		stats[id][2] = targets[id]/TotalArea;
		stats[id][3] = percent(targets[id]-done[id][0],targets[id]);
	}
	print_table("  Id on lines. (Per square degrees) Total, Fibers Used, Available, Percent of observations",stats,latex);
}
