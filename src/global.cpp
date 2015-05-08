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
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "global.h"

// Collecting information from input -------------------------------------------------------------------------------------
// Fast because no big function is called ! Keep it in only 1 function
// For each ~10,000 plate, finds ~25,000 galaxies reachable by the plate,
// Projects them, and for each fiber, finds reachable ones
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, const Feat& F) {
	Time t;
	init_time(t,"# Begin collecting available galaxies");
	List permut = random_permut(F.Nplate);
	double rad = F.PlateRadius*M_PI/180.;
	int jj;
	//omp_set_num_threads(24);
#pragma omp parallel
	{ 	int id = omp_get_thread_num();
		// Collects for each plate ; shuffle order of plates (useless !?)
		for (jj=id; jj<F.Nplate; jj++) { // <- begins at id, otherwise all begin at 0 -> conflict. Do all plates anyway
			int j = permut[jj];
			plate p = P[j];
			// Takes neighboring ~25000 galaxies that can be reached by this plate
			std::vector<int> nbr = T.near(G,p.nhat,rad); // nbr for neighbours
			// Projects thoses galaxies on the focal plane
			Onplates O;
			for (int gg=0; gg<nbr.size(); gg++) {
				int g = nbr[gg];
				struct onplate op = change_coords(G[g],p); 
				op.id = g;
				O.push_back(op);
			}
			debl(id);
			KDtree<struct onplate> kdT(O,2);
			// For each fiber, finds all reachable galaxies
			for (int k=0; k<F.Nfiber; k++) {
				dpair X = pp.coords(k);
				std::vector<int> gals = kdT.near(&(pp.fp[2*k]),0.0,F.PatrolRad);
				for (int g=0; g<gals.size(); g++) {
					struct onplate op = change_coords(G[gals[g]],p); 
					dpair Xg = dpair(op.pos[0],op.pos[1]);
					if (sq(Xg,X)<sq(F.PatrolRad)/*Needed*/) P[j].av_gals[k].push_back(gals[g]);
				}
			}
		}
	} // End parallel
	print_time(t,"# ... took :");
}

void collect_available_tilefibers(Gals& G, const Plates& P, const Feat& F) {
	Time t;
	init_time(t,"# Begin computing available tilefibers");
	for(int j=0; j<F.Nplate; j++) {
		for(int k=0; k<F.Nfiber; k++) {
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
	if (A.find_collision(j,k,g,pp,G,P,F)!=-1) return false;
	if (kind==F.ids.at("SF") && A.nkind(j,k,F.ids.at("SF"),G,P,pp,F)>=F.MaxSF) return false;
	if (kind==F.ids.at("SS") && A.nkind(j,k,F.ids.at("SS"),G,P,pp,F)>=F.MaxSS) return false;
	if (P[j].ipass==F.Npass-1 && !(kind==F.ids.at("ELG") || kind==F.ids.at("SS") || kind==F.ids.at("SF"))) return false; // Only ELG, SF, SS at the last pass
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
		int prio = fprio(g,G,F,A);
		int m = A.nobs(g,G,F);
		bool tfb = as_tf ? !A.is_assigned_tf(j,k) : true;
		if (m>=1 && tfb && ok_assign_g_to_jk(g,j,k,P,G,pp,F,A) && A.is_assigned_jg(j,g,F)==-1 && g!=no_g && !isfound(G[g].id,no_kind)) {
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
			int gb = -1; int bb = -1; int jpb = -1; int kpb = -1; int mb = -1; int pb = 1e3;
			List av_g = P[j].av_gals[k];
			for (int i=0; i<av_g.size(); i++) { // Not shuffled
				int g = av_g[i]; // g : possible galaxy for (j,k)
				if (g!=-1 && g!=no_g) {
					// Is it allowed for jk to take g ?
					if (ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) {
						// What tfs have taken g ? Could they take someone else ?
						Plist tfs = A.chosen_tfs(g,F,begin,next);
						for (int p=0; p<tfs.size(); p++) {
							int jp = tfs[p].f;
							int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
							// FIND BEST JP KP !!!
							int best = find_best(jp,kp,G,P,pp,F,A,false,-1,Null()); // best!=g because !A.assigned_pg(best)

							if (best!=-1 && (A.is_assigned_jg(j,g,F)==-1 || jp==j)) {
								int prio = fprio(best,G,F,A);
								int m = A.nobs(best,G,F);
								if (prio<pb || A.once_obs[g] && m>mb) {
									gb = g; bb = best; jpb = jp; kpb = kp; mb = m; pb = prio;
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
		int gb = -1; int gpb = -1; int bb = -1; int kpb = -1; int mb = -1; int pb = 1e3;
		List av_gals = P[j].av_gals[k];
		for (int i=0; i<av_gals.size(); i++) {
			int g = av_gals[i];
			if (G[g].id==id && A.find_collision(j,k,g,pp,G,P,F)==-1 && A.is_assigned_jg(j,g,F)==-1 && A.nobs(g,G,F)>=1) {
				for (int kkp=0; kkp<fibskind.size(); kkp++) {
					int kp = fibskind[kkp];
					int gp = A.TF[j][kp];
					int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
					if (best!=-1) {
						int prio = fprio(best,G,F,A);
						int m = A.nobs(best,G,F);
						if (prio<pb || A.once_obs[g] && m>mb) {
							gb = g; gpb = gp; bb = best; kpb = kp; mb = m; pb = prio;
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
	int n = next==-1 ? F.Nplate-j0 : next; // Not F.Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	List randPlates = F.Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randFibers = random_permut(F.Nfiber);
		for (int kk=0; kk<F.Nfiber; kk++) { // Fiber
			int k = randFibers[kk];
			assign_fiber(j,k,G,P,pp,F,A);
		}
	}
	str next_str = next==-1 ? "all left" : f(n);
	printf("  %s assignments on %s next plates\n",f(A.na(F,j0,n)).c_str(),next_str.c_str());
	if (next!=1) print_time(t,"# ... took :");
}

void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin improve :");
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-j0 : next;
	int na_start = A.na(F,j0,n);
	List plates = sublist(j0,n,A.order);
	List randPlates = F.Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) for (int k=0; k<F.Nfiber; k++) improve_fiber(j0,n,randPlates[jj],k,G,P,pp,F,A);
	int na_end = A.na(F,j0,n);
	printf("  %s more assignments (%.3f %% improvement)\n",f(na_end-na_start).c_str(),percent(na_end-na_start,na_start));
	if (next!=1) print_time(t,"# ... took :");
}

void improve_from_kind(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, str kind, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin improve "+kind+" :");
	int id = F.ids.at(kind);
	List no_kind; no_kind.push_back(id);
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-A.next_plate : next;
	int na_start(A.na(F,j0,n));
	List plates = sublist(j0,n,A.order);
	List randPlates = F.Randomize ? random_permut(plates) : plates;
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randPetals = random_permut(F.Npetal);
		for (int ppet=0; ppet<F.Npetal; ppet++) {
			int p = randPetals[ppet];
			// Take sublist of fibers assigned to kind, and unassigned ones
			List fibskind = A.fibs_of_kind(id,j,p,G,pp,F);
			List fibsunas = A.fibs_unassigned(j,p,G,pp,F);
			for (int kk=0; kk<fibsunas.size(); kk++) {
				// Take an unassigned tf and its available galaxies
				int k = fibsunas[kk];
				int gb = -1; int gpb = -1; int bb = -1; int kpb = -1; int kkpb = -1; int mb = -1; int pb = 1e3;
				List av_gals = P[j].av_gals[k];
				for (int i=0; i<av_gals.size(); i++) {
					int g = av_gals[i];
					if (G[g].id==id && A.nobs(g,G,F)>=1 && A.is_assigned_jg(j,g,F)==-1 && A.find_collision(j,k,g,pp,G,P,F)==-1) {
						// If the av gal is of the good kind, try to assign it, and improving an other one
						for (int kkp=0; kkp<fibskind.size(); kkp++) {
							int kp = fibskind[kkp];
							int gp = A.TF[j][kp];
							int best = find_best(j,kp,G,P,pp,F,A,false,-1,no_kind);
							if (best!=-1) {
								int prio = fprio(best,G,F,A);
								int m = A.nobs(best,G,F);
								if (prio<pb || A.once_obs[g] && m>mb) {
									gb = g; gpb = gp; bb = best; kpb = kp; kkpb = kkp; mb = m; pb = prio;
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
	int na_end(A.na(F,j0,n));
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
	int na_start(A.na(F,j0,n));
	List to_update;
	// Declare that we've seen those galaxies
	A.update_once_obs(j0,F);
	// Get the list of galaxies to update in the plan
	for (int k=0; k<F.Nfiber; k++) {
		int g = A.TF[j0][k];
		// Only if once_obs, we delete all further assignment. obs!=obs_tmp means that the galaxy is a fake one for example (same priority but different goal)
		if (g!=-1 && A.nobsv_tmp[g]!=A.nobsv[g] && A.once_obs[g]) to_update.push_back(g);
	}
	// Update information on previously seen galaxies
	A.update_nobsv_tmp_for_one(j0,F);
	// Update further in the plan
	for (int gg=0; gg<to_update.size(); gg++) {
		int g = to_update[gg];
		Plist tfs = A.chosen_tfs(g,F,j0+1,n-1); // Begin at j0+1, because we can't change assignment at j0 (already watched)
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
	int na_end(A.na(F,j0,n));
	printf(" %4d unas & %4d replaced\n",cnt,na_end-na_start+cnt); fl();
}

// If no enough SS and SF, remove old_kind an replace to SS-SF (new_kind) on petal (j,p)
void replace(int old_kind, int new_kind, int j, int p, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	int m = A.nkind(j,p,new_kind,G,P,pp,F,true);
	List fibskindd = A.fibs_of_kind(old_kind,j,p,G,pp,F);
	List fibskind = random_permut(fibskindd);
	int Max = new_kind==F.ids.at("SS") ? F.MaxSS : F.MaxSF;
	while (m<Max && fibskind.size()!=0) {
		bool fin(false);
		int k = fibskind[0];
		List av_g = P[j].av_gals[k];
		for (int gg=0; gg<av_g.size() && !fin; gg++) {
			int g = av_g[gg];
			if (G[g].id==new_kind && A.find_collision(j,k,g,pp,G,P,F)==-1 && A.nobs(g,G,F)>=1 && pp.spectrom[k]==p) {
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
	int n = next==-1 ? F.Nplate-j0 : next; // Not F.Nplate-A.next_plate+1
	List plates = sublist(j0,n,A.order);
	List randPlates = F.Randomize ? random_permut(plates) : plates;
	List no_kind; no_kind.push_back(F.ids.at("QSOLy-a")); no_kind.push_back(F.ids.at("QSOTracer")); no_kind.push_back(F.ids.at("FakeQSO")); no_kind.push_back(F.ids.at("SS")); no_kind.push_back(F.ids.at("SF")); 
	List no_reg; no_reg.push_back(F.ids.at("QSOLy-a")); no_reg.push_back(F.ids.at("QSOTracer"));  no_reg.push_back(F.ids.at("LRG")); no_reg.push_back(F.ids.at("ELG")); no_reg.push_back(F.ids.at("FakeQSO")); no_reg.push_back(F.ids.at("Fake LRG"));
	List only_qso; only_qso.push_back(F.ids.at("LRG")); only_qso.push_back(F.ids.at("ELG")); only_qso.push_back(F.ids.at("Fake LRG")); only_qso.push_back(F.ids.at("SS")); only_qso.push_back(F.ids.at("SF"));
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randPetals = random_permut(F.Npetal);
		for (int ppet=0; ppet<F.Npetal; ppet++) {
			// Assign only QSO
			int p = randPetals[ppet];
			List fibs = pp.fibers_of_sp[p];
			List randFibers = random_permut(fibs);
			for (int kk=0; kk<F.Nfbp; kk++) { // Fiber
				int k = randFibers[kk];
				assign_fiber(j,k,G,P,pp,F,A,-1,only_qso);
			}
			// Assign without SS and SF
			for (int kk=0; kk<F.Nfbp; kk++) { // Fiber
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
	printf("  %s assignments on %s next plates\n",f(A.na(F,j0,n)).c_str(),next_str.c_str());
	if (next!=1) print_time(t,"# ... took :");
}

void redistribute_tf(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A, int next) {
	Time t;
	if (next!=1) init_time(t,"# Begin redistribute TF :");
	int j0 = A.next_plate;
	int n = next==-1 ? F.Nplate-A.next_plate : next;
	List plates = sublist(j0,n,A.order);
	List randPlates = F.Randomize ? random_permut(plates) : plates;
	int red(0);
	Table Done = initTable(F.Nplate,F.Nfiber);
	for (int jj=0; jj<n; jj++) {
		int j = randPlates[jj];
		List randFiber = random_permut(F.Nfiber);
		for (int kk=0; kk<F.Nfiber; kk++) {
			int k = randFiber[kk];
			if (Done[j][k]==0) {
				int g = A.TF[j][k];
				if (g!=-1) {
					bool finished = false;
					Plist av_tfs = G[g].av_tfs;
					for (int i=0; i<av_tfs.size() && !finished; i++) {
						int jp = av_tfs[i].f;
						int kp = av_tfs[i].s;
						if (j0<=jp && jp<j0+n && !A.is_assigned_tf(jp,kp) && Done[jp][kp]==0 && ok_assign_g_to_jk(g,jp,kp,P,G,pp,F,A) && A.is_assigned_jg(jp,g,F)==-1) {
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
	print_list("  How many galaxies in range of a fiber :",gals_range_fibers(P,F));

	// 1 Histograms on number of av gals per plate and per fiber
	Cube T = initCube(F.Categories,F.Nplate,F.Nfiber);
	for (int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			List gals = P[j].av_gals[k];
			for (int g=0; g<gals.size(); g++) T[G[gals[g]].id][j][k]++;
		}
	}
	Table hist1;
	for (int id=0; id<F.Categories; id++) hist1.push_back(histogram(T[id],1));
	print_mult_table_latex("Available galaxies (by kind) for a TF",outdir+"avgalhist.dat",hist1,1);

	// 2 Histograms on number of av tfs per galaxy
	Table Tg = initTable(F.Categories,0);
	for (int g=0; g<F.Ngal; g++) {
		int n = G[g].av_tfs.size();
		Tg[G[g].id].push_back(n);
	}
	Table hist2;
	for (int id=0; id<F.Categories; id++) hist2.push_back(histogram(Tg[id],1));
	print_mult_table_latex("Available tile-fibers for a galaxy (by kind)",outdir+"avtfhist.dat",hist2,1);

	// 3 Histogram of number of times (by different plates) reachable galaxies
	List countgals;
	List countgals_nopass;
	for (int g=0; g<F.Ngal; g++) {
		int id = G[g].id;
		List plates;
		for (int i=0; i<G[g].av_tfs.size(); i++) {
			int j = G[g].av_tfs[i].f;
			if (!isfound(j,plates) && P[j].ipass!=F.Npass-1) plates.push_back(j);
		}
		countgals.push_back(plates.size());
	}
	Dtable countstot;
	List h00 = histogram(countgals,1);
	countstot.push_back(percents(h00,sumlist(h00))); 
	print_mult_Dtable_latex("Histogram of percents (by different plates) reachable galaxies (not 5th pass)",outdir+"reachplate.dat",countstot,1);
	
	// By kind
	//Table countgals = initTable(F.Categories,0);
	//Table countgals_nopass = initTable(F.Categories,0);
	//for (int g=0; g<F.Ngal; g++) {
		//int id = G[g].id;
		//List plates; List plates_no;
		//for (int i=0; i<G[g].av_tfs.size(); i++) {
			//int j = G[g].av_tfs[i].f;
			//if (!isfound(j,plates)) {
				//plates.push_back(j);
				//if (P[j].ipass!=F.Npass-1 || F.id("ELG")==id) plates_no.push_back(j);
			//}
		//}
		//countgals[id].push_back(plates.size());
		//countgals_nopass[id].push_back(plates_no.size());
	//}
	//Table countstot; 
	//Table countstot_nopass; 
	//for (int i=0; i<F.Categories; i++) {
		//countstot.push_back(histogram(countgals[i],1)); 
		//countstot_nopass.push_back(histogram(countgals_nopass[i],1));
	//}
	//print_mult_table_latex("Histogram of number of times (by different plates) reachable galaxies",outdir+"reachplate.dat",countstot,1);
	//print_mult_table_latex("Histogram of number of times (by different plates) reachable galaxies without last pass",outdir+"reachplateno.dat",countstot_nopass,1);
}

void display_results(str outdir, const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A, bool latex) {
	printf("# Results :\n");
	// Some petal
	//for (int i=1; i<=10; i++) print_table("Petal of plate "+i2s(i),A.infos_petal(1000*i,5,G,P,pp,F));
	
	// Collision histogram of distances between galaxies
	Dlist coldist;
	for (int j=0; j<F.Nplate; j++) {
		List done = initList(F.Nfiber);
		for (int k=0; k<F.Nfiber; k++) {
			if (done[k] == 0) {
				int c = A.is_collision(j,k,pp,G,P,F);
				if (c!=-1) {
					done[c] = 1;
					struct onplate op = change_coords(G[A.TF[j][k]],P[j]);
					dpair G1 = dpair(op.pos[0],op.pos[1]);
					struct onplate opn = change_coords(G[A.TF[j][c]],P[j]);
					dpair G2 = dpair(opn.pos[0],opn.pos[1]);
					coldist.push_back(norm(G2-G1));
				}
			}
		}
	}
	double intervaldist = 0.01;
	
	Dlist histcoldist = histogram(coldist,intervaldist);
	Dlist redhistcol = percents(histcoldist,sumlist(histcoldist));
	Dtable Dtd; Dtd.push_back(redhistcol); Dtd.push_back(cumulate(redhistcol));
	print_mult_Dtable_latex("Collision histogram of distances between galaxies",outdir+"coldist.dat",Dtd,intervaldist);

	// 0 List of total number of galaxies
	List tots = initList(F.Categories);
	for (int g=0; g<F.Ngal; g++) tots[G[g].id]++;
	
	// 1 Raw numbers of galaxies by id and number of remaining observations
	int MaxObs = max(F.goal);
	Table hist2 = initTable(F.Categories,MaxObs+1);
	Table ob = initTable(F.Categories+1,2*MaxObs+1);

	List maxobs = F.maxgoal();
	for (int g=0; g<F.Ngal; g++) {
		int id = G[g].id;
		int m = A.nobs(g,G,F,false);
		if (!(m>=0 && m<=MaxObs)) { printf("Error obs beyond limits \n"); fl(); }
		if (m>=0 && m<=MaxObs) hist2[id][m]++;
		ob[id][m+MaxObs]++;
	}
	for (int m=0; m<2*MaxObs+1; m++) ob[F.Categories][m] = m-MaxObs;

	Table obsvr = initTable(F.Categories,MaxObs+1);
	for (int i=0; i<F.Categories; i++) {
		int max = F.goal[i];
		for (int j=0; j<=max; j++) obsvr[i][j] = hist2[i][max-j];
	}

	Table with_tots = obsvr;
	for (int i=0; i<F.Categories; i++) {
		int fibs = 0; int obs = 0; int tot =0;
		for (int j=0; j<=MaxObs; j++) tot += obsvr[i][j];
		for (int j=0; j<=MaxObs; j++) fibs += obsvr[i][j]*j;
		for (int j=1; j<=MaxObs; j++) obs += obsvr[i][j];
		with_tots[i].push_back(tot);
		with_tots[i].push_back(fibs);
		with_tots[i].push_back(obs);
	}
	//print_table("  Remaining observations (without negative obs ones)",with_tots,latex,F.kind);
	Dtable obs_per_sqd = ddivide_floor(with_tots,F.TotalArea);

	// 2 Percentages of observation
	Dtable perc = initDtable(F.Categories,2);
	for (int id=0; id<F.Categories; id++) {
		int tot = sumlist(ob[id]);
		int goal = F.goal[id];

		perc[id][0] = percent(tot-ob[id][MaxObs+goal],tot);

		double d(0.0);
		for (int i=0; i<MaxObs; i++) d += ob[id][MaxObs+i]*(goal-i);
		perc[id][1] = percent(d,tot*goal);
	}
	print_table("Obs per sqd and percentages",concatenate(obs_per_sqd,perc),latex,F.kind);

	// 3 Observed galaxies in function of time
	int interval = 10; // Interval of plates for graphic
	Dtable Ttime_scaled = initDtable(F.Categories,0);
	Table Ttime = initTable(F.Categories,0);
	for (int id=0; id<F.Categories; id++) {
		List l;
		int n = 0;
		for (int j=0; j<F.Nplate; j++) {
			for (int p=0; p<F.Npetal; p++) n += A.nkind(j,p,id,G,P,pp,F,true);
			if (j%interval==0) l.push_back(n);
		}
		Ttime[id] = l;
		Ttime_scaled[id] = percents(l,l[l.size()-1]);
	}
	print_mult_Dtable_latex("Observed galaxies in function of time (scaled) (interval 10)",outdir+"time.dat",Ttime_scaled,interval);
	//print_mult_table_latex("Observed galaxies in function of time",,F.kind,Ttime);
	
	/*
	// Lya 1,2,3,4,5, LRG 1,2
	int nk = 9;
	Table Ttim = initTable(nk,0);
	List galaxs = initList(F.Ngal);
	for (int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			int g = A.TF[j][k];
			if (g!=-1) galaxs[g]++;
		}
		if (j%interval==0) {
			List l = initList(9);
			for (int g=0; g<F.Ngal; g++) {
				int n = galaxs[g];
				if (1<=n) {
					if (G[g].id == 0) l[n-1]++;
					if (G[g].id == 2) l[n-1+5]++;
					if (G[g].id == 1) l[n+6]++;
					if (G[g].id == 3) l[n+7]++;
				}
			}
			for (int id=0; id<nk; id++) Ttim[id].push_back(l[id]);
		}
	}
	print_mult_table_latex("Observed galaxies complete (interval 10)",outdir+"time2.dat",Ttim,interval);
	*/

	// 4 Histogramme of percentages of seen Ly-a
	int id = F.ids.at("QSOLy-a");
	int goal = F.goal[id];
	Table Percseen = initTable(goal+1,0);
	for (int g=0; g<F.Ngal; g++) {
		if (G[g].id==id) {
			int n = G[g].av_tfs.size();
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

	// 5 Histogram of time between 2 obs of Ly a
	Table deltas;
	for (int g=0; g<F.Ngal; g++) {
		if (G[g].id == F.ids.at("QSOLy-a")) {
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

	// 6 Free fibers histogram
	Table unused_fbp = A.unused_fbp(pp,F);
	make_square(unused_fbp);
	Table hist0; hist0.push_back(histogram(unused_fbp,1));
	print_mult_table_latex("Number of petals with this many free fiber (interval 1)",outdir+"freefib.dat",hist0,1);

	// 7 Free fibers in function of time (plates)
	List freefibtime = initList(F.Nplate);
	for (int j=0; j<F.Nplate; j++) freefibtime[j] = A.unused_f(j,F);
	Table fft; fft.push_back(freefibtime);
	print_mult_table_latex("Free fibers in function of time (plates)",outdir+"fft.dat",fft);
 
	// 8 Percentage of seen objects as a function of density of objects
	Dcube densities = initDcube(F.Categories+1,0,0);
	for (int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			// For all
			int size = P[j].av_gals[k].size();
			int oc = 0;
			for (int i=0; i<size; i++) if (A.is_assigned_jg(j,P[j].av_gals[k][i])) oc++;
			if (size!=0 && 1<=oc) { 
				double d = percent(oc,size);
				//printf("%f %f %f %d %d %d \n",d,x,invFibArea,size,oc,densities.size()); fl();
				if (size>=densities[F.Categories].size()) densities[F.Categories].resize(size+1);
				densities[F.Categories][size].push_back(d);
			}

			// For kind
			for (int t=0; t<F.Categories; t++) {
				int nkind = 0;
				int ock = 0;
				for (int i=0; i<size; i++) {
					int g = P[j].av_gals[k][i];
					if (G[g].id == t) {
						nkind++;
						if (A.is_assigned_jg(j,g)) ock++;
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
	Dtable densit = initDtable(F.Categories+1,max_row(densities));
	for (int t=0; t<F.Categories+1; t++) for (int i=0; i<densities[t].size(); i++) densit[t][i] = sumlist(densities[t][i])/densities[t][i].size();
	print_mult_Dtable_latex("Perc of seen obj as a fun of dens of objs",outdir+"seendens.dat",densit,1);
	
	// 9 Percentage of TF that have at least 1 QSO av and that are not assigned to a QSO


	// Misc
	// Collision rate
	/*if (Collision)*/ printf("Collision rate : %f \% \n",A.colrate(pp,G,P,F));
	// Histogram of SS
	Table usedSS = A.used_by_kind("SS",G,pp,F);
	print_hist("UsedSS (number of petals)",1,histogram(usedSS,1));

	// Histogram of SF
	Table usedSF = A.used_by_kind("SF",G,pp,F);
	print_hist("UsedSF (number of petals)",1,histogram(usedSF,1));

	// Count
	printf("Count = %d \n",F.Count);

	// Percentage of fibers assigned
	printf("  %s assignments in total (%.4f %% of all fibers)\n",f(A.na(F)).c_str(),percent(A.na(F),F.Nplate*F.Nfiber));
}

void write_FAtile(int j, str outdir, const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A) {
	FILE * FA;
	str s = outdir+"tile"+i2s(j)+".fits";
	FA = fopen(s.c_str(),"w");
	for (int k=0; k<F.Nfiber; k++) {
		int g = A.TF[j][k];
		// k
		fprintf(FA,"%d ",k);
		List av_gals = P[j].av_gals[k];
		// Number of potential galaxies
		fprintf(FA,"%d ",av_gals.size());
		// IDs of potential galaxies
		for (int i=0; i<av_gals.size(); i++) fprintf(FA,"%d ",av_gals[i]);
		// Object type, Target ID, ra, dec, x, y
		if (g!=-1) {
			struct onplate op = change_coords(G[g],P[j]);
			fprintf(FA,"%s %d %f %f %f %f\n",F.kind[G[g].id].c_str(),g,G[g].ra,G[g].dec,op.pos[0],op.pos[1]);
		}
		else fprintf(FA,"-1\n");
	}
	fclose(FA);
}
