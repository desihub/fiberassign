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

// Counts how many time ob should be observed else more
int nobs(int g, const Gals& G, const Feat& F, const Assignment& A) {
	int cnt(F.goal[G[g].id]);
	for (int i=0; i<Npass; i++) if (A.is_assigned_pg(i,g)) cnt--;
	return cnt;
}

inline double plate_dist(const double theta) {
	// Returns the radial distance on the plate (mm) given the angle,
	// theta (radians).  This is simply a fit to the data provided.
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

// Fast because no big function is called ! Keep it in only 1 function
// For each ~10,000 plate, finds ~25,000 galaxies reachable by the plate,
// Projects them, and for each fiber, finds reachable ones
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, Time &t) {
	List permut = random_permut(MaxPlate);
	int i;
	//omp_set_num_threads(24);
#pragma omp parallel
	{ 	int id = omp_get_thread_num();
		double quarter = floor(MaxPlate/(4*omp_get_num_threads()));
		double cnt, avg, std; int max = 0; int min = 1e7; 
		// Collects for each plate ; shuffle order of plates (useless !?)
		for (i=id; i<MaxPlate; i++) { // <- begins at id, otherwise all begin at 0 -> conflict. Do all plates anyway
			int j = permut[i];
			plate p = P[j];
			// Takes neighboring ~25000 galaxies that can be reached by this plate
			std::vector<int> nbr = T.near(G,p.nhat,PlateRadius*M_PI/180.); // nbr for neighbours
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
			for (int k=0; k<MaxFiber; k++) {
				std::vector<int> gals = kdT.near(&(pp.fp[2*k]),PatrolMinRad,PatrolMaxRad);
				P[j].av_gals[k] = initList(gals);
			}
			// Avancement 
			if (cnt==quarter && id==0) {printf("  Thread 0 has done 1/4 of his job"); print_time(t," at");}
		}
		//printf("  Thread %2d finished",id); print_time(t," at");
	} // End parallel
}

void collect_available_tilefibers(Gals& G, const Plates& P) {
	for(int j=0; j<MaxPlate; j++) {
		for(int k=0; k<MaxFiber; k++) {
			for(int m=0; m<P[j].av_gals[k].size(); m++) {
				int i = P[j].av_gals[k][m];
				G[i].av_tfs.push_back(pair(j,k));
			}
		}
	}
}

// (On plate p) finds if there is a collision if fiber k would watch at galaxy g (collision with neighb)
int find_collision(int j, int k, int g, const PP& pp, const Gals& G, const Plates& P, const Assignment& A) {
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
			if (sq(x-xn,y-yn) < sq(collide)) return kn;
		}
	}
	return -1;
}

bool ok_assign_g_to_jk(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Feat& F, const Assignment& A) {
	// No collision, only ELG in last pass
	if (find_collision(j,k,g,pp,G,P,A)!=-1) return false;
	if (P[j].ipass==MaxPass-1 && G[g].kind(F)!="ELG") return false; // Only ELG at the last pass
	//if (G[g].kind(F)!="SS" && G[g].kind(F)!="SF" && A.unused_fbp(j,k,pp)<=MinUnused) return false;
	if (G[g].kind(F)=="SS" && A.nkind(j,k,"SS",G,P,pp,F)>=MaxSS) return false;
	if (G[g].kind(F)=="SF" && A.nkind(j,k,"SF",G,P,pp,F)>=MaxSF) return false;
	return true;
}

bool ok_assign_tot(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Feat& F, const Assignment& A) {
	// No collision, only ELG in last pass
	if (!ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) return false;
	if (A.is_assigned_tf(j,k)) return false;
	if (A.is_assigned_pg(P[j].ipass,g)) return false;
	return true;
}

// Assignment "global" -------------------------------------------------------------------------------------

// Assign fibers naively
void assign_fibers(const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	for (int ipass=0; ipass<MaxPass; ipass++) {
		List randPlates = random_permut(MaxPlate);
		for (int i=0; i<MaxPlate; i++) {
			int j = randPlates[i];
			if (P[j].ipass == ipass) {
				List randFibers = random_permut(MaxFiber);
				for (int m=0; m<MaxFiber; m++) { // Fiber
					int k = randFibers[m];
					int best=-1; float minp=1e30;
					// Need to shuffle available galaxies
					std::vector<int> av_gals = P[j].av_gals[k];
					List randGals = random_permut(av_gals.size());
					for (int n=0; n<av_gals.size(); n++) { 
						int g = av_gals[randGals[n]];
						int m = nobs(g,G,F,A);
						int prio = F.prio[G[g].id];
						//forbid Ly-a in gray time, i.e. pass=5
						//forbid LRGs also 1/3/15 rnc
						if(!A.is_assigned_pg(ipass,g) && m>0 &&((G[g].kind(F)=="ELG")||(P[j].ipass!=Npass-1))) {
							if (prio<minp || (best>=0 && m>nobs(best,G,F,A))) { // !! G[-1] isn't initialized
								int kp = find_collision(j,k,g,pp,G,P,A);
								if (kp==-1) {
									minp = prio;
									best = g;
								}
							}
						}
					}
					// Assign best galaxy to this fiber
					if (best!=-1 && (nobs(best,G,F,A)>0)) A.assign(j,k,best,G,P,pp);
					}}}}
	printf("  %s assignments\n",f(A.na()));
}

// Try to use unused fibers by reassigning some used ones
// Before : (jp,kp) <-> g ; (j,k) & gp free
// After : (j,k) <-> g & (jp,kp) <-> gp
// Prints the number of additional assigned galaxies
void improve(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int na_start = A.na();
	//List randPlates = random_permut(MaxPlate);
	for (int jj=0; jj<MaxPlate; jj++) {
		int j = MaxPlate-jj-1;
		for(int k=0; k<MaxFiber; k++) {
			if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
				bool finished(false);
				std::vector<int> av_g = P[j].av_gals[k];
				for (int i=0; i<av_g.size() && !finished; i++) { // Not shuffled
					int g = av_g[i]; // g : possible galaxy for (j,k)
					// Is it allowed for jk to take g ?
					if (ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) {
						// What tfs have taken g ? Could they take someone else ?
						std::vector<pair> tfs = A.chosen_tfs(g);
						for (int p=0; p<tfs.size() && !finished; p++){
							int jp = tfs[p].f;
							int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
							std::vector<int> av_g2 = P[jp].av_gals[kp];
							for(int m=0; m<av_g2.size() && !finished; m++) {
								int gp = av_g2[m]; // gp : another possibility for (jp,kp)
								if (nobs(gp,G,F,A)>=1 && ok_assign_g_to_jk(gp,jp,kp,P,G,pp,F,A) && !A.is_assigned_pg(P[jp].ipass,g)) {
									// Modify assignment
									A.unassign(jp,kp,g,G,P,pp);
									A.assign(j,k,g,G,P,pp);
									A.assign(jp,kp,gp,G,P,pp);
									finished = true; 
								}
								//if(nobs(gp,G,F,A)>0){
									//A.unassign(jp,kp,g,P);
									//if(ok_assign_g_to_jk(gp,jp,kp,P,G,pp,F,A)){ 
										//Modify assignment
											//A.assign(j,k,g,P);
										//A.assign(jp,kp,gp,P);
										//finished = true; }
									//else {
										//A.assign(jp,kp,g,P);
									//}
								//}
							}}}}}}}
	int na_end(A.na());
	printf("  %s more assignments (%f %% improvement)\n",f(na_end-na_start),percent(na_end-na_start,na_start));
}

// Redistrubte assignments trying to get at least 500 free fibers on each plate/tile
// Redo so we look first at plates with too few free fibers
// Before : (j,k) <-> g, with Sp(k) too much used
// After : (jreassign,kreassign) <-> g & (j,k) free, such that (jreassign,kreassign) comes from most unused (ji,ki)
// Reassign if this is improvement
// The problem must be that we don't recompute unused_fbp at each reassignment 
void redistribute(const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	List Sp = pp.spectrom;
	Table unused_fbp = A.unused_fbp(pp);
	//print_table("unused",unused_fbp);
	List randPlates = random_permut(MaxPlate);
	// Consider all petals with too few free fibers
	for (int jj=0;jj<MaxPlate;jj++) {
		int j = randPlates[jj];
		List randFibers = random_permut(MaxFiber);
		//for(int k=0;k<MaxFiber && (unused_fibers[j]<minUnused+1);++k){
		for (int kk=0; kk<MaxFiber; kk++) {
			int k = randFibers[kk];
			//int nfree = unused_fbp[j][Sp[k]];
			int nfree = A.unused_fbp(j,k,pp);
			int g = A.TF[j][k];
			if (g!=-1 && nfree<MinUnused/* && G[g].id==4*/) { // Only ELG
				// Consider other ways to observe this galaxy
				std::vector<pair> tfs = G[g].av_tfs;
				int jreassign(-1); int kreassign(-1); int mostunused(-1);
				for (int c=0; c<tfs.size(); c++) {
					int jp = tfs[c].f;
					int kp = tfs[c].s;
					//int nfreep = unused_fbp[jp][Sp[kp]];
					int nfreep = A.unused_fbp(jp,kp,pp);
					// Use freest plate and petal
					if (nfreep>mostunused && ok_assign_tot(g,jp,kp,P,G,pp,F,A)) {
						mostunused = nfreep;
						jreassign = jp;
						kreassign = kp; 
					}
				}
				if (mostunused > nfree) {
					//printf("%d %d %d %d %d %d %d\n",mostunused, nfree, j,k,jreassign,kreassign,g);
					A.unassign(j,k,g,G,P,pp);
					A.assign(jreassign,kreassign,g,G,P,pp);
					//unused_fbp[j][Sp[k]]++;
					//unused_fbp[jreassign][Sp[kreassign]]--;		
					redistributions++; 
				}}}}
	printf("  %s redistributions (~%.4f %% redistributed)\n",f(redistributions),percent(redistributions,A.na()));
}

// Assignment "for one" --------------------------------------------------------------------------------

void assign_fibers_for_one(int j, const Gals& G, const Plates& P, const PP& pp, const Feat& F, Assignment& A) {
	List randFibers = random_permut(MaxFiber);
	int as(0);
	for (int kk=0; kk<MaxFiber; kk++) {
		int k = randFibers[kk];
		int best = -1; int mbest = -1; float minp = 1e30;
		List av_gals = P[j].av_gals[k];
		List randGals = random_permut(av_gals.size());
		for (int gg=0; gg<av_gals.size(); gg++) { 
			int g = av_gals[randGals[gg]];
			int prio = G[g].prio(F);
			int m = nobs(g,G,F,A);
			if (ok_assign_tot(g,j,k,P,G,pp,F,A)) {
				if (!A.once_obs(g)) {
					if (m>=1 && prio<minp) {
						best = g;
						minp = prio;
						mbest = m;
					}
				}
				else if (m>=1 && prio<minp || best!=-1 && m>mbest) {
					best = g;
					minp = prio;
					mbest = m;
				}
			}
		}
		// Assign best galaxy to this fiber
		if (best!=-1 && nobs(best,G,F,A)>=1) { A.assign(j,k,best,G,P,pp); as++; }
	}
	printf(" %5s assignments",f(as));
}

void improve_for_one(int j, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int improvement(0);
	List randFibers = random_permut(MaxFiber);
	for (int kk=0; kk<MaxFiber; kk++) {
		int k = randFibers[kk];
		if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
			bool finished(false);
			List av_g = P[j].av_gals[k];
			for (int i=0; i<av_g.size() && !finished; i++) { // Not shuffled
				int g = av_g[i]; // g : possible galaxy for (j,k)
				// Is it allowed for jk to take g ?
				if (ok_assign_g_to_jk(g,j,k,P,G,pp,F,A)) {
					// What tf have taken g on this plate ? Could it take someone else ?
					pair tf = A.PG[P[j].ipass][g];
					if (tf.f != j) printf("\n !!! ERROR, confusion in pass (improve_for_one)\n");
					int kp = tf.s; // (j,kp) currently assigned to galaxy g
					if (!tf.isnull()) {
						List av_g2 = P[j].av_gals[kp];
						for(int m=0; m<av_g2.size() && !finished; m++) {
							int gp = av_g2[m]; // gp : another possibility for (jp,kp)
							if (ok_assign_tot(gp,j,kp,P,G,pp,F,A) && nobs(gp,G,F,A)>=1 /*>=0 if never observed yet*/ && A.unused_fbp(j,k,pp)>MinUnused && A.unused_fbp(j,kp,pp)>MinUnused) {
								// Modify assignment
								A.unassign(j,kp,g,G,P,pp);
								A.assign(j,k,g,G,P,pp);
								A.assign(j,kp,gp,G,P,pp);
								improvement++;
								finished = true; }}}}}}}
	printf(" %5s more",f(improvement));
}
					
void redistribute_for_one(int j, const Gals& G, const Plates&P, const PP& pp, const Feat& F, Assignment& A) {
	int redistributions(0);
	List Sp = pp.spectrom;
	List randFibers = random_permut(MaxFiber);
	List changed_gals = initList(Ngal,0); // Can be opt with av_gals_plate
	for (int kk=0; kk<MaxFiber; kk++) {
		int k = randFibers[kk];
		int g = A.TF[j][k];
		if (g!=-1 && changed_gals[g]==0) { // Don't change if already changed
			// Consider other ways to observe this galaxy
			std::vector<pair> tfs = G[g].av_tfs;
			bool finished(false);
			for (int c=0; c<tfs.size() && !finished; c++) {
				int jp = tfs[c].f;
				int kp = tfs[c].s;
				if (jp==j && ok_assign_tot(g,j,kp,P,G,pp,F,A) && A.unused_fbp(j,kp,pp)>MinUnused) {
					A.unassign(j,k,g,G,P,pp);
					A.assign(j,kp,g,G,P,pp);
					changed_gals[g] = 1;
					redistributions++; 
					finished = true; }}}}
	printf(" %5s redistributions",f(redistributions));
}

// Other useful functions ---------------------------------------------------------------------

// Detemine how many galaxies need to be dropped to guarantee 40 free fibers for each petal
// Count how many free fibers there are beyond 500 in each plate
void print_free_fibers(const PP& pp, const Assignment& A) {
	printf("# Free fibers statistics\n");
	// Prints histogram of free fibers
	Table unused_fbp = A.unused_fbp(pp);
	print_hist("  Number of petals with this many free fiber",5,histogram(unused_fbp,5));
}

void results_on_inputs(const Gals& G, const Plates& P, const Feat& F) {
	printf("# Results on inputs :\n");

	// Print features
	print_list("  Kinds corresponding :",F.kind);
	print_list("  Priorities :",F.prio);
	print_list("  Goals of observations :",F.goal);

	// How many galaxies in range of a fiber ?
	print_list("  How many galaxies in range of a fiber :",gals_range_fibers(P));

	// Stats on number of av gals per plate and per fiber
	if (false) { // Takes a lot of time !
	int avg(0); int std(0); int min(1e7); int max(0);
	int avgf(0); int stdf(0); int minf(1e7); int maxf(0);
	for (int j=0; j<MaxPlate; j++) {
		printf(" %d",j); std::cout.flush();
		// Plates
		List L = P[j].av_gals_plate();
		int n = L.size();
		avg += n; std += n*n;
		if (n>max) max=n;
		if (n<min) min=n;
		// Fibers
		for (int k=0; k<MaxFiber; k++) {
			int nf = P[j].av_gals[k].size();
			avgf += nf; stdf += nf*nf;
			if (nf>maxf) maxf=nf;
			if (nf<minf) minf=nf;
		}
	}
	print_stats("  Available galaxies for a plate :",MaxPlate,avg,std,min,max);
	printf("  %d %d ",avg,std);
	print_stats("  Available galaxies for a fiber :",MaxPlate*MaxFiber,avgf,stdf,minf,maxf);

	// Stats on number of av tfs per galaxy
	int avgg(0); int stdg(0); int ming(1e7); int maxg(0);
	for (int g=0; g<Ngal; g++) {
		int n = G[g].av_tfs.size();
		avgg += n; stdg += n*n;
		if (n>maxg) maxg=n;
		if (n<ming) ming=n;
	}
	print_stats("  Available tile-fibers for a galaxy :",Ngal,avgg,stdg,ming,maxg);
	}
}


void display_results(const Gals& G, const Plates& P, const PP& pp, const Feat& F, const Assignment& A) {
	printf("# Results :\n");

	// Histogram of SS
	Table usedSS = A.used_by_kind("SS",G,pp,F);
	print_hist("  UsedSS (number of petals)",1,histogram(usedSS,1));

	// Histogram of SF
	Table usedSF = A.used_by_kind("SF",G,pp,F);
	print_hist("  UsedSF (number of petals)",1,histogram(usedSF,1));

	// Raw numbers of galaxies by id and number of remaining observations
	Table hist2 = initTable(Categories,MaxObs+1);
	for (int g=0; g<Ngal; g++) {
		int n = nobs(g,G,F,A);
		if (n>=0 && n<=MaxObs) hist2[G[g].id][n]++;
		else printf(" !!! Error in display_result : observation beyond limits\n");
	}
	print_table("  Remaining observations (id on lines, nobs left on rows)",hist2,true);

	// Per sq deg in tex format
	Table done = initTable(Categories,MaxObs+1);
	List fibers_used = initList(Categories);
	List targets = initList(Categories);

	for (int id=0; id<Categories; id++) {
		List hist = initList(hist2[id]);
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
	printf("\n tex format, per sq deg\n");
	printf("id| obsv'd  0         1          2          3          4          5        total fibers used    avail   pct obsvd \n");
	for (int id=0; id<Categories; id++) {
		printf("%2d&",id);
		for (int i=0; i<=F.goal[id]; i++) printf("%10.0f &",done[id][i]/TotalArea);
		for (int i=F.goal[id]+1; i<MaxObs+1; i++) printf("%10d &",0);
		printf("%8.1f &%8.1f &%8.1f &%8.4f \\ \n",(targets[id]-done[id][0])/TotalArea,fibers_used[id]/TotalArea,targets[id]/TotalArea,percent(targets[id]-done[id][0],targets[id]));
	}
}

void plot_freefibers(str s, const Plates& P, const Assignment& A) {
	printf("# Plot free fibers positions\n");
	double pi = 3.1415926535;
	List unused_fibers = A.unused_f();
	FILE * fplot;
	const char * c = s.c_str();
	fplot = fopen(c,"w");
	fprintf(fplot,"j --- ra --- dec --- unused fibers \n");
	for (int j=0; j<MaxPlate; j++) {
		double phi = atan2(P[j].nhat[1],P[j].nhat[0]);
		double theta = acos(P[j].nhat[2]);
		double ra = 180.*phi/pi;
		double dec = 90.*(1.-2.*theta/pi);
		fprintf(fplot,"%d  %2d  %2d  %4d \n",j,ra,dec,unused_fibers[j]);
	}
	fclose(fplot);
}

// how many galaxies have been observed how many times at some point?
// counter(id, no. of observations, times observed that many times)
// allow up to 5 observations
// cannot rely on G.nobs because this is result at end
void time_count(int jmin, int jmax, const Gals& G, Table& nc, List& tss, const Assignment& A){
	int nstep = MaxPlate/ntimes;
	for (int j=jmin; j<jmax; ++j) { // For these tiles
		for (int k=0;k<MaxFiber;++k) { // and all fibers
			int g = A.TF[j][k];
			if (g!=-1) { // If this tile-fiber has an assigned galaxy
				tss[g]++; // tss[g] >= 1
				int times_observed = tss[g];
				int id = G[g].id; 
				// Increase the counts for this id and times_observed
				nc[id][times_observed]++; 
				if(times_observed>=2) nc[id][times_observed-1]--;
				//print_list("\n",initList(nc[id]));
			}
		}
	}
}
// nc should be only a list(Categories)

// time_count counted the number of times a galaxy was observed
void time_line(const Gals& G, const Feat& F, const Assignment& A) {
	//produce numbers for python plot showing what fraction of galaxies of each sort are observed as a function of 'time'
	//what fraction of galaxies has been observed in each category as a function of "time"?
	//output for python plotting

	int curves=Categories+4+1;//3 for Ly-a, 1 for LRG
	int nstep=MaxPlate/ntimes;
	Table nc = initTable(Categories+1,MaxObs+1); // new_count
	Table latest_count = initTable(ntimes+1,curves);
	List tss = initList(Ngal);
	//cannot rely on G.nobs because this is result at end
	//write for python
	std::cout<< "  fraction of galaxies observed, by id, and as function of plates observed\n";
	printf("  suitable for python plot\n");
	printf("x=[");
	for(int i=1;i<ntimes+1; i++) printf("%d,",nstep);
	printf("] \n \n");
	//end python write
	for (int i=0;i<ntimes;++i){
		int m(0);
		time_count(i*nstep,(i+1)*nstep,G,nc,tss,A);
		for(int j=1; j<=Categories; ++j) {
			for (int k=1; k<=F.goal[j]; ++k) {
				m++;
				latest_count[i][m]=nc[j][k];
			}
		}
	}
	for (int m=1; m<curves; ++m) {
		double running_total=0;
		printf("y%d=[",m);
		for(int i=0;i<ntimes;++i) {
			running_total = double(latest_count[i][m])/double(G.size());
			if(i!=0) printf("  , ");
			printf("%.5f",running_total);
		}
		printf("]\n");
	}
	printf("# End of python data\n");
}

// Make list of all conflicts assuming no three-way conflicts
Table conflicts(const Gals& G, const Plates& P, const PP& pp, const Assignment& A) {
	Table Q;
	for(int j=0;j<MaxPlate; j++){
		for(int k=0;k<MaxFiber; k++){
			if(A.is_assigned_tf(j,k)){
				List triplet;
				int g = A.TF[j][k]; 
				int kp = find_collision(j,k,g,pp,G,P,A);
				int gp = A.TF[j][kp];
				if(kp!=-1 && k<kp && gp!=-1) { 
					triplet.push_back(j);triplet.push_back(k);triplet.push_back(kp);
					Q.push_back(triplet);
					printf("!!! Conflict on plate j = %d : (g,k) = (%d,%d) with (gp,kp) = (%d,%d) \n",j,g,k,gp,kp);
				}
			}
		}
	}
	int c(Q.size());
	if (c==0) printf("# No conflict\n");
	else printf("!!! Number of conflicts : %d\n",c);
	return Q;
}
