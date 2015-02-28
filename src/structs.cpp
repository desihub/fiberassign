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
#include        "structs.h"
#include        "macros.h"
#include        "misc.h"
// Galaxies ------------------------------------------------------------------
void galaxy::print_av_tfs() { // Used to debug
	if (av_tfs.size() == 0) printf("Print av tfs : no av tf \n");
	else {
		printf("Available tfs : ");
		for (int g=0; g<av_tfs.size(); g++) av_tfs[g].print_pair();
		printf("\n");
	}
}
// There is a version of this function (to adapt) for ASCII files, ask to Robert Cahn
Gals read_galaxies(const char fname[], int n) {
	// Read galaxies from binary file--format is ra, dec, z, priority and nobs
	// with ra/dec in degrees
	//  1 = Ly-a QSO ; 2 = QSO tracer ; 3 = LRG ; 4 = ELG ; 5 = fake qsos ; 6 = fake lrgs
	//  Reads every n galaxies (to test quicker)
	Gals P;
	std::ifstream fs(fname,std::ios::binary);
	if (!fs) {  // An error occurred opening the file.
		std::cerr << "Unable to open file " << fname << std::endl;
		myexit(1);
	}
	int Nobj;
	fs.read((char *)&Nobj,sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	std::vector<float> ra,dc,zz,pr;
	std::vector<int>   id,no;
	try {
		ra.resize(Nobj);
		dc.resize(Nobj);
		zz.resize(Nobj);
		id.resize(Nobj);
		pr.resize(Nobj);
		no.resize(Nobj);
		P.reserve(Nobj);
	} catch (std::exception& e) {myexception(e);}
	fs.read((char *)&ra[0],Nobj*sizeof(float));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&dc[0],Nobj*sizeof(float));
	// Could fseek over this, but maybe there will be a use for z someday.
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&zz[0],Nobj*sizeof(float));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&id[0],Nobj*sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&pr[0],Nobj*sizeof(float));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&no[0],Nobj*sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.close();
	// reduce time
	for (int i=0; i<Nobj; i++) {
		if (i%n == 0) {
		double theta = (90 - (double)dc[i])*M_PI/180;
		double phi   = ((double)ra[i]     )*M_PI/180;
		struct galaxy Q;
		Q.nhat[0]    = sin(theta)*cos(phi);
		Q.nhat[1]    = sin(theta)*sin(phi);
		Q.nhat[2]    = cos(theta);
		Q.id         = id[i];
		Q.priority   = pr[i];
		Q.nobs       = no[i];
		// I only added this part to the binary input option.
		// Could add a similar thing to ascii input option if ever necessary [L Samushia].
		Q.ra = ra[i];
		Q.dec = dc[i];
		Q.z = zz[i];
		try{P.push_back(Q);}catch(std::exception& e) {myexception(e);}
		}
	}
	return P;
}

// PP ----------------------------------------------------------------------------
// Read the positions of the fibers on each plate.  Assumes the format
// of fiberpos.txt produced by "randomize_fibers".
// need also to get the petal, i.e. spectrometer  rnc 1/16/15  added S
void PP::read_fiber_positions(const char pos_name[]) {
	std::string buf;
	std::ifstream fs(pos_name);
	if (!fs) { // An error occurred opening the file.
		std::cerr << "Unable to open file " << pos_name << std::endl;
		myexit(1);
	}

	// Reserve some storage, since we expect we'll be reading quite a few
	// lines from this file.
	try {fp.reserve(11000);} catch (std::exception& e) {myexception(e);}
	// Skip any leading lines beginning with # or blank lines.
	getline(fs,buf);
	while (fs.eof()==0 && ( (buf[0]=='#') || (buf.size()==0) )) {
		getline(fs,buf);
	}
	while (fs.eof()==0) {
		double x,y; int fiber,positioner,spectro,remove; 
		std::istringstream(buf) >> fiber >> positioner >> spectro >> x >> y;
		try{
			fp.push_back(x);
			fp.push_back(y);
			spectrom.push_back(spectro);  
		} catch(std::exception& e) {myexception(e);}
		getline(fs,buf);
	}
	fs.close();
	// Verification
	if (MaxFiber != fp.size()/2) {
		std::cerr << "Max Fiber does not correspond to the number of fibers in input file" << std::endl;
		myexit(1);
	}
}

PP::PP() {
	try {N.reserve(400000);} catch (std::exception& e) { myexception(e); }
	N.resize(MaxFiber);
}

void PP::get_neighbors() {
	for(int i=0; i<MaxFiber; i++) {
		for(int j=0; j<MaxFiber; j++) {
			if(i!=j) {
				if(sq(fp[2*i]-fp[2*j],fp[2*i+1]-fp[2*j+1]) < sq(NeighbourRad)) {
					N[i].push_back(j); }}}}
}

// plate ---------------------------------------------------------------------------
void plate::print_plate() const {
	printf("  Plate : %d - pass %d\n",idp,ipass); printf("%f %f %f\n",nhat[0],nhat[1],nhat[2]);
	int r = rand() % MaxFiber;
	print_list("Available galaxies of a random fiber :",initList(av_gals[r]));
}

// Plates -----------------------------------------------------------------------------
// Read positions of the plate centers from an ascii file "center_name", and fill in a structure
// There is a version of this function (to adapt) for non ASCII files, ask to Robert Cahn
Plates read_plate_centers(const char center_name[], int modulo) {
	Plates P;
	std::string buf;
	std::ifstream fs(center_name);
	if (!fs) {  // An error occurred opening the file.
		std::cerr << "Unable to open file " << center_name << std::endl;
		myexit(1);
	}
	// Reserve some storage, since we expect we'll be reading quite a few
	// lines from this file.
	try {P.reserve(4000000);} catch (std::exception& e) {myexception(e);}
	// Skip any leading lines beginning with #
	getline(fs,buf);
	while (fs.eof()==0 && buf[0]=='#') { getline(fs,buf); }
	//read lines until we get some blank lines, then skip blank lines, then get real data
	//this conforms to the format of desi-tiles.par
	while(buf.size()!=0){getline(fs,buf);}
	while(buf.size()==0){getline(fs,buf);}				 
	double ra,dec,ebv,airmass,exposefac; int ipass,in_desi,tileid;
	int l = 0;
	while (fs.eof()==0) {
		getline(fs,buf);
		//gymnastics to read line with string in first column
		std::istringstream ss(buf);
		std::string start;
		ss>> start;
		ss>>tileid>>ra >> dec >> ipass>>in_desi>>ebv>>airmass>>exposefac;	
		//only keep those in footprint, fix ra to lie between 0 and 360
		if (in_desi==1) {
			if (ra<   0.) {ra += 360.;}
			if (ra>=360.) {ra -= 360.;}
			if (dec<-90. || dec>90.) {
				std::cout << "DEC="<<dec<<" out of range reading " << center_name<<std::endl;
				myexit(1);
			}
			double theta = (90.0 - dec)*M_PI/180.;
			double phi   = (ra        )*M_PI/180.;
			struct plate Q;
			Q.idp = l;
			l++;
			Q.nhat[0]    = sin(theta)*cos(phi);
			Q.nhat[1]    = sin(theta)*sin(phi);
			Q.nhat[2]    = cos(theta);
			Q.ipass      = ipass-1; // <- be careful, format of input file
			Q.av_gals.resize(MaxFiber); // <- added
			if (l%modulo == 0) {
				try {P.push_back(Q);} catch(std::exception& e) {myexception(e);}
			}
		}
	}
	fs.close();
	return(P);
}

// Assignment -----------------------------------------------------------------------------
Assignment::Assignment() {
	TF = initTable(MaxPlate,MaxFiber,-1);
	PG = initTable_pair(Npass,Ngal); // Doesn't work if defined directly
	na = 0;
}

Assignment::~Assignment() {}

// Assign g with tile/fiber (j,k), and check for duplicates
void Assignment::assign(int j, int k, int g, const Plates& P) {
	// Assign (j,k)
	int ipass = P[j].ipass;
	int q = TF[j][k];
	if (q != -1) printf("### !!! ### DUPLICATE (j,k) = (%d,%d) assigned with g = %d and %d ---> information on first g lost \n",j,k,q,g);
	TF[j][k] = g;
	// Assign (ipass,g)
	pair p = PG[ipass][g];
	if (!p.isnull()) printf("### !!! ### DUPLICATE (ipass,g) = (%d,%d) assigned with (j,k) = (%d,%d) and (%d,%d) ---> information on first (j,k) lost \n",ipass,g,p.f,p.s,j,k);
	PG[ipass][g] = pair(j,k);
	na++;
}

void Assignment::unassign(int j, int k, int g, const Plates& P) {
	TF[j][k] = -1;
	PG[P[j].ipass][g].setnull();
	na--;
}

bool Assignment::is_assigned_pg(int ip, int g) const {
	//printf(("is_assigned"+p2s(ip,g)+siz(PG)).c_str());
	pair p = PG[ip][g];
	return !p.isnull();
}

bool Assignment::is_assigned_tf(int j, int k) const {return (TF[j][k] != -1);}

std::vector<pair> Assignment::chosen_tfs(int g) const {
	std::vector<pair> chosen;
	for(int ip=0; ip<Npass; ip++) {
		pair tf = PG[ip][g];
		if (!tf.isnull()) chosen.push_back(tf);
	}
	return chosen;
}

List Assignment::unused_fibers() const {
	List unused = initList(MaxPlate);
	for(int j=0; j<MaxPlate; j++) {
		for (int k=0; k<MaxFiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j]++;
		}
	}
	return unused;
}

Table Assignment::unused_fibers_by_petal(const PP& pp) const {
	Table unused = initTable(MaxPlate,MaxPetal);
	List Sp = pp.spectrom;
	for(int j=0; j<MaxPlate; j++) {
		for (int k=0; k<MaxFiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j][Sp[k]]++;
		}
	}
	return unused;
}

List Assignment::hist_petal(int binsize_petal, const PP& pp) const {
	List hist = initList(100);
	Table unused = unused_fibers_by_petal(pp);
	for(int j=0;j<MaxPlate; j++){ 
		for(int l=0;l<MaxPetal; l++) 
			hist[unused[j][l]/binsize_petal]++;
	}
	return hist;
}

// Functions ------------------------------------------------------------------------------
// Counts how many time ob should be observed else more
int nobs(int g, const Gals& G, const Assignment& A) {
	int cnt(G[g].nobs);
	for (int i=0; i<Npass; i++) if (A.is_assigned_pg(i,g)) cnt--;
	return cnt;
}

inline double plate_dist(const double theta) {
	// Returns the radial distance on the plate (mm) given the angle,
	// theta (radians).  This is simply a fit to the data provided.
	const double p[4]={8.297e5,-1750.,1.394e4,0.0};
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

// Fast because no big function is called ! Don't change that
// For each ~5700 plate, finds ~25000 galaxies reachable by the plate,
// Projects them, and for each fiber, finds reachable ones
void collect_galaxies_for_all(const Gals& G, const htmTree<struct galaxy>& T, Plates& P, const PP& pp, Time &t) {
	List permut = random_permut(MaxPlate);
	int i;
	//omp_set_num_threads(24);
#pragma omp parallel
	{   int ipass=100; // ??
		int id = omp_get_thread_num();
		double quarter = floor(MaxPlate/(4*omp_get_num_threads()));
		double cnt, avg, std; int max = 0; int min = 1e7; 
		// Collects for each plate ; shuffle order of plates (useless !?)
		for (i=id; i<MaxPlate; i++) { // <- begins at id, otherwise all begin at 0 -> conflict
			int j = permut[i];
			plate p = P[j];
			// Takes neighboring ~25000 galaxies that can be reached by this plate
			std::vector<int> nbr = T.near(G,p.nhat,PlateRadius*M_PI/180.); // nbr for neighbours
			int n = nbr.size();
			// Projects thoses galaxies on the focal plane
			Onplates O;
			for (int k=0; k<n; k++) {
				int g = nbr[k];
				struct onplate op = change_coords(G[g],p); 
				op.id = g;
				O.push_back(op);
			}
			// For each fiber, finds all reachable galaxies
			KDtree<struct onplate> kdT(O,2);
			for (int k=0; k<MaxFiber; k++) {
				std::vector<int> gals = kdT.near(&(pp.fp[2*k]),PatrolMinRad,PatrolMaxRad);
				P[j].av_gals[k] = gals;
			}
			// Stats and avancement 
			cnt++; avg += n; std += n*n;
			if (n>max) max=n;
			if (n<min) min=n;
			if (cnt==quarter && id==0) {printf("  Thread 0 has done 1/4 of his job"); print_time(t," at");}
		}
		// Print stats and time
		avg /= cnt;
		std = sqrt(std/cnt-sq(avg));
		printf("  Thread %2d finished. Galaxies per plate : %7.1f +/- %5.1f [%5d,%5d] -",id,avg,std,min,max);
		print_time(t," at");
	} // End parallel
	// Print how many galaxies in range of a fiber
	List L;
	for (int j=0; j<MaxPlate; j++) {
		for (int k=0; k<MaxFiber; k++) {
			int s = P[j].av_gals[k].size();
			if (s>=L.size()) {L.resize(s+1); L[s] = 0;}
			L[s]++;
		}
	}
	print_list("\n How many galaxies in range of a fiber",L);

	printf("  TEST : Print 100 G[i].priority :\n");
	for (int i=0; i<100; i++) printf(" %d ",G[i].priority); // priority doesn't work
	printf("\n");
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

// Assign fibers naively
// assign_fibers is outside of parallel section because of interaction between plates
void assign_fibers(const Gals& G, const Plates& P, const PP& pp, Assignment& A) {
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
					List randGals = random_permut(av_gals.size()); // <-need to be different each time ?
					for (int n=0; n<av_gals.size(); n++) { 
						int g = av_gals[randGals[n]];
						int m = nobs(g,G,A);
						//forbid Ly-a in gray time, i.e. pass=5
						//forbid LRGs also 1/3/15 rnc
						if(!A.is_assigned_pg(ipass,g) && (m>0) &&((G[g].id==4)||(P[j].ipass!=5))) {
							if (G[g].priority<minp || (best>=0 && m>nobs(best,G,A))) { // !! G[-1] isn't initialized
								int kp = find_collision(j,k,g,pp,G,P,A);
								if (kp==-1) {
									minp = G[g].priority;
									best = g;
								}
							}
						}
					}
					// Assign best galaxy to this fiber
					if (best!=-1 && (nobs(best,G,A)>0)) A.assign(j,k,best,P);
					}}}}
	printf("  %s assignments\n",f(A.na));
}

bool ok_assign_g_to_jk(int g, int j, int k, const Plates& P, const Gals& G, const PP& pp, const Assignment& A) {
	// No collision, g not assigned at the ipass, (j,k) not assigned
	if (A.is_assigned_tf(j,k)) return false;
	if (A.is_assigned_pg(P[j].ipass,g)) return false;
	if (find_collision(j,k,g,pp,G,P,A)!=-1) return false;
	return true;
}

// Try to use unused fibers by reassigning some used ones
// Before : (jp,kp) <-> g ; (j,k) & gp free
// After : (j,k) <-> g & (jp,kp) <-> gp
// Prints the number of additional assigned galaxies
void improve(const Gals& G, const Plates&P, const PP& pp, Assignment& A) {
	int na_start = A.na;
	List randPlates = random_permut(MaxPlate);
	for (int jj=0; jj<MaxPlate; jj++) {
		int j = randPlates[jj];
		for(int k=0; k<MaxFiber; k++) {
			if (!A.is_assigned_tf(j,k)) { // Unused tilefiber (j,k)
				bool finished(false);
				std::vector<int> av_g = P[j].av_gals[k];
				for (int i=0; i<av_g.size() && !finished; i++) { // Not shuffled
					int g = av_g[i]; // g : possible galaxy for (j,k)
					// Is it allowed for jk to take g ?
					if (ok_assign_g_to_jk(g,j,k,P,G,pp,A)) {
						//what tilefibers have taken g?
						std::vector<pair> tfs = A.chosen_tfs(g);
						//could they take something else?
						for (int p=0; p<tfs.size() && !finished; p++){
							int jp = tfs[p].f;
							int kp = tfs[p].s; // (jp,kp) currently assigned to galaxy g
							std::vector<int> av_g2 = P[jp].av_gals[kp];
							for(int m=0; m<av_g2.size() && !finished; m++) {
								int gp = av_g2[m]; // gp : nother possibility for (jp,kp)
								bool ok = (nobs(gp,G,A)>0 && ok_assign_g_to_jk(gp,jp,kp,P,G,pp,A)); 
								if(!finished && ok) {
									// Modify assignment
									A.unassign(jp,kp,g,P);
									A.assign(j,k,g,P);
									A.assign(jp,kp,gp,P);
									finished = true; }}}}}}}}
	printf("  %s more assignments (%f %% improvement)\n",f(A.na-na_start),percent(A.na-na_start,na_start));
}

// Redistrubte assignments trying to get at least 500 free fibers on each plate/tile
// Redo so we look first at plates with too few free fibers
// Before : (j,k) <-> g, with Sp(k) too much used
// After : (jreassign,kreassign) <-> g & (j,k) free, such that (jreassign,kreassign) comes from most unused (ji,ki)
// Reassign if this is improvement
// The problem must be that we don't recompute unused_fbp at each reassignment 
void redistribute(const Gals& G, const Plates&P, const PP& pp, Assignment& A) {
	int redistributions(0);
	List Sp = pp.spectrom;
	Table unused_fbp = A.unused_fibers_by_petal(pp);
	//print_table("unused",unused_fbp);
	List randPlates = random_permut(MaxPlate);
	// Consider all petals with too few free fibers
	for (int jj=0;jj<MaxPlate;jj++) {
		int j = randPlates[jj];
		List randFibers = random_permut(MaxFiber);
		//for(int k=0;k<MaxFiber && (unused_fibers[j]<minUnused+1);++k){
		for (int kk=0; kk<MaxFiber; kk++) {
			int k = randFibers[kk];
			int nfree = unused_fbp[j][Sp[k]];
			int g = A.TF[j][k];
			if (g!=-1 && nfree<MinUnused/* && G[g].id==4*/) { // Only ELG
				// Consider other ways to observe this galaxy
				std::vector<pair> tfs = G[g].av_tfs;
				int jreassign(-1); int kreassign(-1); int mostunused(-1);
				for (int c=0; c<tfs.size(); c++) {
					int jp = tfs[c].f;
					int kp = tfs[c].s;
					int nfreep = unused_fbp[jp][Sp[kp]];
					// Use freest plate and petal
					if (nfreep>mostunused && ok_assign_g_to_jk(g,jp,kp,P,G,pp,A)) {
						mostunused = nfreep;
						jreassign = jp;
						kreassign = kp; 
					}
				}
				if (mostunused > nfree) {
					//printf("%d %d %d %d %d %d %d\n",mostunused, nfree, j,k,jreassign,kreassign,g);
					A.unassign(j,k,g,P);
					A.assign(jreassign,kreassign,g,P);
					unused_fbp[j][Sp[k]]++;
					unused_fbp[jreassign][Sp[kreassign]]--;		
					redistributions++; 
				}}}}
	printf("  %s redistributions (~%.4f %% redistributed)\n",f(redistributions),percent(redistributions,A.na));
}

// Detemine how many galaxies need to be dropped to guarantee 40 free fibers for each petal
// Count how many free fibers there are beyond 500 in each plate
void print_free_fibers(const PP& pp, const Assignment& A) {
	printf("# Free fibers statistics\n");
	// Create histogram with given binsize
	List hist_petal = A.hist_petal(BinSizePet,pp);
	//print_table("  Petals with this many free fiber : binsize ",hist_petal);
	printf("  Petals with this many free fiber : binsize = %d \n",BinSizePet);
	print_hist(hist_petal); // <- sure of what that does ?

	// Consider all petals with too few free fibers
	//List unused_fibers = A.unused_fibers();
	//int minUnused=500;
	Table unused_fibers_by_petal = A.unused_fibers_by_petal(pp);
	// Do not initialize like int a,b = 0; or there are bugs
	int counter(0); int beyond_counter(0); int npetal_upper(0); int npetal_under(0);
	for (int j=0; j<MaxPlate; j++) {
		for (int l=0; l<MaxPetal; l++) {
			int c = unused_fibers_by_petal[j][l]-MinUnused;
			if (c<0) { counter+= -c; npetal_under++; }
			else { beyond_counter+= c; npetal_upper++; }
		}
	}
	printf("  Number of free fibers under %d on each petal %s \n",MinUnused,f(counter));
	printf("  Number of free fibers beyond %d on each petal %s \n",MinUnused,f(beyond_counter));
	printf("  Number of petals which (number of free fibers < %d) : %s (%.4f %%)\n",MinUnused,f(npetal_under),percent(npetal_under,MaxPlate*MaxPetal));
	printf("  Number of petals which (number of free fibers >= %d) : %s (%.4f %%)\n",MinUnused,f(npetal_upper),percent(npetal_upper,MaxPlate*MaxPetal));
}

void display_results(const Gals& G, const List& goal, const Plates& P, const Assignment& A) {
	printf("# Results :\n");
	// Write out the results
	// hist2 not used (previous code)
	Table hist2 = initTable(GalaxyCategories+1,MaxObs+1);
	Table done = initTable(GalaxyCategories+1,MaxObs+1);
	List fibers_used = initList(GalaxyCategories+1);
	List targets = initList(GalaxyCategories+1);
	// Raw numbers of galaxies by id and number of remaining observations
	for (int g=0; g<Ngal; g++) {
		int n = nobs(g,G,A);
		if (n>=0 && n<=MaxObs) hist2[G[g].id][n]++;
	}
	print_table("  Remaining observations (id on lines)",hist2,true);

	for (int id=1; id<GalaxyCategories+1; id++) {
		List hist = initList(hist2[id]);
		for (int i=0; i<=MaxObs; i++) {
			//i here is number of observations lacking  i = goal -done
			//done=goal -i for i<=goal
			//0 for done>goal, i.e. i<0
			if (i<=goal[id]) {
				done[id][i]=hist[goal[id]-i];
				fibers_used[id]+=i*done[id][i];
				targets[id]+=done[id][i];
			}
		}
	}

	printf("\n By priority number, how many targets, how many reached\n");
	for (int pri=1; pri<=6; pri++) {
		int ntot=0,nass=0;
		for (int g=0; g<Ngal; g++) {
			if (fabs(G[g].priority-pri)<0.5) {
				if (nobs(g,G,A)==0) nass++;
				ntot++;
			}
		}
		printf("Priority %2d : %10s of %10s assigned fibers \n",pri,f(nass),f(ntot));
	}
	// Per sq deg in tex format
	double total_area(15789.);
	printf("\n tex format, per sq deg\n");
	printf("id| obsv'd  0         1          2          3          4          5        total fibers used    avail   obsvd \n");
	for (int id=1; id<=GalaxyCategories; id++) {
		printf("%2d&",id);
		for (int i=0; i<=goal[id]; i++) printf("%10.0f&",done[id][i]/total_area);
		for (int i=goal[id]+1; i<MaxObs+1; i++) printf("%11d&",0);
		float percentage = float(targets[id]-done[id][0])/float(targets[id]);
		printf("%10f&%10f&%10f&%10.4f \\ \n",(targets[id]-done[id][0])/total_area,fibers_used[id]/total_area,targets[id]/total_area,percentage);
	}
}

void plot_freefibers(std::string s, const Plates& P, const Assignment& A) {
	printf("# Plot free fibers positions\n");
	double pi = 3.1415926535;
	List unused_fibers = A.unused_fibers();
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
void time_count(int jmin, int jmax, const Gals& G, const List& goal, Table& nc, List& tss, const Assignment& A){
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
// nc should be only a list(galaxy_samples)

// time_count counted the number of times a galaxy was observed
void time_line(const Gals& G, const List& goal, const Assignment& A) {
	//produce numbers for python plot showing what fraction of galaxies of each sort are observed as a function of 'time'
	//what fraction of galaxies has been observed in each category as a function of "time"?
	//output for python plotting
	//10 spots for 

	//  1 = Ly-a QSO
	//  2 = QSO tracer
	//  3 = LRG
	//  4 = ELG
	//  5 = fake qsos
	//  6 = fake lrgs	

	// order the series Ly-a 1,2,3,4,5
	// QSO tracer 6
	// LRG 7,8
	// ELG 9
	int curves=galaxy_samples+4+1;//3 for Ly-a, 1 for LRG
	int nstep=MaxPlate/ntimes;
	Table nc = initTable(galaxy_samples+1,MaxObs+1); // new_count
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
		time_count(i*nstep,(i+1)*nstep,G,goal,nc,tss,A);
		for(int j=1; j<=galaxy_samples; ++j) {
			for (int k=1; k<=goal[j]; ++k) {
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

// Write some files -----------------------------------------------------------------------
// Write very large binary file of P[j].av_gals[k]
void writeTFfile(const Plates& P, std::ofstream TFfile) {
	TFfile.open ("TFout.txt", std::ios::binary);
	for(int j=0;j<MaxPlate;j++){
		for(int k=0;k<MaxFiber;k++){
			TFfile.write(reinterpret_cast<char*>(&j),sizeof(int));TFfile.write(reinterpret_cast<char*>(&k),sizeof(int));
			std::vector<int> gals = P[j].av_gals[k];
			int siz = gals.size();
			TFfile.write(reinterpret_cast<char*>(&siz),sizeof(int));
			for (int q=0;q<siz;q++){
				TFfile.write(reinterpret_cast<char*>(&gals[q]),sizeof(int));
			}

		}

	}
	TFfile.close();
}

void writeGfile(Gals& G, std::ofstream Gfile) {
	Gfile.open ("Gout.bn", std::ios::binary);
	std::cout<<G.size()<<std::endl;
	for(int i=0;i<G.size();i++){
		Gfile.write(reinterpret_cast<char*>(&G[i].ra),sizeof(double));Gfile.write(reinterpret_cast<char*>(&G[i].dec),sizeof(double));
		Gfile.write(reinterpret_cast<char*>(&G[i].z),sizeof(double));Gfile.write(reinterpret_cast<char*>(&G[i].nobs),sizeof(int));	
		Gfile.write(reinterpret_cast<char*>(&G[i].id),sizeof(int));Gfile.write(reinterpret_cast<char*>(&G[i].priority),sizeof(int));		
	}
	Gfile.close();
}

void readGfile(Gals& G, galaxy Gtemp, std::ifstream GXfile) {
	GXfile.open("Gout.bn", std::ios::binary);
	for(int i=0;!GXfile.eof();i++){
		GXfile.read(reinterpret_cast<char*>(&Gtemp.ra),sizeof(double));GXfile.read(reinterpret_cast<char*>(&Gtemp.dec),sizeof(double));
		GXfile.read(reinterpret_cast<char*>(&Gtemp.z),sizeof(double));GXfile.read(reinterpret_cast<char*>(&Gtemp.nobs),sizeof(int));	
		GXfile.read(reinterpret_cast<char*>(&Gtemp.id),sizeof(int));GXfile.read(reinterpret_cast<char*>(&Gtemp.priority),sizeof(int));
		G.push_back(Gtemp);
		if((i/1000000)*1000000==i) std::cout<<i<<std::endl;
	}
	GXfile.close();	
}
