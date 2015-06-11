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
#include        <map>
#include        <stdlib.h>     /* srand, rand */
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "collision.h"

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
// Read galaxies from binary file--format is ra, dec, z, priority and nobs
// with ra/dec in degrees. Priority and nobs information are treated diffenrently now, by Feat, and should be removed here.
//  Reads every n galaxies (to test quicker)
Gals read_galaxies(const Feat& F) {
	Gals P;
	const char* fname;
	fname = F.galFile.c_str();
	std::ifstream fs(fname,std::ios::binary);
	if (!fs) {  // An error occurred opening the file.
		std::cerr << "Unable to open file " << fname << std::endl;
		myexit(1);
	}
	int Nobj;
	fs.read((char *)&Nobj,sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl; myexit(1);}
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
	fs.read((char *)&pr[0],Nobj*sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.read((char *)&no[0],Nobj*sizeof(int));
	if (fs.fail()) {std::cerr<<"Error reading "<<fname<<std::endl;myexit(1);}
	fs.close();
	// reduce time
	for (int i=0; i<Nobj; i++) {
		if (i%F.moduloGal == 0) {
		double theta = (90 - (double)dc[i])*M_PI/180;
		double phi   = ((double)ra[i]     )*M_PI/180;
		struct galaxy Q;
		Q.id         = id[i]-1; // -1 added
		Q.nhat[0]    = sin(theta)*cos(phi);
		Q.nhat[1]    = sin(theta)*sin(phi);
		Q.nhat[2]    = cos(theta);
		Q.ra = ra[i];
		Q.dec = dc[i];
		Q.z = zz[i];
		try{P.push_back(Q);}catch(std::exception& e) {myexception(e);}
		}
	}
	return P;
}


str galaxy::kind(const Feat& F) const {
	return(F.kind[id]);
}

// PP ----------------------------------------------------------------------------
// Read the positions of the fibers on each plate.  Assumes the format
// of fiberpos.txt produced by "F.Randomize_fibers".
// need also to get the petal, i.e. spectrometer  rnc 1/16/15  added S
void PP::read_fiber_positions(const Feat& F) {
	str buf;
	std::ifstream fs(F.fibFile.c_str());

	if (!fs) { // An error occurred opening the file.
		std::cerr << "Unable to open file " << F.fibFile << std::endl;
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
	int i(0);
	int petals_pac[] = {0,1,2,7,8,9};
	List petals_pacL = initList(petals_pac,6);
	List inv = inverse(petals_pacL);
	while (fs.eof()==0) {
		double x,y; int fiber,positioner,spectro,remove; 
		std::istringstream(buf) >> fiber >> positioner >> spectro >> x >> y;
		if (i%F.moduloFiber == 0) {
		if (!F.Pacman || isfound(spectro,petals_pacL)) {
		try{
			fp.push_back(x);
			fp.push_back(y);
			int sp = F.Pacman ? inv[spectro] : spectro;
			spectrom.push_back(sp);  
		} catch(std::exception& e) {myexception(e);}
		}
		}
		getline(fs,buf);
		i++;
	}
	fs.close();
}

PP::PP() {}

void PP::get_neighbors(const Feat& F) {
	N.resize(F.Nfiber);
	for(int i=0; i<F.Nfiber; i++) {
		for (int j=0; j<F.Nfiber; j++) {
			if(i!=j) {
				if(sq(fp[2*i]-fp[2*j],fp[2*i+1]-fp[2*j+1]) < sq(F.NeighborRad)) {
					N[i].push_back(j); }}}}
}

void PP::compute_fibsofsp(const Feat& F) {
	fibers_of_sp.resize(F.Npetal);
	for (int k=0; k<F.Nfiber; k++) fibers_of_sp[spectrom[k]].push_back(k);
}

List PP::fibs_of_same_pet(int k) const {
	return(fibers_of_sp[spectrom[k]]);
}

dpair PP::coords(int k) const {
	return dpair(fp[2*k],fp[2*k+1]);
}

// plate ---------------------------------------------------------------------------
void plate::print_plate(const Feat& F) const {
	printf("  Plate : %d - pass %d\n",idp,ipass); 
	printf("%f %f %f\n",nhat[0],nhat[1],nhat[2]);
	int r = rand() % F.Nfiber;
	print_list("Available galaxies of a random fiber :",av_gals[r]);
}

List plate::av_gals_plate(const Feat& F) const {
	List gals = initList(F.Ngal);
	List L = initList(0);
	for (int k=0; k<F.Nfiber; k++) {
		for (int i=0; i<av_gals[k].size(); i++) {
			if (gals[av_gals[k][i]] == 0) {
				gals[av_gals[k][i]] = 1;
				L.push_back(i);
			}
		}
	}
	return L;
}


// Plates ---------------------------------------------------------------------------
// Read positions of the plate centers from an ascii file "center_name", and fill in a structure
// There is a version of this function (to adapt) for non ASCII files, ask to Robert Cahn
Plates read_plate_centers(const Feat& F) {
	Plates P;
	str buf;
	std::ifstream fs(F.tileFile.c_str());
	if (!fs) {  // An error occurred opening the file.
		std::cerr << "Unable to open file " << F.tileFile << std::endl;
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
		str start;
		ss>> start;
		ss>>tileid>>ra >> dec >> ipass>>in_desi>>ebv>>airmass>>exposefac;	
		//only keep those in footprint, fix ra to lie between 0 and 360
		if (in_desi==1) {
			if (ra<   0.) {ra += 360.;}
			if (ra>=360.) {ra -= 360.;}
			if (dec<-90. || dec>90.) {
				std::cout << "DEC="<<dec<<" out of range reading " << F.tileFile<<std::endl;
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
			Q.av_gals.resize(F.Nfiber); // <- added
			Q.density.resize(F.Nfiber); // <- added
			try {P.push_back(Q);} catch(std::exception& e) {myexception(e);}
		}
	}
	fs.close();
	return(P);
}

List av_gals_of_kind(int kind, int j, int k, const Gals& G, const Plates& P, const Feat& F) {
	List L;
	List av_gals = P[j].av_gals[k];
	for (int gg=0; gg<av_gals.size(); gg++) {
		int g = av_gals[gg];
		if (G[g].id==kind) L.push_back(g);
	}
	return L;
}

// Assignment -----------------------------------------------------------------------------
Assignment::Assignment(const Gals& G, const Feat& F) {
	TF = initTable(F.Nplate,F.Nfiber,-1);
	GL = initPtable(F.Ngal,0); // Doesn't work if defined directly
	order.resize(F.Nplate);
	for (int i=0; i<F.Nplate; i++) order[i] = i;
	next_plate = 0;
	kinds = initCube(F.Nplate,F.Npetal,F.Categories);
	probas = initList(F.Categories);
	once_obs = initList(F.Ngal);
	nobsv = initList(F.Ngal);
	nobsv_tmp = initList(F.Ngal);
	List l = F.maxgoal();
	for (int g=0; g<F.Ngal; g++) {
		nobsv[g] = F.goal[G[g].id];
		nobsv_tmp[g] = l[G[g].id];
	}
	unused = initTable(F.Nplate,F.Npetal,F.Nfbp);
}

Assignment::~Assignment() {}

// Assign g with tile/fiber (j,k), and check for duplicates
void Assignment::assign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp) {
	// Assign (j,k)
	int q = TF[j][k];
	if (q != -1) {
		printf("### !!! ### DUPLICATE (j,k) = (%d,%d) assigned with g = %d and %d ---> information on first g lost \n",j,k,q,g);
		myexit(1);
	}
	TF[j][k] = g;
	// Assign g
	Plist pl = GL[g];
	pair p = pair(j,k);
	int a = isfound(p,pl); // Can be erased for optimization once there is no duplicate
	if (a!=-1) {
		printf("### !!! ### DUPLICATE g = %d assigned with (j,k) = (%d,%d) and (%d,%d) ---> information on first (j,k) lost \n",g,pl[a].f,pl[a].s,j,k);
		erase(a,GL[g]);
		myexit(1); // Can be commented if want to force continuing
	}
	GL[g].push_back(p);
	// Kinds
	kinds[j][pp.spectrom[k]][G[g].id]++;
	// Probas
	probas[G[g].id]++;
	// Nobsv
	nobsv[g]--;
	nobsv_tmp[g]--;

	unused[j][pp.spectrom[k]]--;
}

void Assignment::unassign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp) {
	if (TF[j][k]==-1) printf("### !!! ### TF (j,k) = (%d,%d) gets unassigned but was already not assigned\n",j,k);
	int a = isfound(pair(j,k),GL[g]);
	if (a==-1) printf("### !!! ### Galaxy g = %d gets unassigned but was already not assigned\n",g);

	TF[j][k] = -1;
	if (a!=-1) erase(a,GL[g]);
	kinds[j][pp.spectrom[k]][G[g].id]--;
	probas[G[g].id]--;
	nobsv[g]++;
	nobsv_tmp[g]++;

	unused[j][pp.spectrom[k]]++;
}

void Assignment::verif(const Plates& P, const Gals& G, const PP& pp, const Feat& F) const {
	str qso_lrgA[] = {"QSOLy-a","QSOTracer","FakeQSO","LRG","FakeLRG"}; List qso_lrg = F.init_ids_list(qso_lrgA,5);
	for (int g=0; g<F.Ngal; g++) {
		Plist tfs = GL[g];
		int j0(-1); int j1(-1);
		for (int i=0; i<tfs.size(); i++) {
			pair tf = tfs[i];
			int j0 = j1;
			int j1 = tf.f;
			// Verif on TF
			if (TF[tf.f][tf.s]!=g) { printf("ERROR in verification of correspondance of galaxies !\n"); fl(); }
			// No 2 assignments within an interval of F.InterPlate
			if (j0!=-1 && isfound(G[g].id,qso_lrg) && fabs(j1-j0)<F.InterPlate) { printf("ERROR in verification of F.InterPlate g=%d with j=%d and %d\n",g,j0,j1); fl(); }
		}
	}
	for (int j=0; j<F.Nplate; j++) {
		List gals = initList(F.Ngal);
		for (int k=0; k<F.Nfiber; k++) {
			int g = TF[j][k];
			if (g!=-1) {
				// Verif on GL
				if (isfound(pair(j,k),GL[g])==-1) { printf("ERROR in verification of correspondance of tfs !\n"); fl(); }
				// Verif that a galaxy isn't observed twice
				if (gals[g]==1) printf("ERROR in verification, twice the same galaxy by (%d,%d)\n",j,k);
				else gals[g] = 1;
				// Collision checking
				if (!F.Collision && is_collision(j,k,pp,G,P,F)!=-1) printf("ERROR in verification : collisions\n");
			}
		}
	}
	Table usedSS = used_by_kind("SS",G,pp,F);
	Table usedSF = used_by_kind("SF",G,pp,F);
	for (int j=0; j<F.Nplate; j++) {
		for (int n=0; n<F.Npetal; n++) {
			if (usedSS[j][n]!=F.MaxSS || usedSF[j][n]!=F.MaxSF) printf("ERROR in verification : number of SF or SS\n");
		}
	}
}

int Assignment::is_assigned_jg(int j, int g) const {
	for (int i=0; i<GL[g].size(); i++) if (GL[g][i].f == j) return i;
	return -1;
}

int Assignment::is_assigned_jg(int j, int g, const Gals& G, const Feat& F) const {
	for (int i=0; i<GL[g].size(); i++) if ((/*(F.iftype(G[g].id,"LRG") || F.iftype(G[g].id,"QSO")) &&*/ fabs(j-GL[g][i].f)<F.InterPlate) || j==i) return i; // Makes it crash
	return -1;
}

bool Assignment::is_assigned_tf(int j, int k) const { return (TF[j][k] != -1); }

int Assignment::na(const Feat& F, int begin, int size) const {
	int size1 = (size==-1) ? F.Nplate : size;
	int cnt(0);
	for (int j=begin; j<begin+size1; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (TF[j][k]!=-1) cnt++;
		}
	}
	return cnt;
}

Plist Assignment::chosen_tfs(int g, const Feat& F, int begin, int size0) const {
	int size = (size0==-1) ? F.Nplate : size0;
	Plist chosen;
	if (F.Nplate<begin+size) { printf("ERROR in chosen_tfs - size\n"); fl(); }
	for (int i=0; i<GL[g].size(); i++) {
		pair tf = GL[g][i];
		if (begin<=tf.f && tf.f<begin+size) {
			if (TF[tf.f][tf.s]!=g) { printf("ERROR in chosen_tfs\n"); fl(); }
			chosen.push_back(tf);
		}
	}
	return chosen;
}

Table Assignment::unused_fbp(const PP& pp, const Feat& F) const {
	Table unused = initTable(F.Nplate,F.Npetal);
	List Sp = pp.spectrom;
	for(int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j][Sp[k]]++;
		}
	}
	return unused;
}

List Assignment::unused_f(const Feat& F) const {
	List unused = initList(F.Nplate);
	for(int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j]++;
		}
	}
	return unused;
}

int Assignment::unused_f(int j, const Feat& F) const {
	int unused(0);
	for (int k=0; k<F.Nfiber; k++) if (!is_assigned_tf(j,k)) unused++;
	return unused;
}

int Assignment::unused_fbp(int j, int k, const PP& pp, const Feat& F) const {
	List fibs = pp.fibers_of_sp[pp.spectrom[k]];
	int unused(0);
	for (int i=0; i<fibs.size(); i++) {
		if (!is_assigned_tf(j,fibs[i])) unused++;
	}
	return unused;
}

Table Assignment::used_by_kind(str kind, const Gals& G, const PP& pp, const Feat& F) const {
	Table used = initTable(F.Nplate,F.Npetal);
	for(int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			int g = TF[j][k];
			if (g!=-1 && G[g].kind(F)==kind) used[j][pp.spectrom[k]]++;
		}
	}
	return used;
}

int Assignment::nkind(int j, int k, int kind, const Gals& G, const Plates& P, const PP& pp, const Feat& F, bool pet) const {
	if (!pet) return kinds[j][pp.spectrom[k]][kind];
	else return kinds[j][k][kind];
}

double Assignment::get_proba(int i, const Gals& G, const Feat& F) {
	int tot(0);
	int kind = F.prio[i];
	for (int l=0; l<F.Categories; l++) if (F.prio[l]==kind) tot += probas[l];
	return ((double) probas[i])/((double) tot);
}
	
Table Assignment::infos_petal(int j, int pet, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const {
	Table T;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<fibs.size(); kk++) {
		List L;
		int k = fibs[kk];
		int g0 = TF[j][k];
		L.push_back(g0==-1 ? -1 : G[g0].id);
		List av_gals = P[j].av_gals[k];
		bool lya(false);
		for (int gg=0; gg<av_gals.size(); gg++) {
			int g = av_gals[gg];
			L.push_back(g==g0 ? -3 : -2);
			L.push_back(G[g].id);
			L.push_back(nobs(g,G,F));
			L.push_back(is_assigned_jg(j,g));
			L.push_back(is_assigned_jg(j,g,G,F));
			L.push_back(find_collision(j,k,g,pp,G,P,F));
			if (G[g].id==0) lya = true;
		}
		if (lya) T.push_back(L);
	}
	return T;
}

List Assignment::fibs_of_kind(int kind, int j, int pet, const Gals& G, const PP& pp, const Feat& F) const {
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<F.Nfbp; kk++) {
		int k = fibs[kk];
		int g = TF[j][k];
		if (g!=-1 && G[g].id==kind) L.push_back(k);
	}
	return L;
}

int num_av_gals(int j, int k, const Gals& G, const Plates& P, const Feat& F, const Assignment& A) {
	int cnt = 0;
	for (int i=0; i<P[j].av_gals[k].size(); i++) {
		int g = P[j].av_gals[k][i];
		int id = G[g].id;
		if (isfound(id,F.no_ss_sf)) cnt += A.nobs_time(g,j,G,F);
	}
	return cnt;
}
/*
List Assignment::fibs_of_kind_sorted(int kind, int j, int pet, const Gals& G, const PP& pp, const Feat& F) const {
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<F.Nfbp; kk++) {
		int k = fibs[kk];
		int g = TF[j][k];
		if (g!=-1 && G[g].id==kind) L.push_back(k);
	}
	List num_av_gals;
	for (int k=0; k<L.size(); k++) {
num_av_gals.push_back(P[j].num_av_gals(k,);
A.nobs_time(gg,j,G,F)
	std::sort(L.begin(), L.end());  

}
}
*/

List Assignment::sort_fibs_dens(int j, const List& fibs, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const {
	List num;
	for (int k=0; k<fibs.size(); k++) num.push_back(P[j].density[fibs[k]]);
	List perm = get_permut_sort(num);
	List fibs_sorted;
	for (int k=0; k<num.size(); k++) fibs_sorted.push_back(fibs[perm[k]]);
	return fibs_sorted;
}

List Assignment::fibs_unassigned(int j, int pet, const Gals& G, const PP& pp, const Feat& F) const {
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<F.Nfbp; kk++) {
		int k = fibs[kk];
		if (!is_assigned_tf(j,k)) L.push_back(k);
	}
	return L;
}

int Assignment::nobs(int g, const Gals& G, const Feat& F, bool tmp) const {
	int obs = tmp ? nobsv_tmp[g] : nobsv[g]; // optimization
	return obs;
}

void Assignment::update_nobsv_tmp(const Feat& F) {
	for (int g=0; g<F.Ngal; g++) if (once_obs[g]) nobsv_tmp[g] = nobsv[g];
}

void Assignment::update_nobsv_tmp_for_one(int j, const Feat& F) {
	for (int k=0; k<F.Nfiber; k++) {
		int g = TF[j][k];
		if (g!=-1) nobsv_tmp[g] = nobsv[g];
	}
}

void Assignment::update_once_obs(int j, const Feat& F) {
	for (int k=0; k<F.Nfiber; k++) {
		int g = TF[j][k];
		if (g!=-1) once_obs[g] = 1;
	}
}

int Assignment::nobs_time(int g, int j, const Gals& G, const Feat& F) const {
	int kind = G[g].id;
	int cnt = once_obs[g] ? F.goal[kind] : F.maxgoal(kind);
	for (int i=0; i<GL[g].size(); i++) if (GL[g][i].f<j) cnt--;
	return cnt;
}

// Useful sub-functions -------------------------------------------------------------------------------------------------

int fprio(int g, const Gals& G, const Feat& F, const Assignment& A) {
	if (A.once_obs[g]) return F.priopost[G[g].id];
	else return F.prio[G[g].id];
}

// Returns the radial distance on the plate (mm) given the angle,
// theta (radians).  This is simply a fit to the data provided.
double plate_dist(const double theta) {
	const double p[4] = {8.297e5,-1750.,1.394e4,0.0};
	double rr=0;
	for (int i=0; i<4; i++) rr = theta*rr + p[i];
	return rr;
}

// Returns the x-y position on the plate centered at P for galaxy O.
struct onplate change_coords(const struct galaxy& O, const struct plate& P) {
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
	return obj; 
}


bool collision(dpair O1, dpair G1, dpair O2, dpair G2, const Feat& F) {
	double dist_sq = sq(G1,G2);
	if (dist_sq < sq(F.Collide)) return true;
	if (dist_sq > sq(F.NoCollide)) return false;
	PosP posp(3,3);
	polygon fh1 = F.fh;
	polygon fh2 = F.fh;
	polygon cb1 = F.cb;
	polygon cb2 = F.cb;
	repos_cb_fh(cb1,fh1,O1,G1,posp);
	repos_cb_fh(cb2,fh2,O2,G2,posp);
	if (collision(fh1,fh2)) return true;
	if (collision(cb1,fh2)) return true;
	if (collision(cb2,fh1)) return true;
	return false;
}

// (On plate p) finds if there is a collision if fiber k would watch at galaxy g (collision with neighb)
int Assignment::find_collision(int j, int k, int g, const PP& pp, const Gals& G, const Plates& P, const Feat& F, int col) const {
	bool bol = (col==-1) ? F.Collision : false;
	if (bol) return -1;
	dpair G1 = projection(g,j,G,P);
	for (int i=0; i<pp.N[k].size(); i++) {
		int kn = pp.N[k][i];
		int gn = TF[j][kn];
		if (gn!=-1) {
			dpair G2 = projection(gn,j,G,P);
			bool b = F.Exact ? collision(pp.coords(k),G1,pp.coords(kn),G2,F) : (sq(G1,G2) < sq(F.AvCollide));
			if (b) return kn;
		}
	}
	return -1;
}

bool Assignment::find_collision(int j, int k, int kn, int g, int gn, const PP& pp, const Gals& G, const Plates& P, const Feat& F, int col) const {
	bool bol = (col==-1) ? F.Collision : false;
	if (bol) return false;
	dpair G1 = projection(g,j,G,P);
	dpair G2 = projection(gn,j,G,P);
	return F.Exact ? collision(pp.coords(k),G1,pp.coords(kn),G2,F) : (sq(G1,G2) < sq(F.AvCollide));
}

int Assignment::is_collision(int j, int k, const PP& pp, const Gals& G, const Plates& P, const Feat& F) const {
	int g = TF[j][k];
	if (g!=-1) return find_collision(j,k,g,pp,G,P,F,0);
	else return -1;
}

float Assignment::colrate(const PP& pp, const Gals& G, const Plates& P, const Feat& F, int jend0) const {
	int jend = (jend0==-1) ? F.Nplate : jend0;
	int col = 0;
	for (int j=0; j<jend; j++) {
		List done = initList(F.Nfiber);
		for (int k=0; k<F.Nfiber; k++) {
			if (done[k] == 0) {
				int c = is_collision(j,k,pp,G,P,F);
				if (c!=-1) {
					done[c] = 1;
					col += 2;
				}
			}
		}
	}
	return percent(col,jend*F.Nfiber);
}

dpair projection(int g, int j, const Gals& G, const Plates& P) {
	struct onplate op = change_coords(G[g],P[j]);
	return dpair(op.pos[0],op.pos[1]);
}


pyplot::pyplot(polygon p) {
	pol = p;
}

void pyplot::addtext(dpair p, str s) {
	text.push_back(s);
	textpos.push_back(p);
}

void pyplot::plot_tile(str directory, int j, const Feat& F) const {
	FILE * file;
	str fname = directory+"/tile"+i2s(j)+".py";
	file = fopen(fname.c_str(),"w");

	// Header
	Dlist lims = pol.limits();
	fprintf(file,"from pylab import *\nimport pylab as pl\nimport matplotlib.pyplot as plt\nfrom matplotlib import collections as mc\nax=subplot(aspect='equal')\naxes = plt.gca()\naxes.set_xlim([%f,%f])\naxes.set_ylim([%f,%f])\nax.axis('off')\nax.get_xaxis().set_visible(False)\nax.get_yaxis().set_visible(False)\nset_cmap('hot')\nfig = plt.gcf()\n\n",lims[0],lims[1],lims[2],lims[3]);
	if (j!=-1) fprintf(file,"plt.text(350,-350,'Tile %d',horizontalalignment='center',verticalalignment='center',size=5)\n\n",j);

	// Plot polygon
	for (int i=0; i<pol.elmts.size(); i++) {
		element e = pol.elmts[i];
		if (e.is_seg) {
			if (1<e.segs.size()) {
				fprintf(file,"lines = [[");
				for (int j=0; j<e.segs.size(); j++) fprintf(file,"(%f,%f),",e.segs[j].f,e.segs[j].s);
				if (e.color!='w') fprintf(file,"]]\nlc = mc.LineCollection(lines,linewidths=0.2,color='%c')\nax.add_collection(lc)\n",e.color);
				else fprintf(file,"]]\nlc = mc.LineCollection(lines,linewidths=0.2,color='k',alpha=0.4)\nax.add_collection(lc)\n");
			}
			if (e.segs.size()==1) fprintf(file,"circ=plt.Circle((%f,%f),%f,fill=True,linewidth=0.1,alpha=%f,edgecolor='none',fc='%c')\nfig.gca().add_artist(circ)\n",e.segs[0].f,e.segs[0].s,e.radplot,e.transparency,e.color);
		}
		else {
			if (e.color!='w') fprintf(file,"circ=plt.Circle((%f,%f),%f,fill=False,linewidth=0.2,color='%c')\nfig.gca().add_artist(circ)\n",e.O.f,e.O.s,e.rad,e.color);
			else fprintf(file,"circ=plt.Circle((%f,%f),%f,fill=False,linewidth=0.2,color='k',alpha=0.4)\nfig.gca().add_artist(circ)\n",e.O.f,e.O.s,e.rad);
		}
	}

	// Plot text
	if (text.size()!=textpos.size()) printf("Error sizes pyplot text\n");
	else for (int i=0; i<text.size(); i++) fprintf(file,"plt.text(%f,%f,'%s',horizontalalignment='center',verticalalignment='center',size=1)\n\n",textpos[i].f,textpos[i].s,text[i].c_str());
	
	// Finally
	fprintf(file,"\nfig.savefig('tile%d.pdf',bbox_inches='tight',pad_inches=0,dpi=(300))",j);
	fclose(file);
}
