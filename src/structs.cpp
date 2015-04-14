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
#include	"modules/htmTree.h"
#include	"modules/kdTree.h"
#include        "omp.h"
#include        "macros.h"
#include        "misc.h"
#include        "structs.h"
// Features ------------------------------------------------------------------
Feat::Feat() {
	prio.resize(Categories);
	goal.resize(Categories);
	kind.resize(Categories);
}

int Feat::id(str s) const {
	for (int i=0; i<Categories; i++) if (kind[i]==s) return i;
	std::cout << "ERROR in Feat id(), string not found in kind" << std::endl;
	return -1;
}

int Feat::maxgoal(int kind) const {
	int max(goal[kind]); int prio0= prio[kind];
	for (int i=0; i<Categories; i++) if (prio[i]==prio0 && goal[i]>max) max = goal[i];
	return max;
}

List Feat::maxgoal() const {
	List max = initList(Categories,-1);
	for (int i=0; i<Categories; i++) max[i] = maxgoal(i);
	return max;
}

void Feat::init_ids() {
	for (int i=0; i<Categories; i++) ids[kind[i]] = i;
}

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
Gals read_galaxies(const char fname[], int n) {
	Gals P;
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
		if (i%n == 0) {
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

int galaxy::prio(const Feat& F) const {
	return(F.prio[id]);
}

str galaxy::kind(const Feat& F) const {
	return(F.kind[id]);
}
// PP ----------------------------------------------------------------------------
// Read the positions of the fibers on each plate.  Assumes the format
// of fiberpos.txt produced by "randomize_fibers".
// need also to get the petal, i.e. spectrometer  rnc 1/16/15  added S
void PP::read_fiber_positions(const char pos_name[], int n) {
	str buf;
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
	int i(0);
	while (fs.eof()==0) {
		double x,y; int fiber,positioner,spectro,remove; 
		std::istringstream(buf) >> fiber >> positioner >> spectro >> x >> y;
		if (i%n == 0) {
		try{
			fp.push_back(x);
			fp.push_back(y);
			spectrom.push_back(spectro);  
		} catch(std::exception& e) {myexception(e);}
		}
		getline(fs,buf);
		i++;
	}
	fs.close();
}

PP::PP() {}

void PP::get_neighbors() {
	N.resize(Nfiber);
	for(int i=0; i<Nfiber; i++) {
		for (int j=0; j<Nfiber; j++) {
			if(i!=j) {
				if(sq(fp[2*i]-fp[2*j],fp[2*i+1]-fp[2*j+1]) < sq(NeighborRad)) {
					N[i].push_back(j); }}}}
}

void PP::compute_fibsofsp() {
	fibers_of_sp.resize(Npetal);
	for (int k=0; k<Nfiber; k++) fibers_of_sp[spectrom[k]].push_back(k);
}

List PP::fibs_of_same_pet(int k) const {
	return(fibers_of_sp[spectrom[k]]);
}
// plate ---------------------------------------------------------------------------
void plate::print_plate() const {
	printf("  Plate : %d - pass %d\n",idp,ipass); 
	printf("%f %f %f\n",nhat[0],nhat[1],nhat[2]);
	int r = rand() % Nfiber;
	print_list("Available galaxies of a random fiber :",av_gals[r]);
}

List plate::av_gals_plate() const {
	List gals = initList(Ngal);
	List L = initList(0);
	for (int k=0; k<Nfiber; k++) {
		for (int i=0; i<av_gals[k].size(); i++) {
			if (gals[av_gals[k][i]] == 0) {
				gals[av_gals[k][i]] = 1;
				L.push_back(i);
			}
		}
	}
	return L;
}

// Plates -----------------------------------------------------------------------------
// Read positions of the plate centers from an ascii file "center_name", and fill in a structure
// There is a version of this function (to adapt) for non ASCII files, ask to Robert Cahn
Plates read_plate_centers(const char center_name[], int modulo) {
	Plates P;
	str buf;
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
		str start;
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
			Q.av_gals.resize(Nfiber); // <- added
			if (l%modulo == 0) {
				try {P.push_back(Q);} catch(std::exception& e) {myexception(e);}
			}
		}
	}
	fs.close();
	return(P);
}

List gals_range_fibers(const Plates& P) {
	List L;
	for (int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			int s = P[j].av_gals[k].size();
			if (s>=L.size()) {L.resize(s+1); L[s] = 0;}
			L[s]++;
		}
	}
	return L;
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
	TF = initTable(Nplate,Nfiber,-1);
	GL = initPtable(Ngal,0); // Doesn't work if defined directly
	order.resize(Nplate);
	for (int i=0; i<Nplate; i++) order[i] = i;
	next_plate = 0;
	kinds = initCube(Nplate,Npetal,Categories);
	probas = initList(Categories);
	once_obs = initList(Ngal);
	nobsv = initList(Ngal);
	nobsv_tmp = initList(Ngal);
	List l = F.maxgoal();
	for (int g=0; g<Ngal; g++) {
		nobsv[g] = F.goal[G[g].id];
		nobsv_tmp[g] = l[G[g].id];
	}
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
}

void Assignment::verif(const Plates& P) const {
	for (int g=0; g<Ngal; g++) {
		Plist tfs = GL[g];
		int j0(-1); int j1(-1);
		for (int i=0; i<tfs.size(); i++) {
			pair tf = tfs[i];
			int j0 = j1;
			int j1 = tf.f;
			// Verif on TF
			if (TF[tf.f][tf.s]!=g) { printf("ERROR in verification of correspondance of galaxies !\n"); fl(); }
			// No 2 assignments within an interval of InterPlate
			if (j0!=-1 && fabs(j1-j0)<=InterPlate) { printf("ERROR in verification of interplate g=%d with j=%d and %d\n",g,j0,j1); fl(); }
		}
	}
	for (int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			int g = TF[j][k];
			// Verif on GL
			if (g!=-1 && isfound(pair(j,k),GL[g])==-1) { printf("ERROR in verification of correspondance of tfs !\n"); fl(); }
		}
	}
}

int Assignment::is_assigned_jg(int j, int g) const {
	for (int i=0; i<GL[g].size(); i++) if (GL[g][i].f == j) return i;
	return -1;
}

int Assignment::is_assigned_jg(int j, int g, int InterPlate) const {
	for (int i=0; i<GL[g].size(); i++) {
		if (max(j-InterPlate,0)<=GL[g][i].f && GL[g][i].f<=min(j+InterPlate,Nplate-1)) return i;
	}
	return -1;
}

bool Assignment::is_assigned_tf(int j, int k) const { return (TF[j][k] != -1); }

int Assignment::na(int begin, int size) const {
	int cnt(0);
	for (int j=begin; j<begin+size; j++) {
		for (int k=0; k<Nfiber; k++) {
			if (TF[j][k]!=-1) cnt++;
		}
	}
	return cnt;
}

Plist Assignment::chosen_tfs(int g, int begin, int size) const {
	Plist chosen;
	if (Nplate<begin+size) { printf("ERROR in chosen_tfs - size\n"); fl(); }
	for (int i=0; i<GL[g].size(); i++) {
		pair tf = GL[g][i];
		if (begin<=tf.f && tf.f<begin+size) {
			if (TF[tf.f][tf.s]!=g) { printf("ERROR in chosen_tfs\n"); fl(); }
			chosen.push_back(tf);
		}
	}
	return chosen;
}

Table Assignment::unused_fbp(const PP& pp) const {
	Table unused = initTable(Nplate,Npetal);
	List Sp = pp.spectrom;
	for(int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j][Sp[k]]++;
		}
	}
	return unused;
}

List Assignment::unused_f() const {
	List unused = initList(Nplate);
	for(int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j]++;
		}
	}
	return unused;
}

int Assignment::unused_f(int j) const {
	int unused(0);
	for (int k=0; k<Nfiber; k++) if (!is_assigned_tf(j,k)) unused++;
	return unused;
}

int Assignment::unused_fbp(int j, int k, const PP& pp) const {
	List fibs = pp.fibers_of_sp[pp.spectrom[k]];
	int unused(0);
	for (int i=0; i<fibs.size(); i++) {
		if (!is_assigned_tf(j,fibs[i])) unused++;
	}
	return unused;
}

Table Assignment::used_by_kind(str kind, const Gals& G, const PP& pp, const Feat& F) const {
	Table used = initTable(Nplate,Npetal);
	for(int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			int g = TF[j][k];
			if (g!=-1 && G[g].kind(F)==kind) used[j][pp.spectrom[k]]++;
		}
	}
	return used;
}

int Assignment::nkind(int j, int k, int kind, const Gals& G, const Plates& P, const PP& pp, const Feat& F, bool pet) const {
	if (!pet) return kinds[j][pp.spectrom[k]][kind];
	else return kinds[j][k][kind];
	//List fibers = pp.fibs_of_same_pet(k);
	//int cnt(0);
	//for (int i=0; i<fibers.size(); i++) {
		//int kk = fibers[i];
		//int g = TF[j][kk];
		//if (g!=-1 && G[g].kind(F)==kind) cnt++;
	//}
	//return cnt;
}

double Assignment::get_proba(int i, const Gals& G, const Feat& F) {
	int tot(0);
	int kind = F.prio[i];
	for (int l=0; l<Categories; l++) if (F.prio[l]==kind) tot += probas[l];
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
		for (int gg=0; gg<av_gals.size(); gg++) {
			int g = av_gals[gg];
			L.push_back(g==g0 ? -3 : -2);
			L.push_back(G[g].id);
			L.push_back(nobs(g,G,F));
			L.push_back(is_assigned_jg(j,g));
			L.push_back(is_assigned_jg(j,g,InterPlate));
			L.push_back(find_collision(j,k,g,pp,G,P));
		}
		T.push_back(L);
	}
	return T;
}

List Assignment::fibs_of_kind(int kind, int j, int pet, const Gals& G, const PP& pp, const Feat& F) const {
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<Nfbp; kk++) {
		int k = fibs[kk];
		int g = TF[j][k];
		if (g!=-1 && G[g].id==kind) L.push_back(k);
	}
	return L;
}

List Assignment::fibs_unassigned(int j, int pet, const Gals& G, const PP& pp, const Feat& F) const {
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<Nfbp; kk++) {
		int k = fibs[kk];
		if (!is_assigned_tf(j,k)) L.push_back(k);
	}
	return L;
}

int Assignment::nobs(int g, const Gals& G, const Feat& F, bool tmp) const {
	//int cnt(F.goal[G[g].id]);
	//for (int i=0; i<Npass; i++) if (is_assigned_pg(i,g)) cnt--;
	//return cnt;
	int obs = tmp ? nobsv_tmp[g] : nobsv[g]; // optimization
	return obs;
}

void Assignment::update_nobsv_tmp() {
	for (int g=0; g<Ngal; g++) if (once_obs[g]) nobsv_tmp[g] = nobsv[g];
}

void Assignment::update_nobsv_tmp_for_one(int j) {
	for (int k=0; k<Nfiber; k++) {
		int g = TF[j][k];
		if (g!=-1) nobsv_tmp[g] = nobsv[g];
	}
}

void Assignment::update_once_obs(int j) {
	for (int k=0; k<Nfiber; k++) {
		int g = TF[j][k];
		if (g!=-1) once_obs[g] = 1;
	}
}

// Useful sub-functions -------------------------------------------------------------------------------------------------
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

// (On plate p) finds if there is a collision if fiber k would watch at galaxy g (collision with neighb)
int Assignment::find_collision(int j, int k, int g, const PP& pp, const Gals& G, const Plates& P) const {
	struct onplate op = change_coords(G[g],P[j]);
	double x = op.pos[0];
	double y = op.pos[1];
	for (int i=0; i<pp.N[k].size(); i++) {
		int kn = pp.N[k][i];
		int gn = TF[j][kn];
		if (gn!=-1) {
			struct onplate opn = change_coords(G[gn],P[j]);
			double xn = opn.pos[0];
			double yn = opn.pos[1];
			if (sq(x-xn,y-yn) < sq(Collide)) return kn;
		}
	}
	return -1;
}

// Write some files -----------------------------------------------------------------------
// Write very large binary file of P[j].av_gals[k]. Not checked for a long time, should not work anymore !
void writeTFfile(const Plates& P, std::ofstream TFfile) {
	TFfile.open ("TFout.txt", std::ios::binary);
	for(int j=0;j<Nplate;j++){
		for(int k=0;k<Nfiber;k++){
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
		Gfile.write(reinterpret_cast<char*>(&G[i].z),sizeof(double));
		Gfile.write(reinterpret_cast<char*>(&G[i].id),sizeof(int));
	}
	Gfile.close();
}

void readGfile(Gals& G, galaxy Gtemp, std::ifstream GXfile) {
	GXfile.open("Gout.bn", std::ios::binary);
	for(int i=0;!GXfile.eof();i++){
		GXfile.read(reinterpret_cast<char*>(&Gtemp.ra),sizeof(double));GXfile.read(reinterpret_cast<char*>(&Gtemp.dec),sizeof(double));
		GXfile.read(reinterpret_cast<char*>(&Gtemp.z),sizeof(double));
		GXfile.read(reinterpret_cast<char*>(&Gtemp.id),sizeof(int));
		G.push_back(Gtemp);
		if((i/1000000)*1000000==i) std::cout<<i<<std::endl;
	}
	GXfile.close();	
}
