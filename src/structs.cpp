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
void PP::read_fiber_positions(const char pos_name[]) {
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
	printf("  Plate : %d - pass %d\n",idp,ipass); printf("%f %f %f\n",nhat[0],nhat[1],nhat[2]);
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
	int pass(0);
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
			if (ipass-1<pass) printf("ERROR reading plate centers : passes are not in an inscreasing order !\n");
			pass = ipass-1;
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
// Assignment -----------------------------------------------------------------------------
Assignment::Assignment() {
	TF = initTable(Nplate,Nfiber,-1);
	PG = initTable_pair(Npass,Ngal); // Doesn't work if defined directly
	kinds = initCube(Nplate,Npetal,Categories);
	probas = initDlist(Categories);
}

Assignment::~Assignment() {}

// Assign g with tile/fiber (j,k), and check for duplicates
void Assignment::assign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp) {
	// Assign (j,k)
	int ipass = P[j].ipass;
	int q = TF[j][k];
	if (q != -1) printf("### !!! ### DUPLICATE (j,k) = (%d,%d) assigned with g = %d and %d ---> information on first g lost \n",j,k,q,g);
	TF[j][k] = g;
	// Assign (ipass,g)
	pair p = PG[ipass][g];
	if (!p.isnull()) printf("### !!! ### DUPLICATE (ipass,g) = (%d,%d) assigned with (j,k) = (%d,%d) and (%d,%d) ---> information on first (j,k) lost \n",ipass,g,p.f,p.s,j,k);
	PG[ipass][g] = pair(j,k);
	// Kinds
	kinds[j][pp.spectrom[k]][G[g].id]++;
}

void Assignment::unassign(int j, int k, int g, const Gals& G, const Plates& P, const PP& pp) {
	TF[j][k] = -1;
	PG[P[j].ipass][g].setnull();
	kinds[j][pp.spectrom[k]][G[g].id]--;
}

bool Assignment::is_assigned_pg(int ip, int g) const {
	//printf(("is_assigned"+p2s(ip,g)+siz(PG)).c_str());
	pair p = PG[ip][g];
	return (!p.isnull());
}

bool Assignment::is_assigned_tf(int j, int k) const {return (TF[j][k] != -1);}

int Assignment::na() {
	int cnt(0);
	for (int j=0; j<Nplate; j++) {
		for (int k=0; k<Nfiber; k++) {
			if (TF[j][k]!=-1) cnt++;
		}
	}
	return cnt;
}

std::vector<pair> Assignment::chosen_tfs(int g) const {
	std::vector<pair> chosen;
	for(int ip=0; ip<Npass; ip++) {
		pair tf = PG[ip][g];
		if (!tf.isnull()) chosen.push_back(tf);
	}
	return chosen;
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

int Assignment::unused_f(int j) {
	int unused(0);
	for (int k=0; k<Nfiber; k++) 	if (!is_assigned_tf(j,k)) unused++;
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

bool Assignment::once_obs(int g) {
	for (int i=0; i<Npass; i++) if (!PG[i][g].isnull()) return true;
	return false;
}

int Assignment::nkind(int j, int k, str kind, const Gals& G, const Plates& P, const PP& pp, const Feat& F) const {
	return kinds[j][pp.spectrom[k]][F.id(kind)];
	//List fibers = pp.fibs_of_same_pet(k);
	//int cnt(0);
	//for (int i=0; i<fibers.size(); i++) {
		//int kk = fibers[i];
		//int g = TF[j][kk];
		//if (g!=-1 && G[g].kind(F)==kind) cnt++;
	//}
	//return cnt;
}

void Assignment::update_probas(const Gals& G, const Feat& F) {
	List L = initList(Categories);
	int np = plates_done.size();
	for (int j=0; j<np; j++) {
		for (int k=0; k<Nfiber; k++) {
			int g = TF[j][k];
			if (g!=-1) L[G[g].id]++;
		}
	}
	for (int i=0; i<Categories; i++) {
		int tot(0);
		int kind = F.prio[i];
		for (int l=0; l<Categories; l++) if (F.prio[l]==kind) tot += L[l];
		probas[i] = ((double) L[i])/((double) tot);
	}
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
