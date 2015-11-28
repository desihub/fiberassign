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

// Read galaxies from binary file--format is ra, dec, z, priority and nobs
// with ra/dec in degrees.
//  Reads every n galaxies (to test more quickly)
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
            //if(Q.dec<F.MaxDec && Q.dec>F.MinDec &&Q.ra<F.MaxRa && Q.ra>F.MinRa){
                try{P.push_back(Q);}catch(std::exception& e) {myexception(e);}
            //}
		}
	}
	return P;
}

std::vector<int> count_galaxies(const Gals& G){
    std::vector <int> counter(10,0);
    for (int i=0;i<G.size();i++){
        counter[G[i].id]+=1;
                  }
    return counter;
}


Gals read_galaxies_ascii(const Feat& F)
// Read objects from ascii file--format is ra, dec, priority and nobs
// with ra/dec in degrees.
{
    Gals P;
    std::string buf;
    const char* fname;
    fname= F.galFile.c_str();
    std::ifstream fs(fname);
    if (!fs) {  // An error occurred opening the file.
        std::cerr << "Unable to open file " << fname << std::endl;
        myexit(1);
    }
    // Reserve some storage, since we expect we'll be reading quite a few
    // lines from this file.
    try {P.reserve(4000000);} catch (std::exception& e) {myexception(e);}
    // Skip any leading lines beginning with #
    getline(fs,buf);
    while (fs.eof()==0 && buf[0]=='#') {
        getline(fs,buf);
    }
    int oid=0;
    while (fs.eof()==0) {
        double ra,dec,redshift;  int nobs,priority;
        std::istringstream(buf) >> ra >> dec >> redshift >> priority >> nobs ;
        // use priority as proxy for id
        if (ra<   0.) {ra += 360.;}
        if (ra>=360.) {ra -= 360.;}
        if (dec<=-90. || dec>=90.) {
            std::cout << "DEC="<<dec<<" out of range reading "<<fname<<std::endl;
            myexit(1);
        }
        double theta = (90.0 - dec)*M_PI/180.;
        double phi   = (ra        )*M_PI/180.;
        struct galaxy Q;
        Q.nhat[0]    = cos(phi)*sin(theta);
        Q.nhat[1]    = sin(phi)*sin(theta);
        Q.nhat[2]    = cos(theta);
        Q.id         = priority-1;//priority is proxy for id, starts at zero
        //Q.priority   = priority;
        //Q.nobs       = nobs;
        Q.z = redshift;
        Q.ra = ra;
        Q.dec = dec;


        if (oid%F.moduloGal == 0) {
        try{P.push_back(Q);}catch(std::exception& e) {myexception(e);}

        }

        oid++;
        getline(fs,buf);
    }
    fs.close();
    return(P);
}
//to order galaxies by their priority
bool galaxy_priority(target t1,target t2){return (t1.t_priority<t2.t_priority);}

str galaxy::kind(const Feat& F) const {
	return(F.kind[id]);
}

// targets -----------------------------------------------------------------------
// derived from G, but includes priority and nobs_remain
void make_MTL(const Gals& G, const Feat& F,  MTL& M){
    
    int Nobj=G.size();
    struct target targ;
    int special_count(0);
    int stop_count(0);
    for(int i=0;i<Nobj;++i){
        targ.id=i;
        targ.nhat[0]=G[i].nhat[0];
        targ.nhat[1]=G[i].nhat[1];
        targ.nhat[2]=G[i].nhat[2];
        targ.ra=G[i].ra;
        targ.dec=G[i].dec;
        targ.t_priority=F.prio[G[i].id];
        targ.nobs_remain=F.goal[G[i].id];//needs to be goal prior to knowledge!!
        targ.nobs_done=0;//need to keep track of this, too
        targ.once_obs=0;//changed only in update_plan
        targ.SS=F.SS[G[i].id];
        targ.SF=F.SF[G[i].id];

        targ.lastpass=F.lastpass[G[i].id];
        //make list of priorities
        if(targ.dec<F.MaxDec && targ.dec>F.MinDec &&targ.ra<F.MaxRa && targ.ra>F.MinRa){
            M.push_back(targ);

        }
        int g=M.size()-1;

    }
    
}
void make_MTL_SS_SF(const Gals& G, MTL& Targ, MTL& SStars, MTL& SkyF, Gals& Secret, const Feat& F){
    // Targ contains only galaxy targets
    // SStars contains only standard stars
    // SkyF contains only sky fibers
    int Nobj=G.size();
    struct target targ;
    int special_count(0);
    int stop_count(0);
    for(int i=0;i<Nobj;++i){
        targ.id=i;
        targ.nhat[0]=G[i].nhat[0];
        targ.nhat[1]=G[i].nhat[1];
        targ.nhat[2]=G[i].nhat[2];
        targ.ra=G[i].ra;
        targ.dec=G[i].dec;
        targ.t_priority=F.prio[G[i].id];
        targ.nobs_remain=F.goal[G[i].id];//needs to be goal prior to knowledge!!
        targ.nobs_done=0;//need to keep track of this, too
        targ.once_obs=0;//changed only in update_plan
        targ.SS=F.SS[G[i].id];
        targ.SF=F.SF[G[i].id];
        
        targ.lastpass=F.lastpass[G[i].id];
        //make list of priorities
        if(targ.dec<F.MaxDec && targ.dec>F.MinDec &&targ.ra<F.MaxRa && targ.ra>F.MinRa){
        
            if(targ.SS)SStars.push_back(targ);
            else if(targ.SF)SkyF.push_back(targ);
            else {
                Targ.push_back(targ);
                Secret.push_back(G[i]);
            }
        }
    }
    
}


/*
void write_MTLfile(const Gals& Secret, const MTL& M,const Feat& F){
    FILE * FA;
    str sa=F.MTLfile;
    FA = fopen(sa.c_str(),"w");

    for (int i=0;i<M.size();++i){
        fprintf(FA," %d MartinsMocks %f  %f  %d  %d %d \n",i,M[i].ra,M[i].dec,M[i].nobs_remain,M[i].t_priority,M[i].lastpass);
    }
    fclose(FA);
}
*/
void write_MTL_SS_SFfile(const MTL& Targ, const MTL& SStars,const MTL& SkyF,const Gals& Secret, const Feat& F){
    FILE * FA;
    str sa=F.Targfile;
    FA = fopen(sa.c_str(),"w");
    //str source="MartinsMocks";
    for (int i=0;i<Targ.size();++i){
        fprintf(FA," %d Target %f  %f  %d  %d %d \n",Targ[i].id,Targ[i].ra,Targ[i].dec,Targ[i].nobs_remain,Targ[i].t_priority,Targ[i].lastpass);
    }
    fclose(FA);
    FILE * FB;
    str sb=F.SStarsfile;
    FB = fopen(sb.c_str(),"w");
    //str source="MartinsMocks";
    for (int i=0;i<SStars.size();++i){
        fprintf(FB," %d SStars %f  %f  %d  %d %d \n",SStars[i].id,SStars[i].ra,SStars[i].dec,SStars[i].nobs_remain,SStars[i].t_priority,SStars[i].lastpass);
    }
    fclose(FB);
    FILE * FC;
    str sc=F.SkyFfile;
    FC = fopen(sc.c_str(),"w");
    //str source="MartinsMocks";
    for (int i=0;i<SkyF.size();++i){
        fprintf(FC," %d SkyF %f  %f  %d  %d %d \n",SkyF[i].id,SkyF[i].ra,SkyF[i].dec,SkyF[i].nobs_remain,SkyF[i].t_priority,SkyF[i].lastpass);
    }
    fclose(FC);
    FILE * FD;
    str sd=F.Secretfile;
    FD = fopen(sd.c_str(),"w");
    for (int i=0;i<Secret.size();++i){
        fprintf(FD," %d Secret %f  %f  %d   \n",      i,Secret[i].ra,Secret[i].dec,Secret[i].id);
    }
    fclose(FD);

}

Gals read_Secretfile(str readfile, const Feat&F){
    str s=readfile;
    Gals Secret;
    std::string buf;
    const char* fname;
    fname= s.c_str();
    std::ifstream fs(fname);
    if (!fs) {  // An error occurred opening the file.
        std::cerr << "Unable to open MTLfile " << fname << std::endl;
        myexit(1);
    }
    // Reserve some storage, since we expect we'll be reading quite a few
    // lines from this file.
    try {Secret.reserve(4000000);} catch (std::exception& e) {myexception(e);}
    // Skip any leading lines beginning with #
    getline(fs,buf);
    while (fs.eof()==0 && buf[0]=='#') {
        getline(fs,buf);
    }
    while (fs.eof()==0) {
        double ra,dec;
        int id, i;
        str xname;
        std::istringstream(buf)>> i>>xname>> ra >> dec>> id ;

        if (ra<   0.) {ra += 360.;}
        if (ra>=360.) {ra -= 360.;}
        if (dec<=-90. || dec>=90.) {
            std::cout << "DEC="<<dec<<" out of range reading "<<fname<<std::endl;
            myexit(1);
        }
        struct galaxy Q;
        Q.ra = ra;
        Q.dec = dec;
        Q.id = id;
        Secret.push_back(Q);
        getline(fs,buf);
    }
    return Secret;
}


MTL read_MTLfile(str readfile, const Feat& F, int SS, int SF){
    //str s=F.MTLfile;
    str s=readfile;
    MTL M;
    std::string buf;
    const char* fname;
    fname= s.c_str();
    std::ifstream fs(fname);
    if (!fs) {  // An error occurred opening the file.
        std::cerr << "Unable to open MTLfile " << fname << std::endl;
        myexit(1);
        }
        // Reserve some storage, since we expect we'll be reading quite a few
        // lines from this file.
        try {M.reserve(4000000);} catch (std::exception& e) {myexception(e);}
        // Skip any leading lines beginning with #
        getline(fs,buf);
        while (fs.eof()==0 && buf[0]=='#') {
            getline(fs,buf);
        }
        while (fs.eof()==0) {
            double ra,dec;
            int id, nobs_remain, priority, lastpass;
            str xname;
            std::istringstream(buf)>> id>>xname>> ra >> dec>> nobs_remain >>  priority >>lastpass ;
            //std::istringstream(buf)>> id>> xname>>ra >> dec >>  nobs_remain>> priority;
            if (ra<   0.) {ra += 360.;}
            if (ra>=360.) {ra -= 360.;}
            if (dec<=-90. || dec>=90.) {
                std::cout << "DEC="<<dec<<" out of range reading "<<fname<<std::endl;
                myexit(1);
            }
            double theta = (90.0 - dec)*M_PI/180.;
            double phi   = (ra        )*M_PI/180.;
            struct target Q;
            Q.nhat[0]    = cos(phi)*sin(theta);
            Q.nhat[1]    = sin(phi)*sin(theta);
            Q.nhat[2]    = cos(theta);
            Q.t_priority = priority;//priority is proxy for id, starts at zero
            Q.nobs_remain= nobs_remain;
            Q.nobs_done=0;//need to keep track of this, too
            Q.once_obs=0;//changed only in update_plan
            Q.ra = ra;
            Q.dec = dec;
            Q.id = id;
            Q.lastpass = lastpass;
            Q.SS=SS;
            Q.SF=SF;
            
            if (id%F.moduloGal == 0) {
                try{M.push_back(Q);}catch(std::exception& e) {myexception(e);}
            }
            bool in=false;
            for (int j=0;j<M.priority_list.size();++j){
                if(Q.t_priority==M.priority_list[j]){in=true;
                }
            }
            if(!in){
                M.priority_list.push_back(Q.t_priority);
            }

            getline(fs,buf);
        }
    std::sort(M.priority_list.begin(),M.priority_list.end());
    fs.close();
    return(M);
    }
void assign_priority_class(MTL& M){
    // assign each target to a priority class
    //this needs to be updated
    for(int i=0;i<M.size();++i){
        for(int j=0;j<M.priority_list.size();++j){
            if(M[i].t_priority==M.priority_list[j]){
                M[i].priority_class=j;}
        }
    }
}


// PP ----------------------------------------------------------------------------
// Read the positions of the fibers on each plate.  Assumes the format
// of fiberpos.txt produced by "F.Randomize_fibers".
// need also to get the petal, i.e. spectrometer  rnc 1/16/15  added S
void PP::read_fiber_positions(const Feat& F) {
    std::string buf;
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

void PP::compute_fibsofsp(const Feat& F) {//list of fibers on each spectrometer
	fibers_of_sp.resize(F.Npetal);
	for (int k=0; k<F.Nfiber; k++) fibers_of_sp[spectrom[k]].push_back(k);
}

List PP::fibs_of_same_pet(int k) const {//all fibers on same spectrometer as fiber k
	return(fibers_of_sp[spectrom[k]]);
}

dpair PP::coords(int k) const {// x and y in mm
	return dpair(fp[2*k],fp[2*k+1]);
}

// plate ---------------------------------------------------------------------------

List plate::av_gals_plate(const Feat& F,const MTL& M, const PP& pp) const {//list of galaxies available to plate no repetitions
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
Plates read_plate_centers(const Feat& F) {
	Plates P;
        std::string buf;
	std::ifstream fs(F.tileFile.c_str());
	if (!fs) {  // An error occurred opening the file.
		std::cerr << "Unable to open file " << F.tileFile << std::endl;
		myexit(1);
	}
	// Reserve some storage, since we expect we'll be reading quite a few
	// lines from this file.
	try {P.reserve(4000000);} catch (std::exception& e) {myexception(e);}

	double ra,dec,ebv,airmass,exposefac;
    int ipass,in_desi,tileid;
	int l = 0;
	while (fs.eof()==0) {
		getline(fs,buf);
                if(buf.compare(0, 7, "STRUCT1") != 0) {
                    // std::cout << "Skipping " << buf << std::endl;
                    continue;
                }
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
			Q.tileid = tileid;
                        //std::cout << "TILEID " << tileid << std::endl;
			l++;
			Q.tilera        = ra;
			Q.tiledec       = dec;
			Q.nhat[0]    = sin(theta)*cos(phi);
			Q.nhat[1]    = sin(theta)*sin(phi);
			Q.nhat[2]    = cos(theta);
			Q.ipass      = ipass-1; // <- be careful, format of input file
			Q.av_gals.resize(F.Nfiber); // <- added
			Q.density.resize(F.Nfiber); // <- added
            Q.SS_av_gal.resize(F.Nfbp);
            Q.SF_av_gal.resize(F.Nfbp);
            Q.SS_in_petal.resize(F.Npetal);
            Q.SF_in_petal.resize(F.Npetal);
            for (int i=0;i<F.Npetal;++i){Q.SS_in_petal[i]=0;}
            for (int i=0;i<F.Npetal;++i){Q.SF_in_petal[i]=0;}
            //if(dec<F.MaxDec && dec>F.MinDec &&ra<F.MaxRa && ra>F.MinRa){
                try {P.push_back(Q);} catch(std::exception& e) {myexception(e);
                //}
            }
		}
	}
	fs.close();
	return(P);
}
// Assignment -----------------------------------------------------------------------------
Assignment::Assignment(const MTL& M, const Feat& F) {

    
    TF=initTable(F.Nplate,F.Nfiber,-1);//galaxy assigned to tile-fiber TF[j][k]
	GL = initPtable(F.Ngal,0); //tile-fiber pair for galaxy  GL[g]
    printf("test 1\n");
    order.resize(F.Nplate);
    for (int i=0; i<F.Nplate; i++) order[i] = i;
    printf("test 2\n");
	next_plate = 0;
	kinds = initCube(F.Nplate,F.Npetal,F.Categories);
    printf("test 3\n");
	unused = initTable(F.Nplate,F.Npetal,F.Nfbp);//initialized to number of fibers on a petal
    printf("test 4\n");
    }

Assignment::~Assignment() {}

// Assign g with tile/fiber (j,k), and check for duplicates
void Assignment::assign(int j, int k, int g, MTL& M, Plates& P, const PP& pp) {
	// Assign (j,k)
	int q = TF[j][k];
	if (q != -1) {
		printf("### !!! ### DUPLICATE (j,k) = (%d,%d) assigned with g = %d and %d ---> information on first g lost \n",j,k,q,g);
		myexit(1);
	}
    //DIAGNOSTIC

	// Assign g
    TF[j][k]=g;
	Plist pl = GL[g];//pair list, tf's for this g
	pair p = pair(j,k);
	for(int i=0;i<pl.size();i++){
		if(pl[i].f==j){

		printf("### !!! ### DUPLICATE g = %d assigned with (j,k) = (%d,%d) and (%d,%d) ---> information on first (j,k) lost \n",g,pl[i].f,pl[i].s,j,k);
		
		myexit(1); // Can be commented if want to force continuing
		}
	}
	GL[g].push_back(p);
    M[g].nobs_done++;
    M[g].nobs_remain--;
    if(M[g].SF){
        int q=pp.spectrom[k];
        P[j].SF_in_petal[q]+=1;}
    if(M[g].SS){
        int q=pp.spectrom[k];
        P[j].SS_in_petal[q]+=1;}
	unused[j][pp.spectrom[k]]--;
}

void Assignment::unassign(int j, int k, int g, MTL& M, Plates& P, const PP& pp) {
	if (TF[j][k]==-1) printf("### !!! ### TF (j,k) = (%d,%d) gets unassigned but was already not assigned\n",j,k);
	int a = isfound(pair(j,k),GL[g]);
	if (a==-1) printf("### !!! ### Galaxy g = %d gets unassigned but was already not assigned\n",g);

	TF[j][k] = -1;
	if (a!=-1) erase(a,GL[g]);
    M[g].nobs_done--;
    M[g].nobs_remain++;
    if(M[g].SF){
        int p=pp.spectrom[k];
        P[j].SF_in_petal[p]-=1;}
    if(M[g].SS){
        int p=pp.spectrom[k];
        P[j].SS_in_petal[p]-=1;}

	unused[j][pp.spectrom[k]]++;
}

void Assignment::verif(const Plates& P, const MTL& M, const PP& pp, const Feat& F) const {
	str qso_lrgA[] = {"QSOLy-a","QSOTracer","FakeQSO","LRG","FakeLRG"}; List qso_lrg = F.init_ids_list(qso_lrgA,5);
	for (int g=0; g<F.Ngal; g++) {// make sure observations are separated by at least InterPlate
		Plist tfs = GL[g];
		int j0(-1); int j1(-1);
		for (int i=0; i<tfs.size(); i++) {
			pair tf = tfs[i];
			int j0 = j1;
			int j1 = tf.f;
			// Verif on TF
			if (TF[tf.f][tf.s]!=g) { printf("ERROR in verification of correspondance of galaxies !\n"); fl(); }
			// No 2 assignments within an interval of F.InterPlate
			//if (j0!=-1 && isfound(M[g].id,qso_lrg) && fabs(j1-j0)<F.InterPlate) { printf("ERROR in verification of F.InterPlate g=%d with j=%d and %d\n",g,j0,j1); fl(); }
		}
	}
	for (int j=0; j<F.Nplate; j++) {
		List gals = initList(F.Ngal);
		for (int k=0; k<F.Nfiber; k++) {
			int g = TF[j][k];
			if (g!=-1) {
				// Verif on GL: is it consistent with TF?
				if (isfound(pair(j,k),GL[g])==-1) { printf("ERROR in verification of correspondance of tfs !\n"); fl(); }
				// Verif that a galaxy isn't observed twice
				if (gals[g]==1) printf("ERROR in verification, twice the same galaxy by (%d,%d)\n",j,k);
				else gals[g] = 1;
				// Collision checking
				if (!F.Collision && is_collision(j,k,pp,M,P,F)!=-1) printf("ERROR in verification : collisions\n");
			}
		}
	}
	for (int j=0; j<F.Nplate; j++) {// right number of SS and SF
		for (int n=0; n<F.Npetal; n++){
			if (kinds[j][n][F.ids.at("SS")]!=F.MaxSS || kinds[j][n][F.ids.at("SF")]!=F.MaxSF) printf("ERROR in verification : number of SF or SS\n");
		}
	}
}

int Assignment::is_assigned_jg(int j, int g) const {// is galaxy g assigned on tile j
	for (int i=0; i<GL[g].size(); i++) if (GL[g][i].f == j) return i;
	return -1;
}

int Assignment::is_assigned_jg(int j, int g, const MTL& M, const Feat& F) const { // No occurrence too nearby in tiles
	for (int i=0; i<GL[g].size(); i++) if (fabs(j-GL[g][i].f)<F.InterPlate || j==GL[g][i].f) return i; 
	return -1;
}

bool Assignment::is_assigned_tf(int j, int k) const { return (TF[j][k] != -1); }

int Assignment::na(const Feat& F, int begin, int size) const {//unassigned fibers in tiles begin to begin+size
	int size1 = (size==-1) ? F.Nplate : size;
	int cnt(0);
	for (int j=begin; j<begin+size1; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (TF[j][k]!=-1) cnt++;
		}
	}
	return cnt;
}

Plist Assignment::chosen_tfs(int g, const Feat& F, int begin, int size0) const {//creates list of tile-fibers observing g
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

Table Assignment::unused_fbp(const PP& pp, const Feat& F) const {//table unused fibers on petal  on tile  j
	Table unused = initTable(F.Nplate,F.Npetal);
	List Sp = pp.spectrom;
	for(int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j][Sp[k]]++;
		}
	}
	return unused;
}
//not used
List Assignment::unused_f(const Feat& F) const {//total unused fibers
	List unused = initList(F.Nplate);
	for(int j=0; j<F.Nplate; j++) {
		for (int k=0; k<F.Nfiber; k++) {
			if (!is_assigned_tf(j,k)) unused[j]++;
		}
	}
	return unused;
}

int Assignment::unused_f(int j, const Feat& F) const {//unused fibers on tile j
	int unused(0);
	for (int k=0; k<F.Nfiber; k++) if (!is_assigned_tf(j,k)) unused++;
	return unused;
}
int Assignment::unused_fbp(int j, int k, const PP& pp, const Feat& F) const {//unused fibers on petal contain fiber k, tile j
	List fibs = pp.fibers_of_sp[pp.spectrom[k]];
	int unused(0);
	for (int i=0; i<fibs.size(); i++) {
		if (!is_assigned_tf(j,fibs[i])) unused++;
	}
	return unused;
}

int Assignment::nkind(int j, int k, int kind, const MTL& M, const Plates& P, const PP& pp, const Feat& F, bool pet) const {
	//if pet is false, used petal of fiber k,, if pet is true use petal k        
	if (!pet) return kinds[j][pp.spectrom[k]][kind];
	else return kinds[j][k][kind];
}
List Assignment::fibs_unassigned(int j, int pet, const MTL& M, const PP& pp, const Feat& F) const {//list of unassigned fibers on petal pet
	List L;
	List fibs = pp.fibers_of_sp[pet];
	for (int kk=0; kk<F.Nfbp; kk++) {
		int k = fibs[kk];
		if (!is_assigned_tf(j,k)) L.push_back(k);
	}
	return L;
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
struct onplate change_coords(const struct target& O, const struct plate& P) {
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

// (On plate p) finds if there is a collision if fiber k would observe galaxy g (collision with neighbor)
//  j is in list that runs to F.ONplate since it is used in TF[j][k]
int Assignment::find_collision(int j, int k, int g, const PP& pp, const MTL& M, const Plates& P, const Feat& F, int col) const {//check all neighboring fibers
	bool bol = (col==-1) ? F.Collision : false;
	if (bol) return -1;
	dpair G1 = projection(g,j,M,P);
	for (int i=0; i<pp.N[k].size(); i++) {// i numbers the fibers neighboring fiber k
		int kn = pp.N[k][i];
		int gn = TF[j][kn];
		if (gn!=-1) {
			dpair G2 = projection(gn,j,M,P);
			bool b = F.Exact ? collision(pp.coords(k),G1,pp.coords(kn),G2,F) : (sq(G1,G2) < sq(F.AvCollide));
			if (b) return kn;
		}
	}
	return -1;
}
//probably not used
bool Assignment::find_collision(int j, int k, int kn, int g, int gn, const PP& pp, const MTL& M, const Plates& P, const Feat& F, int col) const {//check two fibers
	bool bol = (col==-1) ? F.Collision : false;
	if (bol) return false;
	dpair G1 = projection(g,j,M,P);
	dpair G2 = projection(gn,j,M,P);
	return F.Exact ? collision(pp.coords(k),G1,pp.coords(kn),G2,F) : (sq(G1,G2) < sq(F.AvCollide));
}

int Assignment::is_collision(int j, int k, const PP& pp, const MTL& M, const Plates& P, const Feat& F) const {//find collision for galaxy g
	int g = TF[j][k];
	if (g!=-1) return find_collision(j,k,g,pp,M,P,F,0);
	else return -1;
}

float Assignment::colrate(const PP& pp, const MTL& M, const Plates& P, const Feat& F, int jend0) const {//rate of collisions
	int jend = (jend0==-1) ? F.Nplate : jend0;
	int col = 0;
	for (int j=0; j<jend; j++) {
		List done = initList(F.Nfiber);
		for (int k=0; k<F.Nfiber; k++) {
			if (done[k] == 0) {
				int c = is_collision(j,k,pp,M,P,F);
				if (c!=-1) {
					done[c] = 1;
					col += 2;
				}
			}
		}
	}
	return percent(col,jend*F.Nfiber);
}

dpair projection(int g, int j, const MTL& M, const Plates& OP) {//x and y coordinates for galaxy observed on plate j
    // USE OLD LIST OF PLATES HERE
	struct onplate op = change_coords(M[g],OP[j]);
    if (op.pos[0]*op.pos[0]+op.pos[1]*op.pos[1]>500.*500.){
        printf("outside positioner range  g  %d  j  %d  x %f  y %f\n",g,j,op.pos[0],op.pos[1]);
    }
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
