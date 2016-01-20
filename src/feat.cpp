#include	<cstdlib>
#include	<cmath>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include        <cstring>
#include	<iomanip>
#include	<string>
#include	<vector>
#include	<algorithm>
#include	<exception>
#include	<sys/time.h>
#include        <map>
#include        <stdlib.h>     /* srand, rand */
#include        "misc.h"
#include        "collision.h"
#include        "feat.h"
// Features ------------------------------------------------------------------
Feat::Feat() { 
	Count = 0; 

}

int Feat::id(str s) const {
	for (int i=0; i<Categories; i++) if (kind[i]==s) return i;
	std::cout << "ERROR in Feat id(), string not found in kind" << std::endl;
	return -1;
}

void Feat::init_ids() {
	for (int i=0; i<Categories; i++) ids[kind[i]] = i;
}

void Feat::init_ids_types() {
	for (int i=0; i<Categories; i++) ids_types[kind[i]] = i;
}

List Feat::init_ids_list(str l[], int size) const {
	List L;
	L.resize(size);
	for (int i=0; i<size; i++) L[i] = ids.at(l[i]);
	return L;
}

void Feat::init_types() {
	for (int i=0; i<type.size(); i++) if (!isfound(type[i],types)) types.push_back(type[i]);
}

void Feat::init_no_ss_sf() {
	//str noA[] = {"LRG","FakeLRG","ELG","QSOLy-a","QSOTracer","FakeQSO"};
    for(int i=0;i<kind.size()-2;++i)no_ss_sf.push_back(i);
	//no_ss_sf = init_ids_list(noA,6);
}

void Feat::init_ss_sf() {
	str A[] = {"SS","SF"}; 
	ss_sf = init_ids_list(A,2);
}
		
bool Feat::iftype(int kind, str type) const {
	return ids_types.at(type)==type[kind];
}

void Feat::readInputFile(const char file[]) {
	const int Mc = 512; // Max chars per line
	char delimiter = ' ';

	std::ifstream fIn;
	fIn.open(file); // open a file
	if (!fIn.good()) myexit(1); // Not found
	while (!fIn.eof()) {
		char buf[Mc];
		fIn.getline(buf,Mc);
		int n = 0; // a for-loop index
		Slist tok = s2vec(buf,delimiter);
		if (2<=tok.size()) {
			if (tok[0]=="galFile") galFile = tok[1];
			if (tok[0]=="tileFile") tileFile= tok[1];
			if (tok[0]=="fibFile") fibFile= tok[1];
			if (tok[0]=="outDir") outDir= tok[1];
			if (tok[0]=="PrintAscii") PrintAscii= s2b(tok[1]);
            if (tok[0]=="PrintFits") PrintFits= s2b(tok[1]);
            
            if (tok[0]=="MTLfile") MTLfile=tok[1];
            if (tok[0]=="Targfile") Targfile=tok[1];
            if (tok[0]=="SStarsfile")SStarsfile=tok[1];
            if (tok[0]=="SkyFfile") SkyFfile=tok[1];
            if (tok[0]=="Secretfile") Secretfile=tok[1];
            
            if (tok[0]=="diagnose") diagnose=s2b(tok[1]);

			if (tok[0]=="kind") {
				Categories = tok.size()-1;
                for (int i=0; i<Categories; i++) kind.push_back(tok[i+1]);
				init_ids();
				init_ss_sf();
				init_no_ss_sf();
			}
			if (tok[0]=="type") {
				Categories = tok.size()-1;
				for (int i=0; i<Categories; i++) type.push_back(tok[i+1]);
				init_types();
				init_ids_types();
			}
            if (tok[0]=="prio")  for (int i=0; i<Categories; i++){prio.push_back(s2i(tok[i+1]));}
			if (tok[0]=="priopost") for (int i=0; i<Categories; i++) priopost.push_back(s2i(tok[i+1]));
			if (tok[0]=="goal") for (int i=0; i<Categories; i++) goal.push_back(s2i(tok[i+1]));
            if (tok[0]=="goalpost") for (int i=0; i<Categories; i++) goalpost.push_back(s2i(tok[i+1]));
            if (tok[0]=="lastpass") for(int i=0; i<Categories;i++)lastpass.push_back(s2i(tok[i+1]));
            if (tok[0]=="SS") for(int i=0; i<Categories;i++)SS.push_back(s2i(tok[i+1]));
            if (tok[0]=="SF") for(int i=0; i<Categories;i++)SF.push_back(s2i(tok[i+1]));
            if (tok[0]=="pass_intervals"){
                int n_intervals=tok.size()-1;
                for (int i=0;i<n_intervals;++i) pass_intervals.push_back(s2i(tok[i+1]));
            }
            
      if (tok[0]=="InterPlate") InterPlate = s2i(tok[1]);
      if (tok[0]=="Randomize") Randomize = s2b(tok[1]);
      if (tok[0]=="Pacman") Pacman = s2b(tok[1]);
      if (tok[0]=="Npass") Npass = s2i(tok[1]);
      if (tok[0]=="MaxSS") MaxSS = s2i(tok[1]);
      if (tok[0]=="MaxSF") MaxSF = s2i(tok[1]);
      if (tok[0]=="PlateRadius") PlateRadius = s2d(tok[1]);
      if (tok[0]=="Analysis") Analysis = s2i(tok[1]);
      if (tok[0]=="InfDens") InfDens = s2b(tok[1]);
      
      if (tok[0]=="TotalArea") TotalArea = s2d(tok[1]);
      if (tok[0]=="invFibArea") invFibArea = s2d(tok[1]);
      if (tok[0]=="moduloGal") moduloGal = s2i(tok[1]);
      if (tok[0]=="moduloFiber") moduloFiber = s2i(tok[1]);
      
      if (tok[0]=="Collision") Collision = s2b(tok[1]);
      if (tok[0]=="Exact") Exact = s2b(tok[1]);
      if (tok[0]=="AvCollide") AvCollide = s2d(tok[1]);
      if (tok[0]=="Collide") Collide = s2d(tok[1]);
      if (tok[0]=="NoCollide") NoCollide = s2d(tok[1]);
      if (tok[0]=="PatrolRad") PatrolRad = s2d(tok[1]);
      if (tok[0]=="NeighborRad") NeighborRad = s2d(tok[1]);
      
      if (tok[0]=="PlotObsTime") PlotObsTime = s2b(tok[1]);
      if (tok[0]=="PlotHistLya") PlotHistLya = s2b(tok[1]);
      if (tok[0]=="PlotDistLya") PlotDistLya = s2b(tok[1]);
      if (tok[0]=="PlotFreeFibHist") PlotFreeFibHist = s2b(tok[1]);
      if (tok[0]=="PlotFreeFibTime") PlotFreeFibTime = s2b(tok[1]);
      if (tok[0]=="PlotSeenDens") PlotSeenDens = s2b(tok[1]);
      if (tok[0]=="Verif") Verif = s2b(tok[1]);
      if (tok[0]=="Ascii") Ascii = s2b(tok[1]);
      if (tok[0]=="PrintGalObs") PrintGalObs = s2i(tok[1]);
      if (tok[0]=="BrightTime") BrightTime = s2b(tok[1]);
      
      if (tok[0]=="MaxDec") MaxDec = s2d(tok[1]);
      if (tok[0]=="MinDec") MinDec = s2d(tok[1]);
      if (tok[0]=="MaxRa") MaxRa = s2d(tok[1]);
      if (tok[0]=="MinRa") MinRa = s2d(tok[1]);
      
    }
  }
  
    
  fIn.close();
}
