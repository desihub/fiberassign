#include    <cstdlib>
#include    <cmath>
#include    <fstream>
#include    <sstream>
#include    <iostream>
#include        <cstring>
#include    <iomanip>
#include    <string>
#include    <vector>
#include    <algorithm>
#include    <exception>
#include    <sys/time.h>
#include        <map>
#include        <stdlib.h>     /* srand, rand */
#include        "misc.h"
#include        "collision.h"
#include        "feat.h"
// Features ------------------------------------------------------------------
Feat::Feat() { 
    Count = 0;
    Categories = 0;
    AvCollide = 3.2;
    Collide = 1.98;
    NoCollide = 7.0;
    PatrolRad = 5.8;
    NeighborRad = 14.05;
    PlateRadius = 1.65;
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
            if (tok[0]=="fibstatusFile") fibstatusFile = tok[1];
            if (tok[0]=="surveyFile") surveyFile= tok[1];
			if (tok[0]=="outDir") outDir= tok[1];
            
            if (tok[0]=="MTLfile") MTLfile=tok[1];
            if (tok[0]=="Targfile") Targfile=tok[1];
            if (tok[0]=="SStarsfile")SStarsfile=tok[1];
            if (tok[0]=="SkyFfile") SkyFfile=tok[1];
            if (tok[0]=="Secretfile") Secretfile=tok[1];
            if (tok[0]=="runDate") runDate=tok[1];

            
      if (tok[0]=="InterPlate") InterPlate = s2i(tok[1]);
      if (tok[0]=="Randomize") Randomize = s2b(tok[1]);
      if (tok[0]=="MaxSS") MaxSS = s2i(tok[1]);
      if (tok[0]=="MaxSF") MaxSF = s2i(tok[1]);
      if (tok[0]=="Analysis") Analysis = s2i(tok[1]);

      
      if (tok[0]=="TotalArea") TotalArea = s2d(tok[1]);
      if (tok[0]=="invFibArea") invFibArea = s2d(tok[1]);
      if (tok[0]=="moduloGal") moduloGal = s2i(tok[1]);
      if (tok[0]=="moduloFiber") moduloFiber = s2i(tok[1]);
      
      if (tok[0]=="Collision") Collision = s2b(tok[1]);
      if (tok[0]=="Exact") Exact = s2b(tok[1]);
      
      if (tok[0]=="Verif") Verif = s2b(tok[1]);
      if (tok[0]=="Ascii") Ascii = s2b(tok[1]);
      if (tok[0]=="PrintGalObs") PrintGalObs = s2i(tok[1]);
      if (tok[0]=="BrightTime") BrightTime = s2b(tok[1]);
      
      
    }
  }
  
    
  fIn.close();
}
