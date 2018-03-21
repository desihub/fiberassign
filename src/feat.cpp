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
#include <string.h>
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
    MaxSS = 10;
    MaxSF = 40;
    InterPlate = 0;
    Collision = false;
    Exact = true;
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
            if (tok[0]=="Targfile") Targfile=tok[1];
            if (tok[0]=="tileFile") tileFile= tok[1];
            if (tok[0]=="fibFile") fibFile= tok[1];
            if (tok[0]=="fibstatusFile") fibstatusFile = tok[1];
            if (tok[0]=="surveyFile") surveyFile= tok[1];
            if (tok[0]=="outDir") outDir= tok[1];
            if (tok[0]=="SStarsfile")SStarsfile=tok[1];
            if (tok[0]=="SkyFfile") SkyFfile=tok[1];
            if (tok[0]=="runDate") runDate=tok[1];
                
    }
  }  
  fIn.close();
}

void Feat::parseCommandLine(int argc, char **argv) {
    int i;
    for (i=1;i<argc;){
        if (!strcmp(argv[i],"-target")){
        i++;
         Targfile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-sky")){
        i++;
        SkyFfile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-star")){
        i++;
            SStarsfile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-survey")){
        i++;
         surveyFile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-outdir")){
        i++;
        outDir = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-tilefile")){
        i++;
        tileFile= str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-fibfile")){
        i++;
        fibFile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-fibstatusfile")){
        i++;
        fibstatusFile = str(argv[i]);
        i++;    
        }else if (!strcmp(argv[i],"-rundate")){
        i++;
        fibstatusFile = str(argv[i]);
        i++;   
        }   
    }
}
