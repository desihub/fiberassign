#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <sys/time.h>
#include <map>
#include <stdlib.h>     /* srand, rand */
#include <string.h>
#include "misc.h"
#include "collision.h"
#include "feat.h"


void Usage (char * ExecName) {
    std::cout << "Usage:  " << ExecName;
    std::cout << " --mtl <target_filename> ";
    std::cout << " --sky <sky_filename> ";
    std::cout << " --stdstar <star_filename> ";
    std::cout << " --surveytiles <surveytiles_filename> ";
    std::cout << " --footprint <tilelist_filename> ";
    std::cout << " --positioners <fiberpos_filename> ";
    std::cout << " --fibstatusfile <fiberstatus_filename> ";
    std::cout << " --outdir <outputdirectory> ";
    std::cout << " --starmask <starbit_mask]";
    std::cout << " [--rundate <YYYY-MM-DD>]" << std::endl;
    exit(0);
}

// Features ------------------------------------------------------------------
Feat::Feat () {
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
    StarMask = 60129542144;
}

void Feat::readInputFile (const char file[]) {
    const int Mc = 512;  // Max chars per line
    char delimiter = ' ';
    std::ifstream fIn;
    fIn.open(file);  // open a file
    if (!fIn.good() ) myexit(1); // Not found
    while (!fIn.eof() ) {
        char buf[Mc];
        fIn.getline(buf, Mc);
        Slist tok = s2vec(buf, delimiter);
        if (2 <= tok.size() ) {
            if (tok[0] == "Targfile") Targfile = tok[1];
            if (tok[0] == "tileFile") tileFile = tok[1];
            if (tok[0] == "fibFile") fibFile = tok[1];
            if (tok[0] == "fibstatusFile") fibstatusFile = tok[1];
            if (tok[0] == "surveyFile") surveyFile = tok[1];
            if (tok[0] == "outDir") outDir = tok[1];
            if (tok[0] == "SStarsfile") SStarsfile = tok[1];
            if (tok[0] == "SkyFfile") SkyFfile = tok[1];
            if (tok[0] == "runDate") runDate = tok[1];
        }
    }
    fIn.close();
}

void Feat::parseCommandLine (int argc, char * * argv) {
    int i;
    for (i = 1; i < argc; ) {
        if (!strcmp(argv[i], "--mtl") ) {
            i++;
            Targfile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--sky") ) {
            i++;
            SkyFfile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--stdstar") ) {
            i++;
            SStarsfile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--surveytiles") ) {
            i++;
            surveyFile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--outdir") ) {
            i++;
            outDir = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--footprint") ) {
            i++;
            tileFile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--positioners") ) {
            i++;
            fibFile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--fibstatusfile") ) {
            i++;
            fibstatusFile = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--rundate") ) {
            i++;
            runDate = str(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--starmask") ) {
            i++;
            StarMask = atol(argv[i]);
            i++;
        } else if (!strcmp(argv[i], "--help") ) {
            Usage(argv[0]);
        } else if (!strcmp(argv[i], "-h") ) {
            Usage(argv[0]);
        } else {
            fprintf (stderr, "\nUnrecognized option: %s\n\n", argv[i]);
            Usage(argv[0]);
        }
    }
}
