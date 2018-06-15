#ifndef FEAT_H
#define FEAT_H

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <exception>
#include <sys/time.h>
#include <map>
#include "misc.h"
#include "collision.h"

// Feat --------------------------------------------------
// F for features
class Feat {
    public:
        // Set by features file
        str tileFile;
        str fibFile;
        str outDir;
        str Targfile;
        str SStarsfile;
        str SkyFfile;
        str surveyFile;
        str fibstatusFile;
        str runDate;
        // bit in desi_target that signals a standard star
        long StarMask;
        int InterPlate;
        int MaxSS;
        int MaxSF;
        double PlateRadius;
        int PrintGalObs;
        // True when we take collisions into account
        bool Collision;
        bool Exact;
        double AvCollide;
        double Collide;
        double NoCollide;
        double PatrolRad;
        double NeighborRad;
        bool Verif;
        bool Ascii;
        bool BrightTime;
        int Categories;
        // Set after reading other input files
        int Nplate;
        int NUsedplate;
        int Ngal;
        int NSStars;
        int NSkyF;
        int Ntarg;
        int Nfiber;
        int Npetal;
        // Number of fibers by petals
        int Nfbp;
        // Memorizes geometry of central body and fiber holder
        polygon cb;
        polygon fh;
        // Permit to count some things
        int Count;
        // Practical used lists
        List no_ss_sf;
        List ss_sf;
        // Methods
        Feat ();
        void readInputFile (const char fname[]);
        void parseCommandLine (int argc, char * * argv);
        int id (str s) const;
        // Init ids member
        void init_ids ();
        // Same
        void init_ids_types ();
        List init_ids_list (str l[], int size) const;
        void init_types ();
        void init_no_ss_sf ();
        void init_ss_sf ();
        // Whether kind is of type "type"
        bool iftype (int kind, str type) const;
};

#endif
