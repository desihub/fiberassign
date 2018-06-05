#ifndef FEAT_H
#define FEAT_H

#include    <cstdlib>
#include    <cmath>
#include    <fstream>
#include    <sstream>
#include    <iostream>
#include    <iomanip>
#include    <string>
#include    <vector>
#include    <algorithm>
#include    <exception>
#include    <sys/time.h>
#include        <map>
#include        "misc.h"
#include        "collision.h"

// Feat --------------------------------------------------
class Feat { // F for features
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


    long StarMask; //bit in desi_target that signals a standard star
    
    int InterPlate; 
    int MaxSS;
    int MaxSF;
    double PlateRadius;
    int PrintGalObs;
    
    bool Collision; // True when we take collisions into account
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
	int Nfbp; // Number of fibers by petals

    // Memorizes geometry of central body and fiber holder
    polygon cb;
    polygon fh;

    // Permit to count some things
    int Count;

    // Practical used lists
    List no_ss_sf;
    List ss_sf;

    // Methods
    Feat();
    void readInputFile(const char fname[]);
    void parseCommandLine(int argc, char ** argv);
    int id(str s) const;

    void init_ids(); // Init ids member
    void init_ids_types(); // Same
    List init_ids_list(str l[], int size) const;
    void init_types();
    void init_no_ss_sf();
    void init_ss_sf();
    bool iftype(int kind, str type) const; // Whether kind is of type "type"
};
#endif
