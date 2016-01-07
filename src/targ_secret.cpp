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
#include	"modules/htmTree.h"
#include	"modules/kdTree.h"
#include        "misc.h"
#include        "feat.h"
#include        "structs.h"
#include        "collision.h"
#include        "global.h"
//reduce redistributes, updates  07/02/15 rnc
int main(int argc, char **argv) {
	//// Initializations ---------------------------------------------
	srand48(1234); // Make sure we have reproducability
	check_args(argc);
	Time t, time; // t for global, time for local
	init_time(t);
	Feat F;

	// Read parameters file //
	F.readInputFile(argv[1]);
	printFile(argv[1]);
    std::cout.flush();
	// Read galaxies
	Gals G, Secret;
    MTL Targ, SStars, SkyF;
    if(F.Ascii){
        printf("to read galaxies\n");
        std::cout.flush();
        G=read_galaxies_ascii(F);}
    else{
        G = read_galaxies(F);
    }
	F.Ngal = G.size();
    
	printf("# Read %s galaxies from %s \n",f(F.Ngal).c_str(),F.galFile.c_str());
        std::cout.flush();
    std::vector<int> count;
    count=count_galaxies(G);
    printf(" Number of galaxies by type, QSO-Ly-a, QSO-tracers, LRG, ELG, fake QSO, fake LRG, SS, SF\n");
    std::cout.flush();
    for(int i=0;i<count.size();i++){
        if(count[i]>0)printf (" type %d number  %d  \n",i, count[i]);
        std::cout.flush();
    }
    // make MTL
    //read_MTL(G,F,Secret,Targ);
    //write_MTLfile(Secret,Targ,F);
    make_Targ_Secret(G,Targ,Secret,F);
    write_Targ_Secret(Targ,Secret,F);
	return(0);
  
}
