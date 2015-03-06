#ifndef MACROS_H
#define MACROS_H

#define		PlateRadius	1.65	// In degrees
#define		PatrolMaxRad	 6.0	// In mm
#define		PatrolMinRad	 0.0	// In mm
//rnc 08/04/14
#define		MaxFiber	5000	// Maxm available fibers per plate
#define     MaxPetal    10      //10 petals in focal plane

#define     MaxPass      5
#define     MaxObs       5  //goal for Ly-a
#define     Categories   8 
#define     collide  2.1 // stay-clear for fiber collisions in mm
#define     ntimes 10 //number of intervals (on tiles) for time plot
#define     ASCIICENTERS 1//read tile locations in sky in ASCII format

#define     Npass 5 //number of pass
#define     NeighbourRad 11. //number of pass
#define     MinUnused 50 // Number of minimum unused fibers by petals
#define     TotalArea 15789. 
#define     MaxSS 10
#define     MaxSF 40

extern int MaxPlate;
extern int Ngal;
extern int MaxPrio;
#endif
