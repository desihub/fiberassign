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
#include        "feat.h"
// Features ------------------------------------------------------------------
Feat::Feat() {
	prio.resize(Categories);
	priopost.resize(Categories);
	goal.resize(Categories);
	kind.resize(Categories);
	Count = 0;
}

int Feat::id(str s) const {
	for (int i=0; i<Categories; i++) if (kind[i]==s) return i;
	std::cout << "ERROR in Feat id(), string not found in kind" << std::endl;
	return -1;
}

int Feat::maxgoal(int kind) const {
	int max(goal[kind]); int prio0= prio[kind];
	for (int i=0; i<Categories; i++) if (prio[i]==prio0 && goal[i]>max) max = goal[i];
	return max;
}

List Feat::maxgoal() const {
	List max = initList(Categories,-1);
	for (int i=0; i<Categories; i++) max[i] = maxgoal(i);
	return max;
}

void Feat::init_ids() {
	for (int i=0; i<Categories; i++) ids[kind[i]] = i;
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
		if (tok[0]=="galFile") galFile = tok[1];
		if (tok[0]=="tileFile") tileFile= tok[1];
		if (tok[0]=="fibFile") fibFile= tok[1];

		if (tok[0]=="kind") {
			Categories = tok.size()-1;
			for (int i=0; i<Categories; i++) kind.push_back(tok[i+1]);
			init_ids();
		}
		if (tok[0]=="prio") for (int i=0; i<Categories; i++) prio.push_back(s2i(tok[i+1]));
		if (tok[0]=="priopost") for (int i=0; i<Categories; i++) priopost.push_back(s2i(tok[i+1]));
		if (tok[0]=="goal") for (int i=0; i<Categories; i++) goal.push_back(s2i(tok[i+1]));

		if (tok[0]=="InterPlate") InterPlate = s2i(tok[1]);
		if (tok[0]=="Randomize") Randomize = s2b(tok[1]);
		if (tok[0]=="Pacman") Pacman = s2b(tok[1]);
		if (tok[0]=="Npass") Npass = s2i(tok[1]);
		if (tok[0]=="MaxSS") MaxSS = s2i(tok[1]);
		if (tok[0]=="MaxSF") MaxSF = s2i(tok[1]);
		if (tok[0]=="PlateRadius") PlateRadius = s2d(tok[1]);

		if (tok[0]=="NeighborRad") NeighborRad = s2d(tok[1]);
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
	}
	fIn.close();
}
