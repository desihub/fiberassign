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
#include        <unistd.h> // getpid()
#include	"../macros.h"
#include	"../misc.h"
#include	"collide.h"

double pi2;

int main(int argc, char **argv) {
	srand(time(NULL)+getpid());
	pi2 = 6.283185;
	int times = 10000000;
	PosP posp(1.339,0.990,3.0,3.0,10.392);

	Dplist cases;
	int cnt = 0;
	for (int i=0; i<times; i++) {
		Schema sch = hexag(posp);
		int col = colliding(0,sch,posp);
		if (col!=-1) {
			cnt++;
			cases.push_back(dpair(sch.dist_fib(0,col,posp),cartesian_weight(sch.poss[0],sch.poss[col])));
		} 
	}
	printf("%f \n",percent(cnt,times));
	double interval = 0.01;
	
	Dlist hist = histogram(cases,interval);
	Dtable Dt; Dt.push_back(percents(hist,sumlist(hist)));
	print_mult_Dtable_latex("collisions","doc/figs/col.dat",Dt,interval);

	return(0);
}
