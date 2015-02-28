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
#include        "macros.h"
#include        "misc.h"
// pair ------------------------------------------------------
pair::pair() {f = -1; s = -1;}
pair::pair(int a, int b) {f = a; s = b;}
void pair::setnull() {f = -1; s = -1;}
bool pair::isnull() {
	if (f==-1 && s==-1) return true;
	return false;
}
void pair::print_pair() {printf("(%d,%d) ",f,s); std::cout.flush();}

// List ------------------------------------------------------
List initList(int l, int val) {
	List L;
	L.resize(l,val);
	return (L);
}

List initList(std::vector<int> l) {
	int n(l.size());
	List L;
	L.resize(n);
	for (int i=0; i<n; i++) L[i] = l[i];
	return L;
}

List random_permut(int n) { // Return a random permutation
	List L;
	try { L.resize(n); } 
	catch(std::exception& e) {myexception(e);}
	for (int i=0; i<n; i++) L[i]=i;
	std::random_shuffle(L.begin(),L.end());
	return(L);
}

void print_list(std::string s, const List& L) {
	printf(s.c_str()); printf("\n"); int n = L.size();
	//if (n==0) { // doesn't want to be compiled... mystery
		//printf("   ! Empty list\n"); 
		//return;
	//}
	for (int i=0; i<n; i++) printf("%5d : %s \n",i,f(L[i]));
	printf("\n");
}
bool isnull(const List& L) { // Test if the list is null
	for (int i=0; i<L.size(); i++) if (L[i]!=0) return false;
	return true;
}

int sumlist(const List& L) { // Sum of the list
	int s(0);
	for (int i=0; i<L.size(); i++) s += L[i];
	return s;
}
void print_hist(List hist_petal) {
	for (int i=0;i<hist_petal.size();i++) {
		printf("%5d , ",hist_petal[i]);
		if(((i+1)/10)*10==i+1) printf("\n"); //<- ?
	}
}

// Table -----------------------------------------------------
Table initTable(int l, int c, int val) {
	Table T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,val);
	return(T);
}

std::vector<std::vector<pair> > initTable_pair(int l, int c) {
	std::vector<std::vector<pair> > T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,pair());
	return(T);
}

void print_table(std::string s, const Table& T, bool b) {
	printf(s.c_str()); if (b) printf(" with total"); printf("\n");
	int l = T.size();
	if (l==0) {
		printf("   ! Empty table");
		return;
	}
	int c = T[0].size();
	for (int i=0; i<l; i++) {
		if (T[i].size()!=c) {
			printf("   ! Table has different number of columns depending to the lines\n");
			return;
		}
	}
	std::string s0(" ",10);
	const char* space = s0.c_str();
	printf(space);
	for (int j=0; j<c; j++) printf("%11d",j);
	printf("\n");
	int begin = (isnull(initList(T[0]))) ? 1 : 0;
	for (int i=begin; i<l; i++) {
		printf("%4d",i);
		for (int j=0; j<c; j++) printf("%11s",f(T[i][j]));
		if (b) printf("%11s",f(sumlist(initList(T[i]))));
		printf("\n");
	}
	printf("\n");
}

void print_table(std::string s, const std::vector<std::vector<pair> >& T) {
	printf(s.c_str());
	int l = T.size();
	if (l==0) {
		printf("   ! Empty table");
		return;
	}
	int c = T[0].size();
	for (int i=0; i<l; i++) {
		if (T[i].size()!=c) {
			printf("   ! Table has different number of columns depending to the lines\n");
			return;
		}
	}
	std::string s0(" ",10);
	const char* space = s0.c_str();
	printf(space);
	for (int j=0; j<c; j++) printf("%11d",j);
	printf("\n");
	for (int i=0; i<l; i++) {
		printf("%4d",i);
		for (int j=0; j<c; j++) printf("%11s",p2s(T[i][j]).c_str());
		printf("\n");
	}
	printf("\n");
}

// Time ------------------------------------------------------
Time::Time() {
	s = get_time();
	e = get_time();
}

double get_time() { // Get wall time
	struct timeval time;
	if (gettimeofday(&time,NULL))
		return(0);
	else
		return((double)time.tv_sec + (double)time.tv_usec * .000001);
}

double time_diff(Time t) {
	return (t.e - t.s);
}

void init_time(Time& t, const std::string s) {
	t.s = get_time();
	if (s.size() != 0) std::cout << s << std::endl;
}

void init_time_at(Time& time, const std::string s, Time& t) {
	time.s = get_time();
	std::cout << s;
	print_time(t," at");
}

void print_time(Time& t, const std::string s) {
	t.e = get_time();
	double d = time_diff(t);
	if (d<60.) std::cout << s << " " << std::setprecision(3) << d << " s\n" << std::flush;
	else {
		int mn = floor(d/60);
		double sc = d-(mn*60);
	       	std::cout << s << " " << mn << " mn " << std::setprecision(3) << sc << " s\n" << std::flush;
	}
}

// Errors -----------------------------------------------------
void myexit(const int flag) {
	std::cout.flush();
	std::cerr.flush();
	exit(flag);
}

void myexception(const std::exception& e) {
	std::cout<<"Exception: "<<e.what()<<std::endl;
	std::cout.flush();
	std::cerr.flush();
	exit(1);
}

void myexception(std::exception& e) {
	std::cout<<"Exception: "<<e.what()<<std::endl;
	std::cout.flush();
	std::cerr.flush();
	exit(1);
}

void error(std::string s) {
	std::cout << "Error : " << s << std::endl;
	std::cout.flush();
	std::cerr.flush();
	exit(1);
}
// To String --------------------------------------------------
const char* f(int i) { // int 1526879 -> const char* 1,526,879
	std::stringstream ss;
	ss << i;
	std::string s = ss.str();
	std::string str("");
	int n = s.size();
	for (int l=0; l<n; l++) {
		str.push_back(s[l]);
		if ((n-1-l)%3==0 && n-l-1!=0) str.push_back(',');
	}
	return str.c_str();
}

std::string i2s(int i) {
	std::stringstream ss;
	ss << i;
	std::string s = ss.str();
	return s;
}

std::string p2s(int j, int k) {
	std::stringstream ss;
	ss << "(" << j << "," << k << ")";
	std::string s = ss.str();
	return s;
}

std::string p2s(pair p) {
	return(p2s(p.f,p.s));
}

std::string siz(const Table& T) {
	std::stringstream ss;
	ss << " size " << T.size() << "x" << T[0].size();
	return ss.str();
}
std::string siz(const std::vector<std::vector<pair> >& T) {
	std::stringstream ss;
	ss << " size " << T.size() << "x" << T[0].size();
	return ss.str();
}
// Other ------------------------------------------------------
void check_args(int n) { // Check if the arguments of the executable are right
	if (n != 5) {
		std::cerr << "Usage: assign <ObjFile> <PlateFile> <FiberFile> <OutFile>" << std::endl;
		myexit(1);
	}
}

int sq(int n) { return(n*n); }
double sq(double n) { return(n*n); }
double sq(double a,double b) { return(a*a + b*b); }
double percent(int a, int b) { return(a*100./b); }
