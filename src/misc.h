#ifndef MISC_H
#define MISC_H

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
#include        <map>

typedef std::string str;
typedef std::map<str,int> Smap;
// pair ------------------------------------------------------
class pair {
	public:
	int f, s;

	pair();
	pair(int a, int b);
	void setnull();
	bool isnull() const;
	void print_pair() const;
};

class dpair {
	public:
	double f, s;

	dpair();
	dpair(double a, double b);
	void print() const;
	bool isnull() const;
	/*dpair pair_rand(double a, double b);*/
};
dpair operator+(dpair const& a, dpair const& b);
dpair operator-(dpair const& a, dpair const& b);
dpair operator-(dpair const& a, double b);

// List ------------------------------------------------------
class List : public std::vector<int> {};
class Dlist : public std::vector<double> {};
class Plist : public std::vector<pair> {};
class Dplist : public std::vector<dpair> {};
class Slist : public std::vector<str> {};
List Null();
Slist Snull();
List initList(int l, int val = 0);
List initList(std::vector<int> l);
List initList(int l[], int size);
Dlist initDlist(int l, double val = 0.0);
Slist initList(str l[], int size);
void print_list(str s, const List& L);
void print_list(str s, Slist L);
void print_list_line(const List& L);
void print_Dlist(str s, const Dlist& L);
void print_Plist(const Plist& L, str s="");
bool isnull(const List& L); // Test if the list is null
bool isnull(const Slist& L); // Test if the list is null ("")
int sumlist(const List& L); // Sum of the list
double sumlist(const Dlist& L);
int max(const List& L);
Dlist percents(const List& L, int n); // Percents list of L/n
Dlist percents(const Dlist& L, double n);
bool isfound(int n, const List& L); // True iff i is found in L
int isfound(pair p, const Plist& L); // Return the position of where it's found, -1 if not
List values(const List& L); // Different taken values
List random_permut(int n); // Return a random permutation
List random_permut(const List& L);
void print_hist(str s, int i, List hist_petal, bool latex=false);
void erase(int i, List& L); // Erase i th element
void erase(int i, Plist& L); // Erase i th element
List complementary(int size, const List& L); // Complementary list (provided L is a subset of [0,size-1])
List sublist(int begin, int size, const List& L);
void switch_elmts(int a, int b, List& L); // Switch elements a and b
List sort(const List& L); // Sort by increasing order
List inverse(const List& L);
List cumulate(const List& L); // Cumulated sum
Dlist cumulate(const Dlist& L);
Dlist division(const Dlist& L, double d);

// Table -----------------------------------------------------
class Table : public std::vector<List> {};
class Dtable : public std::vector<Dlist> {};
class Ptable : public std::vector<Plist> {};
Table initTable(int l, int c, int val = 0); 
Table initTable(const Table& T); 
Ptable initPtable(int l, int c);
Dtable initDtable(int l, int c, double val=0.0);
void verif(const Table& T); // Verifies it's square
void print_table(str s, const Table& T, bool latex=false, Slist labels=Snull());
void print_table(str s, const Dtable& T, bool latex=false, Slist labels=Snull());
void print_table(str s, const Ptable& T);
// Other constr wouldn't work if we would define a constructor in class Table
List histogram(const Table& T, int interval);
List histogram(const List& L, int interval);
Dlist histogram(const Dlist& L, double interval);
Dlist histogram(const Dplist& L, double interval);
void print_mult_table_latex(str s, str ss, Table T, int multX=1);
void print_mult_Dtable_latex(str s, str ss, Dtable T, double multX=1);
Table with_tot(const Table& T); // Add a row with sum of lines
int max_row(const Table& T); // Max number of rows
int max_row(const Dtable& T);
List max_on_row(const Table& T); // List of max on each row
void make_square(Table& T); // Make it square, filling with zeros
Dtable divide(const Table& T, double d);
Table divide_floor(const Table& T, double d);
Dtable ddivide_floor(const Table& T, double d);
Dtable concatenate(const Dtable& T1, const Dtable& T2); // Horizontal concatenation T1:T2

// Cube ------------------------------------------------------
class Cube : public std::vector<Table> {};
Cube initCube(int l, int c, int d, int val = 0); 
class Dcube : public std::vector<Dtable> {};
Dcube initDcube(int l, int c, int d, double val = 0.0); 
int max_row(const Dcube& C); // Max size of 2nd dim

// Time ------------------------------------------------------
struct Time { // Tstart, Tend (t.s & t.e)
	double s, e;
	Time();
};
double get_time();
double time_diff(Time t);
void print_time(Time& t, const str s = "# Elapsed time");
void init_time(Time& t, const str s = "");
void init_time_at(Time& time, const str s, Time& t);

// Errors -----------------------------------------------------
void myexit(const int i);
void myexception(const std::exception& e);
void myexception(std::exception& e); // Intel compiler needs this version too ...
void error(str s);
void deb(int a);
void deb(str a);
void deb(int a, int b);
void debl(int a);
void debl(str a);
void debl(double a);
/*void debl(bool a);*/
void deb(double a);

// Conversions ------------------------------------------------
str f(int i); // int 1526879 -> const char* 1,526,879
str f(double i);
str i2s(int i);
str p2s(pair p);
str p2s(int j, int k);
int s2i(str s);
bool s2b(str s);
double s2d(str s);
str siz(const Table& T);
str siz(const Ptable& T);
str format(int size, str s);
str erase_spaces(str s); // Same str without spaces

// Parse ------------------------------------------------------
Slist s2vec(str const& s, char const delimiter);
void printFile(const char file[]);

// Other ------------------------------------------------------
void check_args(int n); 
int max(int a, int b);
int min(int a, int b);
void print_stats(str s, int cnt, int std, int avg, int min, int max);
int sq(int n); // n*n
double sq(double n);
double sq(double a,double b); // a*a+b*b
double sq(dpair p);
double sq(dpair p, dpair q); // |vec(pq)|²
double scalar_prod(dpair p, dpair q, dpair d); // vec(pq).vec(pd)
double norm(double a,double b); // sqrt(a*a+b*b)
double norm(dpair p);
double percent(int a, int b);
void fl();
double fRand(double fMin, double fMax);
dpair cartesian(double r, double theta);
dpair cartesian(dpair X);
dpair polar(dpair X);
dpair cos_sin_angle(dpair P); // Cos and sin of the angle xOP
dpair sum_angles(dpair p1, dpair p2); // cos(t+a) = cos t * cos a + ... and sin(t+a) = ...
double dist(dpair c1, dpair c2);
int sign(double a);
#endif
