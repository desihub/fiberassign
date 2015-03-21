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

typedef std::string str;
// pair ------------------------------------------------------
class pair {
	public:
	int f, s;

	pair();
	pair(int a, int b);
	void setnull();
	bool isnull();
	void print_pair();
};

// List ------------------------------------------------------
class List : public std::vector<int> {};
class Dlist : public std::vector<double> {};
class Plist : public std::vector<pair> {};
List initList(int l, int val = 0);
List initList(std::vector<int> l);
List initList(int l[]);
Dlist initDlist(int l, double val = 0.0);
std::vector<str> initList(str l[]);
void print_list(str s, const List& L);
void print_list(str s, std::vector<str> L);
void print_list_line(const List& L);
void print_Dlist(str s, const Dlist& L);
bool isnull(const List& L); // Test if the list is null
int sumlist(const List& L); // Sum of the list
int max(const List& L);
bool isfound(int n, const List& L); // True iff i is found in L
List values(const List& L); // Different taken values
List random_permut(int n); // Return a random permutation
List random_permut(const List& L);
void print_hist(str s, int i, List hist_petal, bool latex=false);
void erase(int i, List& L); // Erase i th element
List complementary(int size, const List& L); // Complementary list (provided L is a subset of [0,size-1])
List sublist(int begin, int size, const List& L);

// Table -----------------------------------------------------
class Table : public std::vector<List> {};
class Dtable : public std::vector<Dlist> {};
class Ptable : public std::vector<Plist> {};
Table initTable(int l, int c, int val = 0); 
Table initTable(const Table& T); 
Ptable initTable_pair(int l, int c);
Dtable initTable_double(int l, int c, double val=0.0);
void verif(const Table& T); // Verifies it's square
void print_table(str s, const Table& T, bool latex=false);
void print_table(str s, const Dtable& T, bool latex=false);
void print_table(str s, const Ptable& T);
// Other constr wouldn't work if we would define a constructor in class Table
List histogram(const Table& T, int interval);
Table with_tot(const Table& T); // Add a row with sum of lines
int max_row(const Table& T); // Max number of rows
List max_on_row(const Table& T); // List of max on each row

// Cube ------------------------------------------------------
class Cube : public std::vector<Table> {};
Cube initCube(int l, int c, int d, int val = 0); 

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
void deb(int a, int b);

// To String --------------------------------------------------
str f(int i); // int 1526879 -> const char* 1,526,879
str f(double i);
str i2s(int i);
str p2s(pair p);
str p2s(int j, int k);
str siz(const Table& T);
str siz(const std::vector<std::vector<pair> >& T);
str format(int size, str s);

// Other ------------------------------------------------------
void check_args(int n); 
void print_stats(str s, int cnt, int std, int avg, int min, int max);
int sq(int n); // n*n
double sq(double n);
double sq(double a,double b); // a*a+b*b
double percent(int a, int b);
void fl();
#endif
