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
List initList(int l, int val = 0);
List initList(std::vector<int> l);
void print_list(std::string s, const List& L);
bool isnull(const List& L); // Test if the list is null
int sumlist(const List& L); // Sum of the list
List random_permut(int n); // Return a random permutation
void print_hist(List hist_petal);

// Table -----------------------------------------------------
class Table : public std::vector<std::vector <int> > {};
Table initTable(int l, int c, int val = 0); 
void print_table(std::string s, const Table& T, bool b = false);
void print_table(std::string s, const std::vector<std::vector<pair> >& T);
// Other constr wouldn't work if we would define a constructor in class Table
std::vector<std::vector<pair> > initTable_pair(int l, int c);

// Time ------------------------------------------------------
struct Time { // Tstart, Tend (t.s & t.e)
	double s, e;
	Time();
};
double get_time();
double time_diff(Time t);
void print_time(Time& t, const std::string s = "# Elapsed time");
void init_time(Time& t, const std::string s = "");
void init_time_at(Time& time, const std::string s, Time& t);

// Errors -----------------------------------------------------
void myexit(const int i);
void myexception(const std::exception& e);
void myexception(std::exception& e); // Intel compiler needs this version too ...
void error(std::string s);

// To String --------------------------------------------------
const char* f(int i); // int 1526879 -> const char* 1,526,879
std::string i2s(int i);
std::string p2s(pair p);
std::string p2s(int j, int k);
std::string siz(const Table& T);
std::string siz(const std::vector<std::vector<pair> >& T);

// Other ------------------------------------------------------
void check_args(int n); 
int sq(int n); // n*n
double sq(double n);
double sq(double a,double b); // a*a+b*b
double percent(int a, int b);
#endif
