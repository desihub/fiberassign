#ifndef MISC_H
#define MISC_H

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
};

// Add coords of 2 points
dpair operator+(dpair const& a, dpair const& b);

dpair operator-(dpair const& a, dpair const& b);

dpair operator-(dpair const& a, double b);

// Gives cartesian coordinates of a point defined by polar coords
dpair cartesian(double r, double theta);

// Gives cartesian coordinates of a point defined by polar coords (in a dpair)
dpair cartesian(dpair X);

// Cos and sin of the angle xOP
dpair cos_sin_angle(dpair P);

// cos(t+a) = cos t * cos a + ... and sin(t+a) = ...
dpair sum_angles(dpair p1, dpair p2);

// Distance between 2 points
double dist(dpair c1, dpair c2);

double norm(dpair p);

// Norm squared of a point
double sq(dpair p);

// |AB|²
double sq(dpair A, dpair B);

// vec(AB).vec(AC)
double scalar_prod(dpair A, dpair B, dpair C);

// List ------------------------------------------------------

class List : public std::vector<int> {};

class Dlist : public std::vector<double> {};

class Plist : public std::vector<pair> {};

class Dplist : public std::vector<dpair> {};

class Slist : public std::vector<str> {};

List Null();

Slist Snull();

// Initialize a list of size l, with values val
List initList(int l, int val = 0);

Dlist initDlist(int l, double val = 0.0);

// Transform l into the List type
List initList(std::vector<int> l);

// Initialize a list from array l which is of size size
List initList(int l[], int size);

Slist initList(str l[], int size);

// Print list L, printing the message s before
void print_list(str s, const List& L);

void print_list(str s, Slist L);

void print_Dlist(str s, const Dlist& L);

void print_Plist(const Plist& L, str s="");

// Print on a line
void print_list_line(const List& L);

// Test if the list is null
bool isnull(const List& L);

// Test if the list is null ("")
bool isnull(const Slist& L);

// Sum of the list
int sumlist(const List& L);

double sumlist(const Dlist& L);

// Max of the list
int max(const List& L);

// List of 100*L[i]/n
Dlist percents(const List& L, int n);

Dlist percents(const Dlist& L, double n);

// True iff n is found in L
bool isfound(int n, const List& L);

// Return the position of where p is found, -1 if not
int isfound_pos(double n, const Dlist& L);

// Return the position of where p is found, -1 if not
int isfound(pair p, const Plist& L);

bool isfound(str n, const Slist& L);

// List of different taken values (ex 1,1,5,2,5 gives 1,5,2)
List values(const List& L);

// Return a random permutation of [0,n-1]
List random_permut(int n);

// Random permut of L
List random_permut(const List& L);

// Print histogram hist with interval i, printing s before, in latex format if wished
void print_hist(str s, int i, List hist, bool latex=false);

// Erase i th element
void erase(int i, List& L);

// Erase i th element
void erase(int i, Plist& L);

 // Complementary list (provided L is a subset of [0,size-1])
List complementary(int size, const List& L);

// Return sublist
List sublist(int begin, int size, const List& L);

// Switch elements a and b in L
void switch_elmts(int a, int b, List& L);

// Sort L by increasing order
List sort(const List& L);

// Gives the inverse map defined by L
List inverse(const List& L);

// Cumulated sum of L (integral)
List cumulate(const List& L);

Dlist cumulate(const Dlist& L);

// L[i]/d
Dlist division(const Dlist& L, double d);

// Add elements of l to L
void addlist(List& L, const List& l);

// L but with slight variations of elements
Dlist var(const List& L);

// Get permutation that leads to sort this list
List get_permut_sort(const List& l);


// Table -----------------------------------------------------

class Table : public std::vector<List> {};

class Dtable : public std::vector<Dlist> {};

class Ptable : public std::vector<Plist> {};

Table initTable(int l, int c, int val = 0);

Table initTable(const Table& T);

Ptable initPtable(int l, int c);

Dtable initDtable(int l, int c, double val=0.0);

// Verifies it's square
void verif(const Table& T);

void print_table(str s, const Table& T, bool latex=false, Slist labels=Snull());

void print_table(str s, const Dtable& T, bool latex=false, Slist labels=Snull());

void print_table(str s, const Ptable& T);

// Other constr wouldn't work if we would define a constructor in class Table

// Build histogram of interval interval from values in T
List histogram(const Table& T, int interval);

// Build histogram of interval interval from values in L
List histogram(const List& L, int interval);

Dlist histogram(const Dlist& L, double interval);

Dlist histogram(const Dplist& L, double interval);

// Print, in latex format, data from table T in file ss, writing s on cout
void print_mult_table_latex(str s, str ss, Table T, int multX=1);

void print_mult_Dtable_latex(str s, str ss, Dtable T, double multX=1);

// Add a row with sum of lines
Table with_tot(const Table& T);

// Max number of rows
int max_row(const Table& T);

int max_row(const Dtable& T);

// List of max on each row
List max_on_row(const Table& T);

// Make it square, filling with zeros
void make_square(Table& T);

Dtable divide(const Table& T, double d);

Dtable divide(const Dtable& T, double d);

Dtable mult(const Dtable& T, double d);

// floor(T[i][j]/d)
Table divide_floor(const Table& T, double d);

// floor(T[i][j]/d)
Dtable ddivide_floor(const Table& T, double d);

// Horizontal concatenation T1:T2
Dtable concatenate(const Dtable& T1, const Dtable& T2);


// Cube ------------------------------------------------------

class Cube : public std::vector<Table> {};

Cube initCube(int l, int c, int d, int val = 0);

class Dcube : public std::vector<Dtable> {};

Dcube initDcube(int l, int c, int d, double val = 0.0);

// Max size of 2nd dim
int max_row(const Dcube& C);


// Time ------------------------------------------------------
struct Time {
    // Tstart, Tend (t.s & t.e)
    double s, e;
    Time();
};

// Time from the hardware
double get_time();

// Difference e-s
double time_diff(Time t);

// Update e and print elapsed time
void print_time(Time& t, const str s = "# Elapsed time");

// Init t (init s)
void init_time(Time& t, const str s = "");

// Init time, printing global time t
void init_time_at(Time& time, const str s, Time& t);


// Errors -----------------------------------------------------

void myexit(const int i);

void myexception(const std::exception& e);

// Intel compiler needs this version too
void myexception(std::exception& e);  ...

void error(str s);

// Used for debug, just prints a
void deb(int a);

void deb(str a);

void deb(int a, int b);

void deb(double a);

// Debug : debl prints without skipping a line
void debl(int a);

void debl(str a);

void debl(double a);


// Conversions ------------------------------------------------

// int 1526879 -> str 1,526,879
str f(int i);

// double 1523.5412 -> str 1,523.54
str f(double i);

// int to string
str i2s(int i);

str d2s(double i);

str p2s(pair p);

str p2s(int j, int k);

int s2i(str s);

bool s2b(str s);

double s2d(str s);

// like in prinft("%5s") where 5 = size
str format(int size, str s);

// Same str without spaces
str erase_spaces(str s);


// Parse ------------------------------------------------------

 // "a b cd" gives vec of "a"; "b"; "cd" if delimiter = " "
Slist s2vec(str const& s, char const delimiter);

// Print the file
void printFile(const char file[]);


// Other ------------------------------------------------------

// Check number of arguments
void check_args(int n);

int max(int a, int b);

int min(int a, int b);

// n*n
int sq(int n);

double sq(double n);

// a*a+b*b
double sq(double a,double b);

// sqrt(a*a+b*b)
double norm(double a,double b);

double percent(int a, int b);

// cout.flush()
void fl();

// Generates a random number between fMin and fMax
double fRand(double fMin, double fMax);

#endif
