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
	return L;
}

List initList(std::vector<int> l) {
	int n(l.size());
	List L;
	L.resize(n);
	for (int i=0; i<n; i++) L[i] = l[i];
	return L;
}

List initList(int l[]) {
	List L;
	int size = sizeof(l);
	L.resize(size);
	for (int i=0; i<size; i++) L[i] = l[i];
	return L;
}

Dlist initDlist(int l, double val) {
	Dlist L;
	L.resize(l,val);
	return L;
}

std::vector<str> initList(str l[]) {
	std::vector<str> L;
	int size = sizeof(l);
	L.resize(size);
	for (int i=0; i<size; i++) L[i] = l[i];
	return L;
}

List random_permut(int n) {
	List L;
	try { L.resize(n); } 
	catch(std::exception& e) {myexception(e);}
	for (int i=0; i<n; i++) L[i]=i;
	std::random_shuffle(L.begin(),L.end());
	return L;
}

List random_permut(const List& L) {
	int n = L.size();
	List perm = random_permut(n);
	List l;
	try { l.resize(n); } 
	catch(std::exception& e) {myexception(e);}
	for (int i=0; i<n; i++) l[i] = L[perm[i]];
	return l;
}

void print_list(str s, const List& L) {
	printf("%s \n",s.c_str());
	int n = L.size();
	//if (n==0) { // doesn't want to be compiled... mystery
		//printf("   ! Empty list\n"); 
		//return;
	//}
	for (int i=0; i<n; i++) printf("%5d : %s \n",i,f(L[i]).c_str());
	printf("\n");
}

void print_list(str s, std::vector<str> L) {
	printf(s.c_str()); printf("\n"); int n = L.size();
	//if (n==0) { // doesn't want to be compiled... mystery
		//printf("   ! Empty list\n"); 
		//return;
	//}
	for (int i=0; i<n; i++) printf("%5d : %s \n",i,L[i].c_str());
	printf("\n");
}

void print_list_line(const List& L) {
	printf("(");
	for (int i=0; i<L.size(); i++) {
		str s = i==L.size()-1 ? "" : ",";
		printf("%d%s",L[i],s.c_str());
	}
	printf(") "); std::cout.flush();
}

void print_Dlist(str s, const Dlist& L) {
	printf("%s \n",s.c_str());
	int n = L.size();
	for (int i=0; i<n; i++) printf("%5d : %.3f\n",i,L[i]);
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

int max(const List& L) {
	int m(-1e8);
	for (int i=0; i<L.size(); i++) if (L[i]>m) m = L[i];
	return m;
}

bool isfound(int n, const List& L) {
	if (L.size()==0) return false;
	for (int i=0; i<L.size(); i++) if (L[i]==n) return true;
	return false;
}

List values(const List& L) {
	List l;
	for (int i=0; i<L.size(); i++) {
		int n = L[i];
		if (!isfound(n,l)) l.push_back(n);
	}
	return l;
}

void print_hist(str s, int i, List hist_petal, bool latex) {
	str et = latex ? " & " : " | ";
	str slash = latex ? " \\\\ \n" : "\n";

	str rrr(10,'r');
	if (latex) printf("\\begin{table}[h]\\begin{center} \n \\caption{%s (interval %d)} \n \\begin{tabular}{%s}",s.c_str(),i,rrr.c_str());
	else printf("%s (interval %d)\n",s.c_str(),i);

	int size = hist_petal.size();
	for (int i=0; i<size; i++) {
		str s = ((i+1)%10==0 || i==size-1) ? slash : et;
		printf("%5d %s",hist_petal[i],s.c_str());
	}
	if (latex) printf("\\end{tabular}\\end{center}\\end{table}");
	std::cout << std::endl;
}

void erase(int i, List& L) { L.erase (L.begin()+i); }

List complementary(int size, const List& L) {
	List l = initList(size);
	for (int i=0; i<size; i++) l.push_back(i);
	for (int i=0; i<L.size(); i++) erase(L[i],l);
	return l;
}

List sublist(int begin, int size, const List& L) {
	List l = initList(size);
	for (int i=0; i<size; i++) l[i] = L[begin+i];
	return l;
}

// Table -----------------------------------------------------
Table initTable(int l, int c, int val) {
	Table T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,val);
	return T;
}

Ptable initTable_pair(int l, int c) {
	Ptable T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,pair());
	return T;
}

Dtable initTable_double(int l, int c, double val) {
	Dtable T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,val);
	return T;
}

void verif(const Table& T) {
	int l = T.size();
	if (l==0) {
		printf("   ! Empty table\n");
		return;
	}
	int c = T[0].size();
	for (int i=0; i<l; i++) {
		if (T[i].size()!=c) {
			printf("   ! Table has different number of columns depending to the lines\n");
			return;
		}
	}
}

int max_row(const Table& T) {
	int max(0);
	for (int i=0; i<T.size(); i++) {
		int c = T[i].size();
		if (c>max) max = c;
	}
	return max;
}

void print_table(str s, const Table& T, bool latex) {
	int SIZE = 7;
	int MAXSIZE = 10;

	bool square(true);
	int l = T.size();
	int cc = T[0].size();
	for (int i=0; i<l; i++) if (T[i].size()!=cc) square = false;

	List maxRow = initList(cc);
	if (square) {
		Table sizestr = initTable(l,cc);
		for (int i=0; i<l; i++) {
			for (int j=0; j<cc; j++) {
				int sz = f(T[i][j]).size()+1;
				sizestr[i][j] = sz<MAXSIZE ? sz : MAXSIZE;
			}
		}
		maxRow = max_on_row(sizestr);
	}
	else maxRow = initList(cc,7);

	str et = latex ? " &" : "|";
	str slash = latex ? " \\\\ \n" : "\n";

	str rrr(T[0].size(),'r');
	if (latex) printf("\n\\begin{table}[h]\\begin{center} \n\\caption{%s} \n\\begin{tabular}{%s}",s.c_str(),rrr.c_str());
	else printf(s.c_str());
	printf("\n");

	str space(6,' ');
	printf(space.c_str());
	int max = max_row(T);
	for (int j=0; j<max; j++) {
		str s = (j==max-1) ? slash : et;
		printf("%s%s",format(maxRow[j],f(j)).c_str(),s.c_str());
	}
	for (int i=0; i<l && i<10; i++) {
		printf("%3d  %s",i,et.c_str());
		int c = T[i].size();
		for (int j=0; j<c; j++) {
			str s = (j==c-1) ? slash : et;
			str num = T[i][j]==-1 ? "" : f(T[i][j]).c_str();
			printf("%s%s",format(maxRow[j],num).c_str(),s.c_str());
		}
	}
	if (latex) printf("\\end{tabular}\\end{center}\\end{table}");
	printf("\n");
}

void print_table(str s, const Dtable& T, bool latex) {
	//verif(T);
	int l = T.size();
	int c = T[0].size();

	str et = latex ? " & " : " | ";
	str slash = latex ? " \\\\ \n" : "\n";

	str rrr(T[0].size(),'r');
	if (latex) printf("\n\\begin{table}[h]\\begin{center} \n\\caption{%s} \n\\begin{tabular}{%s}",s.c_str(),rrr.c_str());
	else printf(s.c_str());
	printf("\n");

	str space(" ",10);
	printf(space.c_str());
	for (int j=0; j<c; j++) {
		str s = (j==c-1) ? slash : et;
		printf("%11d %s",j,s.c_str());
	}
	for (int i=0; i<l; i++) {
		printf("%4d %s",i,et.c_str());
		for (int j=0; j<c; j++) {
			str s = (j==c-1) ? slash : et;
			printf("%11s %s",f(T[i][j]).c_str(),s.c_str());
		}
	}
	if (latex) printf("\\end{tabular}\\end{center}\\end{table}");
	printf("\n");
}

void print_table(str s, const Ptable& T) {
	printf(s.c_str()); printf("\n");
	int l = T.size();
	if (l==0) {
		printf("   ! Empty table\n");
		return;
	}
	int c = T[0].size();
	for (int i=0; i<l; i++) {
		if (T[i].size()!=c) {
			printf("   ! Table has different number of columns depending to the lines\n");
			return;
		}
	}
	str s0(" ",10);
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

List histogram(const Table& T, int interval) {
	List hist; int l = T.size(); int c = T[0].size();
	for(int i=0; i<l; i++) { 
		for(int j=0; j<c; j++) {
			int n = floor(T[i][j]/interval);
			if (n>=hist.size()) { hist.resize(n+1); hist[n] = 0;}
			hist[n]++;
		}
	}
	return hist;
}

Table with_tot(const Table& T) {
	Table M = T;
	for (int i=0; i<M.size(); i++) M[i].push_back(sumlist(M[i]));
	return M;
}

List max_on_row(const Table& T) {
	verif(T);
	int l = T.size();
	int c = T[0].size();
	List L = initList(c,-1e7);
	for (int i=0; i<l; i++) {
		for (int j=0; j<c; j++) {
			int a = T[i][j];
			if (a>L[j]) L[j] = a;
		}
	}
	return L;
}

// Cube ------------------------------------------------------
Cube initCube(int l, int c, int d, int val) {
	Cube C;
	C.resize(l);
	for (int i=0; i<l; i++) {
		C[i].resize(c);
		for (int j=0; j<c; j++) C[i][j].resize(d,val);
	}
	return C;
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

void init_time(Time& t, const str s) {
	t.s = get_time();
	if (s.size() != 0) std::cout << s << std::endl;
}

void init_time_at(Time& time, const str s, Time& t) {
	time.s = get_time();
	std::cout << s;
	print_time(t," at");
}

void print_time(Time& t, const str s) {
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

void error(str s) {
	std::cout << "Error : " << s << std::endl;
	std::cout.flush();
	std::cerr.flush();
	exit(1);
}

void deb(int a) { // debug
	std::cout << " " << a << " " << std::endl;
}

void deb(int a, int b) { // debug
	std::cout << " " << a << "," << b << " " << std::endl;
}

// To String --------------------------------------------------
str f(int i) { // int 1526879 -> const char* 1,526,879
	std::stringstream ss;
	ss << i;
	str s = ss.str();
	str str("");
	int n = s.size();
	for (int l=0; l<n; l++) {
		str.push_back(s[l]);
		if ((n-1-l)%3==0 && n-l-1!=0) str.push_back(',');
	}
	return str.c_str();
}

str f(double i) { // Be careful : double 1523.5412 -> 1,523.54 (max 6 numbers)
	std::stringstream ss;
	ss << i;
	str s = ss.str();
	str str("");
	int n = s.size();
	bool dot(true);
	int dec(0);
	int dot_pos(-1);
	for (int a=0; a<n; a++) if (s[a]=='.') dot_pos = a;
	int bef(n-dot_pos);
	for (int l=0; l<n; l++) {
		if (s[l]=='.') dot = false;
		if (dot) str.push_back(s[l]);
		else if (dec <= 3) { dec++; str.push_back(s[l]); }
		if ((n-bef-1-l)%3==0 && n-bef-l-1!=0 && dot) str.push_back(',');
	}
	return str.c_str();
}

str i2s(int i) {
	std::stringstream ss;
	ss << i;
	str s = ss.str();
	return s;
}

str p2s(int j, int k) {
	std::stringstream ss;
	ss << "(" << j << "," << k << ")";
	str s = ss.str();
	return s;
}

str p2s(pair p) {
	return(p2s(p.f,p.s));
}

str siz(const Table& T) {
	std::stringstream ss;
	ss << " size " << T.size() << "x" << T[0].size();
	return ss.str();
}

str siz(const std::vector<std::vector<pair> >& T) {
	std::stringstream ss;
	ss << " size " << T.size() << "x" << T[0].size();
	return ss.str();
}

str format(int size, str s) {
	str S = s; str sp = "";
	int n = s.size();
	if (size>n) for (int h=0; h<size-n; h++) sp.push_back(' ');
	return (sp+s);
}

// Other ------------------------------------------------------
void check_args(int n) { // Check if the arguments of the executable are right
	if (n != 5) {
		std::cerr << "Usage: assign <ObjFile> <PlateFile> <FiberFile> <OutFile>" << std::endl;
		myexit(1);
	}
}

void print_stats(str s, int cnt, int avg, int std, int min, int max) {
	double avg1 = ((double)avg)/((double) cnt);
	double std1 = sqrt(((double) std)/((double) cnt)-sq(avg1)); 
	printf("%s %.1f +/- %.1f [%d,%d]",s.c_str(),avg1,std1,min,max);
	std::cout << std::endl;
}

int sq(int n) { return(n*n); }
double sq(double n) { return(n*n); }
double sq(double a, double b) { return(a*a + b*b); }
double percent(int a, int b) { return(a*100./b); }
