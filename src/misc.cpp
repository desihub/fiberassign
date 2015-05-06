#include	<cstdlib>
#include	<cmath>
#include        <math.h>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include        <stdio.h>
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
bool pair::isnull() const {
	if (f==-1 && s==-1) return true;
	return false;
}
void pair::print_pair() const {printf("(%d,%d) ",f,s); fl();}

// dpair ------------------------------------------------------
dpair::dpair() {f = 0; s = 0;}
dpair::dpair(double a, double b) {f = a; s = b;}
dpair operator+(dpair const& a, dpair const& b) {
    return dpair(a.f+b.f,a.s+b.s);
}
dpair operator-(dpair const& a, dpair const& b) {
    return dpair(a.f-b.f,a.s-b.s);
}
dpair operator-(dpair const& a, double b) {
    return dpair(a.f-b,a.s-b);
}
void dpair::print() const {printf("(%f,%f) ",f,s); fl();}
//void dpair::pair_rand(double a, double b) { return dpair(1,2); }

// List ------------------------------------------------------
List Null() { List L; return L; }

Slist Snull() { Slist L; return L; }

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

List initList(int l[], int size) {
	List L;
	L.resize(size);
	for (int i=0; i<size; i++) L[i] = l[i];
	return L;
}

Dlist initDlist(int l, double val) {
	Dlist L;
	L.resize(l,val);
	return L;
}

Slist initList(str l[], int size) {
	Slist L;
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

void print_list(str s, Slist L) {
	printf("%s \n",s.c_str());
	int n = L.size();
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

void print_Plist(const Plist& L, str s) {
	printf("%s \n",s.c_str());
	int n = L.size();
	for (int i=0; i<n; i++) printf("%5d : (%d,%d)\n",i,L[i].f,L[i].s);
	printf("\n");
}

bool isnull(const List& L) { // Test if the list is null
	for (int i=0; i<L.size(); i++) if (L[i]!=0) return false;
	return true;
}

bool isnull(const Slist& L) {
	for (int i=0; i<L.size(); i++) if (L[i]!="") return false;
	return true;
}

int sumlist(const List& L) { // Sum of the list
	int s(0);
	for (int i=0; i<L.size(); i++) s += L[i];
	return s;
}

double sumlist(const Dlist& L) {
	double d = 0;
	for (int l=0; l<L.size(); l++) d += L[l];
	return d;
}

int max(const List& L) {
	int m(-1e8);
	for (int i=0; i<L.size(); i++) if (L[i]>m) m = L[i];
	return m;
}

Dlist percents(const List& L, int n) {
	Dlist Dl;
	for (int i=0; i<L.size(); i++) {
		Dl.push_back(percent(L[i],n));
	}
	return Dl;
}

Dlist percents(const Dlist& L, double n) {
	Dlist Dl;
	for (int i=0; i<L.size(); i++) {
		Dl.push_back(percent(L[i],n));
	}
	return Dl;
}

bool isfound(int n, const List& L) {
	if (L.size()==0) return false;
	for (int i=0; i<L.size(); i++) if (L[i]==n) return true;
	return false;
}

int isfound(pair p, const Plist& L) {
	int j = p.f; int k = p.s;
	if (L.size()==0) return -1;
	for (int i=0; i<L.size(); i++) if (L[i].f==j && L[i].s==k) return i;
	return -1;
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
	str slash = latex ? "\\\\ \n" : "\n";

	str rrr(10,'r');
	if (latex) printf("\\begin{table}[H]\\begin{center} \n \\caption{%s (interval %d)} \n \\begin{tabular}{%s}\n",s.c_str(),i,rrr.c_str());
	else printf("%s (interval %d)\n",s.c_str(),i);

	int size = hist_petal.size();
	for (int i=0; i<size; i++) {
		str s = ((i+1)%10==0 || i==size-1) ? slash : et;
		printf("%4s %s",f(hist_petal[i]).c_str(),s.c_str());
	}
	if (latex) printf("\\end{tabular}\\end{center}\\end{table} \n");
	fl();
}

void print_mult_table_latex(str s, str ss, Table T, int multX) {
	printf("# %s \n  Output data in %s \n\n",s.c_str(),ss.c_str());
	FILE* pFile;
	pFile = fopen(ss.c_str(),"w");
	if (!pFile) {
		printf("File opening failed \n"); fl();
		myexit(1);
	}

	fprintf(pFile,"x ");
	for (int i=0; i<T.size(); i++) fprintf(pFile,"%d ",i);
	fprintf(pFile,"\n");
	int maxrow = max_row(T);
	for (int j=0; j<maxrow; j++) {
		fprintf(pFile,"%d ",j*multX);
		for (int i=0; i<T.size(); i++) {
			if (j<T[i].size()) fprintf(pFile,"%d ",T[i][j]);
			else fprintf(pFile,"0 ");
		}
		fprintf(pFile,"\n");
	}
	fclose(pFile);
}

void print_mult_Dtable_latex(str s, str ss, Dtable T, double multX) {
	printf("# %s \n  Output data in %s \n\n",s.c_str(),ss.c_str());
	FILE* pFile;
	pFile = fopen(ss.c_str(),"w");
	if (!pFile) {
		printf("File opening failed \n"); fl();
		myexit(1);
	}

	fprintf(pFile,"x ");
	for (int i=0; i<T.size(); i++) fprintf(pFile,"%d ",i);
	fprintf(pFile,"\n");
	int maxrow = max_row(T);
	for (int j=0; j<maxrow; j++) {
		fprintf(pFile,"%f ",j*multX);
		for (int i=0; i<T.size(); i++) {
			if (j<T[i].size()) fprintf(pFile,"%f ",T[i][j]);
			else fprintf(pFile,"0 ");
		}
		fprintf(pFile,"\n");
	}
	fclose(pFile);
}

void erase(int i, List& L) { L.erase (L.begin()+i); }

void erase(int i, Plist& L) { L.erase (L.begin()+i); }

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

void switch_elmts(int a, int b, List& L) {
	if (L.size()<=a || L.size()<=b) { printf("Error in switch elmts\n"); fl(); }
	int tmp = L[a];
	L[a] = L[b];
	L[b] = tmp;
}

List sort(const List& L) {
	List l = L;
	std::sort(l.begin(), l.end());  
	return l;
}

List inverse(const List& L) {
	List l = initList(max(L));
	for (int i=0; i<L.size(); i++) l[L[i]] = i;
	return l;
}

List cumulate(const List& L) {
	List l = initList(L.size());
	int a = 0;
	for (int i=0; i<L.size(); i++) {
		a += L[i];
		l[i] = a;
	}
	return l;
}

Dlist division(const Dlist& L, double d) {
	Dlist Dl;
	for (int i=0; i<L.size(); i++) Dl.push_back(L[i]/d);
	return Dl;
}

// Table -----------------------------------------------------
Table initTable(int l, int c, int val) {
	Table T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,val);
	return T;
}

Ptable initPtable(int l, int c) {
	Ptable T;
	T.resize(l);
	for (int i=0; i<l; i++) T[i].resize(c,pair());
	return T;
}

Dtable initDtable(int l, int c, double val) {
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

int max_row(const Dtable& T) {
	int max(0);
	for (int i=0; i<T.size(); i++) {
		int c = T[i].size();
		if (c>max) max = c;
	}
	return max;
}


void print_table(str s, const Table& T, bool latex, Slist labels) {
	int SIZE = 7;
	int MAXSIZE = 10;
	int MAXLINE = 50;
	int STANDARDSIZE = 1;

	bool square(true);
	int l = T.size();
	if (MAXLINE<l) printf("Cautious, table is of size %d, only %d will be displayed \n",l,MAXLINE);
	int cc = T[0].size();
	for (int i=0; i<l; i++) if (T[i].size()!=cc) square = false;

	int max = max_row(T);
	List maxRow = initList(max);
	if (square) {
		Table sizestr = initTable(l,max);
		for (int i=0; i<l; i++) {
			for (int j=0; j<max; j++) {
				int sz = f(T[i][j]).size()+1;
				sizestr[i][j] = sz<MAXSIZE ? sz : MAXSIZE;
			}
		}
		maxRow = max_on_row(sizestr);
	}
	else maxRow = initList(max,STANDARDSIZE);

	str et = latex ? " &" : "|";
	str slash = latex ? " \\\\ \n" : "\n";

	str rrr(T[0].size(),'r');
	printf(s.c_str());
	printf("\n");

	bool labelsB = isnull(labels);
	str space(6,' ');
	printf(space.c_str());
	for (int j=0; j<max; j++) {
		str s = (j==max-1) ? slash : et;
		printf("%s%s",format(maxRow[j],f(j)).c_str(),s.c_str());
	}
	for (int i=0; i<l && i<MAXLINE; i++) {
		str pre = labelsB ? f(i) : format(11,labels[i]);
		printf("%3s  %s",pre.c_str(),et.c_str());
		int c = T[i].size();
		if (c==0) printf("\n");
		for (int j=0; j<c; j++) {
			str s = (j==c-1) ? slash : et;
			str num = T[i][j]==-1 ? "" : f(T[i][j]).c_str();
			printf("%s%s",format(maxRow[j],num).c_str(),s.c_str());
		}
	}
	printf("\n");
}

void print_table(str s, const Dtable& T, bool latex, Slist labels) {
	int SIZE = 7;
	int MAXSIZE = 10;
	int MAXLINE = 50;
	int STANDARDSIZE = 1;

	bool square(true);
	int l = T.size();
	if (MAXLINE<l) printf("Cautious, table is of size %d, only %d will be displayed \n",l,MAXLINE);
	int cc = T[0].size();
	for (int i=0; i<l; i++) if (T[i].size()!=cc) square = false;

	int max = max_row(T);
	List maxRow = initList(max);
	if (square) {
		Table sizestr = initTable(l,max);
		for (int i=0; i<l; i++) {
			for (int j=0; j<max; j++) {
				int sz = f(T[i][j]).size()+1;
				sizestr[i][j] = sz<MAXSIZE ? sz : MAXSIZE;
			}
		}
		maxRow = max_on_row(sizestr);
	}
	else maxRow = initList(max,STANDARDSIZE);

	str et = latex ? " &" : "|";
	str slash = latex ? " \\\\ \n" : "\n";

	str rrr(T[0].size(),'r');
	printf(s.c_str());
	printf("\n");

	bool labelsB = isnull(labels);
	str space(6,' ');
	printf(space.c_str());
	for (int j=0; j<max; j++) {
		str s = (j==max-1) ? slash : et;
		printf("%s%s",format(maxRow[j],f(j)).c_str(),s.c_str());
	}
	for (int i=0; i<l && i<MAXLINE; i++) {
		str pre = labelsB ? f(i) : format(11,labels[i]);
		printf("%3s  %s",pre.c_str(),et.c_str());
		int c = T[i].size();
		if (c==0) printf("\n");
		for (int j=0; j<c; j++) {
			str s = (j==c-1) ? slash : et;
			str num = T[i][j]==-1 ? "" : f(T[i][j]).c_str();
			printf("%s%s",format(maxRow[j],num).c_str(),s.c_str());
		}
	}
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
	List hist; int l = T.size();
	for(int i=0; i<l; i++) {
		int c = T[i].size();
		for(int j=0; j<c; j++) {
			int a = T[i][j];
			if (a!=-1) {
				if (a<0) { printf("Error in print_hist, neg\n"); fl(); }
				int n = floor(a/interval);
				//printf("%d %d %d %d-",n,a,i,j);
				if (n>=hist.size()) { hist.resize(n+1); hist[n] = 0;}
				hist[n]++;
			}
		}
	}
	return hist;
}

List histogram(const List& L, int interval) {
	List hist; int l = L.size();
	for(int i=0; i<l; i++) {
		int a = L[i];
		if (a!=-1) {
			if (a<0) { printf("Error in print_hist, neg\n"); fl(); }
			int n = floor(a/interval);
			if (n>=hist.size()) { hist.resize(n+1); hist[n] = 0;}
			hist[n]++;
		}
	}
	return hist;
}

Dlist histogram(const Dlist& L, double interval) {
	Dlist hist; int l = L.size();
	for(int i=0; i<l; i++) {
		double a = L[i];
		if (a!=-1) {
			if (a<0) { printf("Error in print_hist, neg\n"); fl(); }
			int n = floor(a/interval);
			if (n>=hist.size()) { hist.resize(n+1); hist[n] = 0;}
			hist[n]++;
		}
	}
	return hist;
}

Dlist histogram(const Dplist& L, double interval) {
	Dlist hist; int l = L.size();
	for(int i=0; i<l; i++) {
		double a = L[i].f;
		if (a!=-1) {
			if (a<0) { printf("Error in print_hist, neg\n"); fl(); }
			int n = floor(a/interval);
			if (n>=hist.size()) { hist.resize(n+1); hist[n] = 0;}
			hist[n] += L[i].s;
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

void make_square(Table& T) {
	int max = max_row(T);
	for (int i=0; i<T.size(); i++) {
		if (T[i].size()<max) T[i].resize(max);
	}
}

Dtable divide(const Table& T, double d) {
	Dtable t;
	for (int i=0; i<T.size(); i++) {
		Dlist l;
		for (int j=0; j<T[i].size(); j++) {
			l.push_back(T[i][j]/d);
		}
		t.push_back(l);
	}
	return t;
}

Table divide_floor(const Table& T, double d) {
	Table t;
	for (int i=0; i<T.size(); i++) {
		List l;
		for (int j=0; j<T[i].size(); j++) {
			l.push_back(T[i][j]/d);
		}
		t.push_back(l);
	}
	return t;
}

Dtable ddivide_floor(const Table& T, double d) {
	Dtable t;
	for (int i=0; i<T.size(); i++) {
		Dlist l;
		for (int j=0; j<T[i].size(); j++) {
			l.push_back(floor(T[i][j]/d));
		}
		t.push_back(l);
	}
	return t;
}

Dtable concatenate(const Dtable& T1, const Dtable& T2) {
	if (T1.size()!=T2.size()) {
		printf("Error in concatenation : different sizes \n");
		myexit(1);
	}
	Dtable T3 = T1;
	for (int i=0; i<T1.size(); i++) {
		for (int j=0; j<T2[i].size(); j++) {
			T3[i].push_back(T2[i][j]);
		}
	}
	return T3;
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

Dcube initDcube(int l, int c, int d, double val) {
	Dcube C;
	C.resize(l);
	for (int i=0; i<l; i++) {
		C[i].resize(c);
		for (int j=0; j<c; j++) C[i][j].resize(d,val);
	}
	return C;
}

int max_row(const Dcube& C) {
	int max(0);
	for (int i=0; i<C.size(); i++) {
		int c = C[i].size();
		if (c>max) max = c;
	}
	return max;
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
	std::cout << " " << a << std::endl;
}

void deb(int a, int b) { // debug
	std::cout << " " << a << "," << b << " " << std::endl;
}

void debl(int a) { // debug
	std::cout << " " << a << std::flush;
}

void deb(double a) { // debug
	std::cout << " " << a << std::endl;
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
	//for (int a=0; a<n; a++) printf("%c",s[a]);
	for (int a=0; a<n; a++) if (s[a]=='.') dot_pos = a;
	if (dot_pos==-1) dot_pos = n;
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

str siz(const Ptable& T) {
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

str erase_spaces(str s) {
	str ss = "";
	for (int h=0; h<s.size(); h++) if (s[h]!=' ') ss.push_back(s[h]);
	return ss;
}

// Other ------------------------------------------------------
void check_args(int n) { // Check if the arguments of the executable are right
	if (n != 7) {
		std::cerr << "Usage: assign <ObjFile> <PlateFile> <FiberFile> <OutFile> <Modulo galaxies> <Modulo fibers>" << std::endl;
		myexit(1);
	}
}

int max(int a, int b) {
	if (a>b) return a;
	return b;
}

int min(int a, int b) {
	if (a>b) return b;
	return a;
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
double sq(dpair p) { return p.f*p.f + p.s*p.s; }
double sq(dpair p, dpair q) { return sq(p.f-q.f) + sq(p.s-q.s); }
double scalar_prod(dpair p, dpair q, dpair d) {
	return (q.f-p.f)*(d.f-p.f)+(q.s-p.s)*(d.s-p.s);
}
double norm(double a, double b) { return(sqrt(a*a + b*b)); }
double norm(dpair p) { return(sqrt(p.f*p.f + p.s*p.s)); }
double percent(int a, int b) { return(a*100./b); }
void fl() { std::cout.flush(); }

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

dpair cartesian(double r, double theta) {
	return dpair(r*cos(theta),r*sin(theta));
}

dpair cartesian(dpair X) {
	return dpair(X.f*cos(X.s),X.f*sin(X.s));
}

dpair polar(dpair X) {
	return dpair(norm(X.f,X.s),atan(X.s/X.f));
}

dpair sum_angles(dpair t, dpair a) {
	return dpair(t.f*a.f-t.s*a.s,t.s*a.f+t.f*a.s);
}

dpair cos_sin_angle(dpair P) {
	double r_inv = 1/norm(P);
	return dpair(P.f*r_inv,P.s*r_inv);
}

double dist(dpair c1, dpair c2) {
	return norm(c1.f-c2.f,c1.s-c2.s);
}
