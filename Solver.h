#pragma once
#include<string>

struct Perem {
	double ro;
	double u;
	double p;
};

struct consv {
	double U1;
	double U2;
	double U3;
};

struct potoc {
	double F1;
	double F2;
	double F3;
};

class Solver
{
protected:
	Perem kr_l;
	Perem kr_r;
	double x0;
	int N, K;

	double h, tau, max_T;

	double max_x;
	double min_x;

public:
	int metod;
	Solver();
	void read_file(std::string name_file);//сделать их виртуальнами
	void virtual solve() {}
	void print(std::string name_metod, double *ro, double *u, double *p, int n1);
	void convert(std::string name_metod);
	void token(std::string str);
	~Solver();
};

