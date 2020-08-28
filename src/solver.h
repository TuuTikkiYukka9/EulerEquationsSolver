#pragma once
#include<string>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

struct Variables {
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

struct EquationSystem
{
	Variables leftBoundaryCondition;
	Variables rightBoundaryCondition;
	double x0;
	double maxX;
	double minX;
};

class Solver
{
protected:
	Variables kr_l;
	Variables kr_r;
	double x0;
	int N, K;

	double h, tau, max_T;

	double max_x;
	double min_x;

public:
	void virtual solve() {}
	void init(EquationSystem &eq);
};

