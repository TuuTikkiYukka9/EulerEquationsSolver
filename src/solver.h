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

struct ConservativeVariables {
	double U1;
	double U2;
	double U3;
};

struct Flow {
	double F1;
	double F2;
	double F3;
};

struct —omputationalGrid {
	int numberOfTimeSplits;
	int numberOfXSplits;
};

struct EquationSystem
{
	Variables leftBoundaryCondition;
	Variables rightBoundaryCondition;
	double x0;
	double maxX;
	double minX;
};

struct Response {
	bool success;
	std::string message;
};

class Solver
{
protected:
	Variables bcLeft;
	Variables bcRight;
	double x0;
	int N, K;

	double h, tau, maxT;

	double maxX;
	double minX;

public:
	void virtual solve() {}
	Response virtual init—omputationalGrid(const —omputationalGrid &grid, double maxTime) { return Response { false, "" }; }
	void init(const EquationSystem &eq);
};

