#pragma once
#include<string>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "array.h"

template <typename T>
struct Variables {
	T ro;
	T u;
	T p;
};

template <typename T>
struct Discontinuity {
	T left;
	T right;
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
	Variables<double> leftBoundaryCondition;
	Variables<double> rightBoundaryCondition;
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
	Variables<double> bcLeft;
	Variables<double> bcRight;
	double x0;
	int N, K;

	double h, tau, maxT;

	double maxX;
	double minX;

public:
	Variables<Array<double>*> virtual solve() { return { nullptr, nullptr, nullptr }; }
	Response virtual init—omputationalGrid(const —omputationalGrid &grid, double maxTime) { return Response { false, "" }; }
	void init(const EquationSystem &eq);
};

