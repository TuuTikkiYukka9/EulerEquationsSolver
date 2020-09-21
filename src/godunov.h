#pragma once
#include "Solver.h"
#include "array.h"

class Godunov : public Solver
{
	double rightD;
	double leftD;
public:
	Godunov();
	~Godunov();

	Variables solver_riman(Variables U_L, Variables U_R);
	void solve();
	Response init—omputationalGrid(const —omputationalGrid &grid, double maxTime);
	void solver(Array<double> &ro, Array<double> &u, Array<double> &p, int n1);
	void setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, Variables &left, Variables &right);
	void setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p, double x0, Variables &left, Variables &right);

	double max(double a, double b) {
		if (a>b) return a;
		else return b;
	}
	double min(double a, double b) {
		if (b == 0) {
			return -1;
		}
		else if (b < a) {
			return b;
		}
		else {
			return a;
		}
	}
};

