#pragma once
#include "solver.h"
#include "array.h"

class Godunov : public Solver
{
	double rightD;
	double leftD;
public:
	Godunov();
	~Godunov();

	Variables<double> solver_riman(Variables<double> U_L, Variables<double> U_R);
	Variables<Array<double>*> solve();
	Response init—omputationalGrid(const —omputationalGrid &grid, double maxTime);
	void solver(Array<double> &ro, Array<double> &u, Array<double> &p, int n1);
	void setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, Variables<double> &left, Variables<double> &right);
	void setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p, double x0, Variables<double> &left, Variables<double> &right);

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

