#pragma once
#include "Solver.h"

class Godunov : public Solver
{
	double D_R;
	double D_L;
public:
	Godunov();
	~Godunov();

	Perem solver_riman(Perem U_L, Perem U_R);
	void solve();
	void solver(double *ro, double *u, double *p, int n1);

	double max(double *u, int n1) {
		if (n1 == 0) return 0;//тут исключение
		double m = u[0];
		for (int i = 1;i<n1;i++) {
			if (u[i]>m) m = u[i];
		}
		return m;
	}

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

