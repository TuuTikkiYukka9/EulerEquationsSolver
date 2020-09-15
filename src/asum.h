#pragma once
#include "Solver.h"
#include "array.h"

class ASUM : public Solver
{
private:
	double calcSpeedOfSound(Variables var);
	double calcMachNumber(Variables var);
	double m_plus(Variables var);
	double m_minus(Variables var);
	double p_plus(Variables var);
	double p_minus(Variables var);
	double calcBorderMachNumber(Variables &left, Variables &right);
	double calcBorderPressure(Variables &left, Variables &right);
	potoc calcFlow(Variables var, double borderMachNumber, double borderPressure);
	potoc calcFlow(Variables &left, Variables &right);
public:
	ASUM();
	~ASUM();
	void solve();
	void setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, Variables &left, Variables &right);
	void setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p, double x0, Variables &left, Variables &right);

};

