#pragma once
#include "Solver.h"
#include "array.h"

class ASUM : public Solver
{
private:
	double calcSpeedOfSound(Variables var);
	double calcMachNumber(Variables var);
	double mPlus(Variables var);
	double mMinus(Variables var);
	double pPlus(Variables var);
	double pMinus(Variables var);
	double calcBorderMachNumber(Variables &left, Variables &right);
	double calcBorderPressure(Variables &left, Variables &right);
	Flow calcFlow(Variables var, double borderMachNumber, double borderPressure);
	Flow calcFlow(Variables &left, Variables &right);
	void writeVectorOfConservedAndVectorOfFlows(ConservativeVariables& U, Flow& F, Variables var);
	ConservativeVariables calculateVectorOfConservativeValues(ConservativeVariables& lastU, Flow& F_plus, Flow& F_minus);
	void calculateValues(Flow& F, double& ro, double& u, double& p, ConservativeVariables& U);
public:
	ASUM();
	~ASUM();
	void solve();
	Response init—omputationalGrid(const —omputationalGrid &grid, double maxTime);
	void setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, Variables &left, Variables &right);
	void setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p, double x0, Variables &left, Variables &right);

};

