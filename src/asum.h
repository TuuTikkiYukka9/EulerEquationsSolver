#pragma once
#include "solver.h"
#include "array.h"

class ASUM : public Solver
{
private:
	double calcSpeedOfSound(const Variables<double>& var);
	double calcMachNumber(const Variables<double>& var);
	double mPlus(const Variables<double>& var);
	double mMinus(const Variables<double>& var);
	double pPlus(const Variables<double>& var);
	double pMinus(const Variables<double>& var);
	double calcBorderMachNumber(const Variables<double> &left, const Variables<double> &right);
	double calcBorderPressure(const Variables<double> &left, const Variables<double> &right);
	Flow calcFlow(const Variables<double>& var, const double& borderMachNumber, const double& borderPressure);
	Flow calcFlow(const Variables<double> &left, const Variables<double> &right);
	//void writeVectorOfConservedAndVectorOfFlows(ConservativeVariables& U, Flow& F, Variables<double> var);
	ConservativeVariables ASUM::getVectorOfConserved(const Variables<double>& var);
	Flow ASUM::getVectorOfFlows(const Variables<double>& var);
	ConservativeVariables calculateVectorOfConservativeValues(const ConservativeVariables& lastU, const Flow& F_plus, const Flow& F_minus);
	Variables<double> ASUM::calculateVariables(const ConservativeVariables& U);
	Flow ASUM::calculateFlow(const Variables<double>& calculatedVariables, const ConservativeVariables& U);
	Variables<Array<double>*> ASUM::getInitialConditions(const int &arrayLength, const double &x0, const Variables<double> &left, const Variables<double> &right);
	Variables<Array<double>*> ASUM::getBoundaryConditions(const int &arrayLength, const Variables<double> &left, const Variables<double> &right);
public:
	ASUM();
	~ASUM();
	Variables<Array<double>*> solve();
	Response init—omputationalGrid(const —omputationalGrid &grid, double maxTime);
	void setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, Variables<double> &left, Variables<double> &right) {};
	void setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p, double x0, Variables<double> &left, Variables<double> &right) {};

};

