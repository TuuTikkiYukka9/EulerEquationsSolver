#include "ASUM.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "solutionwriter.h"
#include "matrix.h"

using namespace std;
ASUM::ASUM()
{
	N = 0;
	K = 0;
	maxT = 0;
	h = 0;
	tau = 0;
}

Response ASUM::initСomputationalGrid(const СomputationalGrid &grid, double maxTime) {
	const double dh = (maxX - minX) / grid.numberOfXSplits;
	const double dT = maxTime / grid.numberOfTimeSplits;
	const double courantNumber = dT / (2 * dh);
	if (courantNumber > 1.0) return Response{ false, "Courant number should not be more than 1!" };
	N = grid.numberOfXSplits;
	K = grid.numberOfTimeSplits;
	maxT = maxTime;
	h = dh;
	tau = dT;
	return Response{ true, "" };
}


void ASUM::solve() {

	if (maxT == 0) return;

	const int nk = (K + 1);
	const int n1 = N + 1;

	Matrix<double> ro(nk, n1);
	Matrix<double> u(nk, n1);
	Matrix<double> p(nk, n1);

	Matrix<ConservativeVariables> U(nk, n1);
	Matrix<Flow> F(nk, n1);
	Matrix<ConservativeVariables> U_drob(nk, n1);
	Matrix<Flow> F_drob(nk, n1);

	ro.SetAllElements(0);
	p.SetAllElements(0);
	u.SetAllElements(0);

	Variables left, right;
	Flow F_plus, F_minus;

	
	//начальное условие
	setInitialConditions(ro[0], u[0], p[0], x0, bcLeft, bcRight);

	//граничные условия
	for (int k = 0; k < p.height(); k++) {
		setBoundaryConditions(ro[k], u[k], p[k], bcLeft, bcRight);
	}

	//Предствим уравнение как вектор консервативных величин и вектор потоков
	for (int k = 0; k < p.height(); k++) {
		for (int i = 0; i < p.width(); i++) {
			writeVectorOfConservedAndVectorOfFlows(U[k][i], F[k][i], { ro[k][i], u[k][i], p[k][i] });
		}
	}

	for (int k = 0; k < p.height() - 1; k++) {
		for (int i = 1; i < p.width() - 1; i++) {
			//_________________________AUSM_______________________________//	

			for (int iter = 0; iter < 2; iter++) {
				int j = i + iter;
					
				left = { ro[k][j - 1], u[k][j - 1], p[k][j - 1] };
				right = { ro[k][j], u[k][j], p[k][j] };

				if (iter == 0) {
					F_minus = calcFlow(left, right);
				}
				else {
					F_plus = calcFlow(left, right);
				}
			}

			U[k + 1][i] = calculateVectorOfConservativeValues(U[k][i], F_plus, F_minus);
			//_____________________________end AUSM__________________________________//

		}

		for (int i = 1; i < p.width() - 1; i++) {
			calculateValues(F[k + 1][i], ro[k + 1][i], u[k + 1][i], p[k + 1][i], U[k + 1][i]);
		}
	}

	//----------------------Testing----------------------------------------//
	(new SolutionWriter())->write("AUSM", ro[nk - 1], u[nk - 1], p[nk - 1]);
}


void ASUM::calculateValues(Flow& outputF, double& outputRo, double& outputU, double& outputP, ConservativeVariables& U) {
	outputRo = U.U1;
	outputU = U.U2 / U.U1;
	outputP = (1.4 - 1.0)*(U.U3 - (pow(U.U2, 2)) / (2.0*U.U1));

	outputF.F1 = U.U2;
	outputF.F2 = outputP + outputRo * pow(outputU, 2);
	outputF.F3 = outputU * (outputP + U.U3);
}


ConservativeVariables ASUM::calculateVectorOfConservativeValues(ConservativeVariables& lastU, Flow& F_plus, Flow& F_minus) {
	ConservativeVariables U;
	U.U1 = lastU.U1 - (tau / h)*(F_plus.F1 - F_minus.F1);
	U.U2 = lastU.U2 - (tau / h)*(F_plus.F2 - F_minus.F2);
	U.U3 = lastU.U3 - (tau / h)*(F_plus.F3 - F_minus.F3);
	return U;
}


void ASUM::writeVectorOfConservedAndVectorOfFlows(ConservativeVariables& outputU, Flow& outputF, Variables var) {
	outputU.U1 = var.ro;
	outputU.U2 = var.ro * var.u;
	outputU.U3 = (var.p / (1.4 - 1.0)) + 0.5*var.ro * pow(var.u, 2);

	outputF.F1 = var.ro * var.u;
	outputF.F2 = var.p + var.ro * pow(var.u, 2);
	outputF.F3 = var.u * (var.p + outputU.U3);
}

double ASUM::calcSpeedOfSound(Variables var){
	double gamma = 1.4;
	return sqrt(gamma * var.p / var.ro);
}

double ASUM::calcMachNumber(Variables var) {
	return var.u / calcSpeedOfSound(var);
}

double ASUM::mPlus(Variables var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		   (0.25*pow((calcMachNumber(var) + 1), 2)) : 
		   (0.5*(calcMachNumber(var) + fabs(calcMachNumber(var))));
}
double ASUM::mMinus(Variables var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		   (-0.25*pow((calcMachNumber(var) - 1), 2)) : 
		   (0.5*(calcMachNumber(var) - fabs(calcMachNumber(var))));
}
double ASUM::pPlus(Variables var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 + calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) + fabs(calcMachNumber(var))) / calcMachNumber(var));
}
double ASUM::pMinus(Variables var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 - calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) - fabs(calcMachNumber(var))) / calcMachNumber(var));
}


double ASUM::calcBorderMachNumber(Variables &left, Variables &right){
	return mPlus(left) + mMinus(right);
}
double ASUM::calcBorderPressure(Variables &left, Variables &right) {
	return pPlus(left) + pMinus(right);
}
Flow ASUM::calcFlow(Variables var, double borderMachNumber, double borderPressure)
{
	return {
		borderMachNumber * var.ro * calcSpeedOfSound(var) + 0,
		borderMachNumber * var.ro * var.u * calcSpeedOfSound(var) + borderPressure,
		borderMachNumber * (var.p + (var.p / 0.4 + (var.ro * pow(var.u, 2) / 2))) * calcSpeedOfSound(var) + 0
	};
}
Flow ASUM::calcFlow(Variables &left, Variables &right)
{
	double borderMachNumber = calcBorderMachNumber(left, right);
	double borderPressure = calcBorderPressure(left, right);
	return (borderMachNumber >= 0) ?
		   calcFlow(left, borderMachNumber, borderPressure) :
		   calcFlow(right, borderMachNumber, borderPressure);
}

void ASUM::setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p,
	Variables &left, Variables &right) {
	ro.setFirst(left.ro);
	u.setFirst(left.u);
	p.setFirst(left.p);

	ro.setLast(right.ro);
	u.setLast(right.u);
	p.setLast(right.p);
}

void ASUM::setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p,
	double x0, Variables &left, Variables &right) {
	
	double currentCoordinate;
	bool isLeft;

	double borderMachNumber;

	 for (int i = 0; i < p.length(); i++) {
		currentCoordinate = minX + (h*i);
		if (currentCoordinate == x0) {
			borderMachNumber = calcBorderMachNumber(left, right);
		}

		isLeft = currentCoordinate < x0 || (currentCoordinate == x0 && borderMachNumber >= 0);

		p[i] = isLeft ? left.p : right.p;
		ro[i] = isLeft ? left.ro : right.ro;
		u[i] = isLeft ? left.u: right.u;
	}
}

ASUM::~ASUM()
{
}
