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
	minX = 0.0;
	maxX = 1.0;

	double Kur = 2;
	while (Kur>1.0) {
		cout << "Enter the number of partitions for x:";
		cin >> N;
		cout << "Enter the number of partitions for t:";
		cin >> K;
		cout << "Enter max T:";
		cin >> maxT;
		h = (maxX - minX) / N;
		tau = (maxT) / K;
		Kur = (tau) / (2 * h);
		cout << "chislo_Kuranta:" << Kur << "\n";
	}

}


void ASUM::solve() {

	const int nk = (K + 1);
	const int n1 = N + 1;

	Matrix<double> ro(nk, n1);
	Matrix<double> u(nk, n1);
	Matrix<double> p(nk, n1);

	Matrix<consv> U(nk, n1);
	Matrix<potoc> F(nk, n1);
	Matrix<consv> U_drob(nk, n1);
	Matrix<potoc> F_drob(nk, n1);

	ro.SetAllElements(0);
	p.SetAllElements(0);
	u.SetAllElements(0);

	Variables left, right;
	potoc F_plus, F_minus;

	
	//начальное условие
	setInitialConditions(ro[0], u[0], p[0], x0, bcLeft, bcRight);

	//граничные условия
	for (int k = 0; k < p.height(); k++) {
		setBoundaryConditions(ro[k], u[k], p[k], bcLeft, bcRight);
	}

	for (int k = 0; k < p.height(); k++) {
		for (int i = 0; i < p.width(); i++) {
			U[k][i].U1 = ro[k][i];
			U[k][i].U2 = ro[k][i] * u[k][i];
			U[k][i].U3 = (p[k][i] / (1.4 - 1.0)) + 0.5*ro[k][i] * pow(u[k][i], 2);

			F[k][i].F1 = ro[k][i] * u[k][i];
			F[k][i].F2 = p[k][i] + ro[k][i] * pow(u[k][i], 2);
			F[k][i].F3 = u[k][i] * (p[k][i] + U[k][i].U3);
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

			U[k + 1][i].U1 = U[k][i].U1 - (tau / h)*(F_plus.F1 - F_minus.F1);
			U[k + 1][i].U2 = U[k][i].U2 - (tau / h)*(F_plus.F2 - F_minus.F2);
			U[k + 1][i].U3 = U[k][i].U3 - (tau / h)*(F_plus.F3 - F_minus.F3);
			//_____________________________end AUSM__________________________________//

		}

		for (int i = 1; i < p.width() - 1; i++) {
			ro[k + 1][i] = U[k + 1][i].U1;
			u[k + 1][i] = U[k + 1][i].U2 / U[k + 1][i].U1;
			p[k + 1][i] = (1.4 - 1.0)*(U[k + 1][i].U3 - (pow(U[k + 1][i].U2, 2)) / (2.0*U[k + 1][i].U1));

			F[k + 1][i].F1 = U[k + 1][i].U2;
			F[k + 1][i].F2 = p[k + 1][i] + ro[k + 1][i] * pow(u[k + 1][i], 2);
			F[k + 1][i].F3 = u[k + 1][i] * (p[k + 1][i] + U[k + 1][i].U3);
		}
	}

	//----------------------Testing----------------------------------------//
	(new SolutionWriter())->write("AUSM", ro[nk - 1], u[nk - 1], p[nk - 1]);
}

double ASUM::calcSpeedOfSound(Variables var){
	double gamma = 1.4;
	return sqrt(gamma * var.p / var.ro);
}

double ASUM::calcMachNumber(Variables var) {
	return var.u / calcSpeedOfSound(var);
}

double ASUM::m_plus(Variables var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		     (0.25*pow((calcMachNumber(var) + 1), 2)) : 
		   (0.5*(calcMachNumber(var) + fabs(calcMachNumber(var))));
}
double ASUM::m_minus(Variables var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		    (-0.25*pow((calcMachNumber(var) - 1), 2)) : 
		  (0.5*(calcMachNumber(var) - fabs(calcMachNumber(var))));
}
double ASUM::p_plus(Variables var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 + calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) + fabs(calcMachNumber(var))) / calcMachNumber(var));
}
double ASUM::p_minus(Variables var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 - calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) - fabs(calcMachNumber(var))) / calcMachNumber(var));
}


double ASUM::calcBorderMachNumber(Variables &left, Variables &right){
	return m_plus(left) + m_minus(right);
}
double ASUM::calcBorderPressure(Variables &left, Variables &right) {
	return p_plus(left) + p_minus(right);
}
potoc ASUM::calcFlow(Variables var, double borderMachNumber, double borderPressure)
{
	potoc flow;
	flow.F1 = borderMachNumber * var.ro * calcSpeedOfSound(var) + 0;
	flow.F2 = borderMachNumber * var.ro * var.u * calcSpeedOfSound(var) + borderPressure;
	flow.F3 = borderMachNumber * (var.p + (var.p / 0.4 + (var.ro * pow(var.u, 2) / 2))) * calcSpeedOfSound(var) + 0;
	return flow;
}
potoc ASUM::calcFlow(Variables &left, Variables &right)
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
