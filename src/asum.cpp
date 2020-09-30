#include "asum.h"
#include "matrix.h"


ASUM::ASUM() {
	N = 0;
	K = 0;
	maxT = 0;
	h = 0;
	tau = 0;
}

Response ASUM::initСomputationalGrid(const СomputationalGrid &grid, const double maxTime) {
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


Variables<Array<double>*> ASUM::solve() {

	if (maxT == 0) return { nullptr, nullptr, nullptr };

	const int nk = (K + 1);
	const int n1 = N + 1;

	Matrix<double> ro(nk, n1);
	Matrix<double> u(nk, n1);
	Matrix<double> p(nk, n1);
	Array<double> *ResultRo = new Array<double>(n1);
	Array<double> *ResultU = new Array<double>(n1);
	Array<double> *ResultP = new Array<double>(n1);

	Matrix<ConservativeVariables> U(nk, n1);
	Matrix<Flow> F(nk, n1);
	Matrix<ConservativeVariables> U_drob(nk, n1);
	Matrix<Flow> F_drob(nk, n1);

	ro.SetAllElements(0);
	p.SetAllElements(0);
	u.SetAllElements(0);

	Variables<double> left, right;
	Flow F_plus, F_minus;
	Variables<double> calculatedVariables;
	//начальное условие

	Variables<Array<double>*> initValues = getInitialConditions(p[0].length(), x0, bcLeft, bcRight);
	ro[0] = (*initValues.ro);
	u[0] = (*initValues.u);
	p[0] = (*initValues.p);
	delete initValues.ro;
	delete initValues.u;
	delete initValues.p;

	//граничные условия
	for (int k = 1; k < p.height(); k++) {
		ro[k].setFirst(bcLeft.ro);
		u[k].setFirst(bcLeft.u);
		p[k].setFirst(bcLeft.p);

		ro[k].setLast(bcRight.ro);
		u[k].setLast(bcRight.u);
		p[k].setLast(bcRight.p);
	}
	//Предствим уравнение как вектор консервативных величин и вектор потоков
	for (int k = 0; k < p.height(); k++) {
		for (int i = 0; i < p.width(); i++) {
			U[k][i] = getVectorOfConserved({ ro[k][i], u[k][i], p[k][i] });
			F[k][i] = getVectorOfFlows({ ro[k][i], u[k][i], p[k][i] });
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
			calculatedVariables = calculateVariables(U[k + 1][i]);
			F[k + 1][i] = calculateFlow(calculatedVariables, U[k + 1][i]);
			ro[k + 1][i] = calculatedVariables.ro;
			u[k + 1][i] = calculatedVariables.u;
			p[k + 1][i] = calculatedVariables.p;
		}
	}

	(*ResultRo) = ro[nk - 1];
	(*ResultU) = u[nk - 1];
	(*ResultP) = p[nk - 1];
	p[nk - 1].print();
	return { ResultRo, ResultU, ResultP };
}

Variables<Array<double>*> ASUM::getBoundaryConditions(const int &arrayLength, const Variables<double> &left, const Variables<double> &right) {
	Variables<Array<double>*> result;
	result.p = new Array<double>(arrayLength, 0);
	result.ro = new Array<double>(arrayLength, 0);
	result.u = new Array<double>(arrayLength, 0);

	result.ro->setFirst(left.ro);
	result.u->setFirst(left.u);
	result.p->setFirst(left.p);

	result.ro->setLast(right.ro);
	result.u->setLast(right.u);
	result.p->setLast(right.p);
	return result;
}

Variables<Array<double>*> ASUM::getInitialConditions(const int &arrayLength, const double &x0,
	const Variables<double> &left, const Variables<double> &right) {

	Variables<Array<double>*> result;
	result.p = new Array<double>(arrayLength, 0);
	result.ro = new Array<double>(arrayLength, 0);
	result.u = new Array<double>(arrayLength, 0);

	double currentCoordinate;
	bool isLeft;

	double borderMachNumber;
	
	for (int i = 0; i < arrayLength; i++) {
		currentCoordinate = minX + (h*i);
		if (currentCoordinate == x0) {
			borderMachNumber = calcBorderMachNumber(left, right);
		}

		isLeft = currentCoordinate < x0 || (currentCoordinate == x0 && borderMachNumber >= 0);
		(*result.p)[i] = isLeft ? left.p : right.p;
		(*result.ro)[i] = isLeft ? left.ro : right.ro;
		(*result.u)[i] = isLeft ? left.u : right.u;
	}
	result.p->print();
	return result;
}

ConservativeVariables ASUM::getVectorOfConserved(const Variables<double>& var) {
	return {
		var.ro,
		var.ro * var.u,
	    (var.p / (1.4 - 1.0)) + 0.5 * var.ro * pow(var.u, 2)
	};
}

Flow ASUM::getVectorOfFlows(const Variables<double>& var) {
	return {
		var.ro * var.u,
		var.p + var.ro * pow(var.u, 2),
		var.u * (var.p + (var.p / (1.4 - 1.0)) + 0.5*var.ro * pow(var.u, 2))
	};
}

ConservativeVariables ASUM::calculateVectorOfConservativeValues(const ConservativeVariables& lastU, 
	                                                            const Flow& F_plus, const Flow& F_minus) {
	return { lastU.U1 - (tau / h)*(F_plus.F1 - F_minus.F1),
			 lastU.U2 - (tau / h)*(F_plus.F2 - F_minus.F2),
			 lastU.U3 - (tau / h)*(F_plus.F3 - F_minus.F3) } ;
}

double ASUM::calcSpeedOfSound(const Variables<double>& var){
	const double gamma = 1.4;
	return sqrt(gamma * var.p / var.ro);
}

double ASUM::calcMachNumber(const Variables<double>& var) {
	return var.u / calcSpeedOfSound(var);
}

double ASUM::mPlus(const Variables<double>& var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		   (0.25*pow((calcMachNumber(var) + 1), 2)) : 
		   (0.5*(calcMachNumber(var) + fabs(calcMachNumber(var))));
}

double ASUM::mMinus(const Variables<double>& var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		   (-0.25*pow((calcMachNumber(var) - 1), 2)) : 
		   (0.5*(calcMachNumber(var) - fabs(calcMachNumber(var))));
}

double ASUM::pPlus(const Variables<double>& var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 + calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) + fabs(calcMachNumber(var))) / calcMachNumber(var));
}

double ASUM::pMinus(const Variables<double>& var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5*var.p*(1 - calcMachNumber(var)) :
		   0.5*var.p*((calcMachNumber(var) - fabs(calcMachNumber(var))) / calcMachNumber(var));
}

double ASUM::calcBorderMachNumber(const Variables<double> &left, const Variables<double> &right){
	return mPlus(left) + mMinus(right);
}

double ASUM::calcBorderPressure(const Variables<double> &left, const Variables<double> &right) {
	return pPlus(left) + pMinus(right);
}

Flow ASUM::calcFlow(const Variables<double>& var, const double &borderMachNumber, const double &borderPressure)
{
	return {
		borderMachNumber * var.ro * calcSpeedOfSound(var) + 0,
		borderMachNumber * var.ro * var.u * calcSpeedOfSound(var) + borderPressure,
		borderMachNumber * (var.p + (var.p / 0.4 + (var.ro * pow(var.u, 2) / 2))) * calcSpeedOfSound(var) + 0
	};
}

Flow ASUM::calcFlow(const Variables<double> &left, const Variables<double> &right)
{
	const double borderMachNumber = calcBorderMachNumber(left, right);
	const double borderPressure = calcBorderPressure(left, right);
	return (borderMachNumber >= 0) ?
		   calcFlow(left, borderMachNumber, borderPressure) :
		   calcFlow(right, borderMachNumber, borderPressure);
}

Variables<double> ASUM::calculateVariables(const ConservativeVariables& U) {
	return {
	    U.U1,
		U.U2 / U.U1,
		(1.4 - 1.0) * (U.U3 - (pow(U.U2, 2)) / (2.0 * U.U1))
	};
}

Flow ASUM::calculateFlow(const Variables<double>& calculatedVariables, const ConservativeVariables& U) {
	return {
		U.U2,
		calculatedVariables.p + calculatedVariables.ro * pow(calculatedVariables.u, 2),
		calculatedVariables.u * (calculatedVariables.p + U.U3)
	};
}

ASUM::~ASUM() {
}