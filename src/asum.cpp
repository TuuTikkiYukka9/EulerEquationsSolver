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

	const int nmbOfTimeSplits = K + 1;
	const int nmbOfXSplits = N + 1;

	Matrix<Variables<double>> res(nmbOfTimeSplits, nmbOfXSplits, {0, 0, 0});

	Matrix<ConservativeVariables> U(nmbOfTimeSplits, nmbOfXSplits, {0, 0, 0});
	Matrix<Flow> F(nmbOfTimeSplits, nmbOfXSplits, {0, 0, 0});
	Flow F_plus, F_minus;
	
	res[0] = getInitialConditions(res.width(), x0, bcLeft, bcRight);

	//граничные условия
	for (int k = 1; k < res.height(); k++) {
		res[k].setFirst(bcLeft);
		res[k].setLast(bcRight);
	}
	//Предствим уравнение как вектор консервативных величин и вектор потоков
	for (int k = 0; k < res.height(); k++) {
		for (int i = 0; i < res.width(); i++) {
			U[k][i] = getVectorOfConserved(res[k][i]); 
			F[k][i] = getVectorOfFlows(res[k][i]);
		}
	}

	for (int k = 0; k < res.height() - 1; k++) {
		for (int i = 1; i < res.width() - 1; i++) {
			//_________________________AUSM_______________________________//	

			for (int iter = 0; iter < 2; iter++) {
				int j = i + iter;

				if (iter == 0) {
					F_minus = calcFlow(res[k][j - 1], res[k][j]);
				}
				else {
					F_plus = calcFlow(res[k][j - 1], res[k][j]);
				}
			}

			U[k + 1][i] = calculateVectorOfConservativeValues(U[k][i], F_plus, F_minus);
			//_____________________________end AUSM__________________________________//

		}

		for (int i = 1; i < res.width() - 1; i++) {
			res[k + 1][i] = calculateVariables(U[k + 1][i]);
			F[k + 1][i] = calculateFlow(res[k + 1][i], U[k + 1][i]);
		}
	}

	Array<double> *ResultRo = new Array<double>(nmbOfXSplits);
	Array<double> *ResultU = new Array<double>(nmbOfXSplits);
	Array<double> *ResultP = new Array<double>(nmbOfXSplits);
	for (int i = 0;i< ResultP->length(); i++){
		(*ResultRo)[i] = res[nmbOfTimeSplits - 1][i].ro;
		(*ResultU)[i] = res[nmbOfTimeSplits - 1][i].u;
		(*ResultP)[i] = res[nmbOfTimeSplits - 1][i].p;
	}
	return { ResultRo, ResultU, ResultP };
}

Array<Variables<double>> ASUM::getInitialConditions(const int &arrayLength, const double &x0,
	const Variables<double> &left, const Variables<double> &right) {

	Array<Variables<double>> result(arrayLength, { 0, 0, 0 });

	double currentCoordinate;
	bool isLeft;

	double borderMachNumber;
	
	for (int i = 0; i < arrayLength; i++) {
		currentCoordinate = minX + (h * i);
		if (currentCoordinate == x0) {
			borderMachNumber = calcBorderMachNumber(left, right);
		}

		isLeft = currentCoordinate < x0 || (currentCoordinate == x0 && borderMachNumber >= 0);
		result[i] = isLeft ? left : right;

	}
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
		   (0.25 * pow((calcMachNumber(var) + 1), 2)) : 
		   (0.5 * (calcMachNumber(var) + fabs(calcMachNumber(var))));
}

double ASUM::mMinus(const Variables<double>& var) {
	return (fabs(calcMachNumber(var)) <= 1) ? 
		   (-0.25 * pow((calcMachNumber(var) - 1), 2)) : 
		   (0.5 * (calcMachNumber(var) - fabs(calcMachNumber(var))));
}

double ASUM::pPlus(const Variables<double>& var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5 * var.p * (1 + calcMachNumber(var)) :
		   0.5 * var.p * ((calcMachNumber(var) + fabs(calcMachNumber(var))) / calcMachNumber(var));
}

double ASUM::pMinus(const Variables<double>& var) {
	return fabs(calcMachNumber(var)) <= 1 ? 
		   0.5 * var.p * (1 - calcMachNumber(var)) :
		   0.5 * var.p * ((calcMachNumber(var) - fabs(calcMachNumber(var))) / calcMachNumber(var));
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