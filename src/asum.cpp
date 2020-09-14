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

	double a_L, a_R, M_L, M_R;

	Matrix<consv> U(nk, n1);
	Matrix<potoc> F(nk, n1);
	Matrix<consv> U_drob(nk, n1);
	Matrix<potoc> F_drob(nk, n1);

	ro.SetAllElements(0);
	p.SetAllElements(0);
	u.SetAllElements(0);

	double m_plus, m_minus;
	double p_plus, p_minus;

	double p_L = bcLeft.p, ro_L = bcLeft.ro, u_L = bcLeft.u;//вот это надо считывать из файла
	double p_R = bcRight.p, ro_R = bcRight.ro, u_R = bcRight.u;
	double x_0 = x0;
	double currentCoordinate;
	bool isLeft;

	//начальное условие
	for (int i = 0; i < p.height(); i++) {
		currentCoordinate = minX + (h*i);
		if (currentCoordinate == x_0) {
			a_L = sqrt(1.4*p_L / ro_L);
			M_L = u_L / a_L;
			m_plus = (fabs(M_L) <= 1) ? (0.25*pow((M_L + 1), 2)) : (0.5*(M_L + fabs(M_L)));
			a_R = sqrt(1.4*p_R / ro_R);
			M_R = u_R / a_R;
			m_minus = (fabs(M_R) <= 1) ? (-0.25*pow((M_R - 1), 2)) : (0.5*(M_R - fabs(M_R)));
		}
		
		isLeft = currentCoordinate < x_0 || (currentCoordinate == x_0 && (m_plus + m_minus) >= 0);
		
		p[0][i] = isLeft ? p_L : p_R;
		ro[0][i] = isLeft ? ro_L : ro_R;
		u[0][i] = isLeft ? u_L : u_R;
	}

	//граничные условия
	for (int k = 0; k < p.height(); k++) {
		ro[k].setFirst(ro_L);
		u[k].setFirst(u_L);
		p[k].setFirst(p_L);
		ro[k].setLast(ro_R);
		u[k].setLast(u_R);
		p[k].setLast(p_R);
	}

	for (int k = 0; k < p.height(); k++) {
		for (int i = 0;i<p.width();i++) {
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
				struct potoc F_plus, F_minus, calculetedFlow;

				for (int iter = 0; iter<2;iter++) {
					int j = i + iter;
					p_L = p[k][j - 1];//ro
					ro_L = ro[k][j - 1];// ro8u
					u_L = u[k][j - 1];//E

					p_R = p[k][j];//ro
					ro_R = ro[k][j];// ro8u
					u_R = u[k][j];//E

					a_L = sqrt(1.4*p_L / ro_L);
					M_L = u_L / a_L;
					m_plus = fabs(M_L) <= 1 ? 0.25*pow((M_L + 1), 2) : 0.5*(M_L + fabs(M_L));
					p_plus = fabs(M_L) <= 1 ?  0.5*p_L*(1 + M_L) : 0.5*p_L*((M_L + fabs(M_L)) / M_L);
					

					a_R = sqrt(1.4*p_R / ro_R);
					M_R = u_R / a_R;
					m_minus = fabs(M_R) <= 1 ? -0.25*pow((M_R - 1), 2) : 0.5*(M_R - fabs(M_R));
					p_minus = fabs(M_R) <= 1 ? 0.5*p_R*(1 - M_R) : 0.5*p_R*((M_R - fabs(M_R)) / M_R);;
					

					if ((m_plus + m_minus) >= 0) {
						calculetedFlow.F1 = (m_plus + m_minus)*ro_L*a_L + 0;
						calculetedFlow.F2 = (m_plus + m_minus)*ro_L*u_L*a_L + (p_plus + p_minus);
						calculetedFlow.F3 = (m_plus + m_minus)*(p_L + (p_L / 0.4 + (ro_L*pow(u_L, 2) / 2)))*a_L + 0;
					}
					else {
						calculetedFlow.F1 = (m_plus + m_minus)*ro_R*a_R + 0;
						calculetedFlow.F2 = (m_plus + m_minus)*ro_R*u_R*a_R + (p_plus + p_minus);
						calculetedFlow.F3 = (m_plus + m_minus)*(p_R + (p_R / 0.4 + (ro_R*pow(u_R, 2) / 2)))*a_R + 0;
					}

					if (iter == 0) {
						F_minus = calculetedFlow;
					}
					else {
						F_plus = calculetedFlow;
					}
				}


				U[k + 1][i].U1 = U[k][i].U1 - (tau / h)*(F_plus.F1 - F_minus.F1);
				U[k + 1][i].U2 = U[k][i].U2 - (tau / h)*(F_plus.F2 - F_minus.F2);
				U[k + 1][i].U3 = U[k][i].U3 - (tau / h)*(F_plus.F3 - F_minus.F3);
				//_____________________________end AUSM__________________________________//

		}

		for (int i = 1;i< p.width() - 1;i++) {
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


ASUM::~ASUM()
{
}
