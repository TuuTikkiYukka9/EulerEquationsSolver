#include "Godunov.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "solutionwriter.h"
#include "array.h"


using namespace std;

//-------------------------------------------------------------------------
/*Не должен конструктор ничего печатать, начальные и краевые услови§ то же надо передавать 4 вектора по 3 элемента, и еще положение разрыва х0, до звуковое течение*/
Godunov::Godunov(void)
{
	minX = 0.0;
	maxX = 1.0;

	cout << "Enter the number of partitions for x:";
	cin >> N;
	cout << "Enter max T:";
	cin >> maxT;

	h = (maxX - minX) / N;
	tau = 0.411 * h;
	K = (int)(maxT / tau);

}

Variables Godunov::solver_riman(Variables U_L, Variables U_R) {
	double gamma = 1.4;
	struct Variables U;
	double cL = 1;
	double cR = 1;
	double Ro;

	double aL = pow(gamma*(U_L.p / U_L.ro), 0.5);
	double aR = pow(gamma*(U_R.p / U_R.ro), 0.5);
	double P0 = (aR*U_R.ro*U_L.p + aL*U_L.ro*U_R.p + aR*U_R.ro*aL*U_L.ro*(U_L.u - U_R.u)) / (aR*U_R.ro + aL*U_L.ro);
	double U0 = (aR*U_R.ro*U_R.u + aL*U_L.ro*U_L.u + U_L.p - U_R.p) / (aR*U_R.ro + aL*U_L.ro);

	double P_buf = 3;
	double U_buf = 3000;

	double eps = 0.000001;//eps=0.00000000001;
	while ((fabs(P0 - P_buf)) >= eps || (fabs(U0 - U_buf)) >= eps) {
		P_buf = P0;
		U_buf = U0;

		if (P0<U_L.p) {
			cL = aL * U_L.ro*(gamma - 1)*(1 - (P0 / U_L.p)) / (2 * gamma*(1 - pow(P0 / U_L.p, (gamma - 1) / (2 * gamma))));
			//cout<<"veer voln razr vlevo \n";
			leftD = U_L.u - cL / U_L.ro;
		}
		else if (P0>U_L.p) {
			cL = sqrt(U_L.ro*((gamma + 1)*P0 + (gamma - 1)*U_L.p) / 2.0);
			//cout<<"Udarnaia volna vlevo\n";
			leftD = U_L.u - aL;
		}


		if (P0<U_R.p) {
			cR = aR * U_R.ro*(gamma - 1)*(1 - P0 / U_R.p) / (2 * gamma*(1 - pow(P0 / U_R.p, (gamma - 1) / (2 * gamma))));
			//cout<<"veer voln razr vpravo\n";
			rightD = U_R.u + cR / U_R.ro;
		}
		else if (P0>U_R.p) {
			cR = sqrt(U_R.ro*((gamma + 1)*P0 + (gamma - 1)*U_R.p) / 2.0);
			//cout<<"Udarnaia volna vpravo \n";
			rightD = U_R.u + aR;
		}

		P0 = (cR * U_L.p + cL *U_R.p + cR * cL * (U_L.u - U_R.u)) / (cR + cL);
		U0 = (cL * U_L.u + cR * U_R.u + U_L.p - U_R.p) / (cR + cL);
	}

	//Ччитаем ро
	Ro = (U_L.ro + U_R.ro) / 2;
	if (U0<0) {
		if (P0>U_R.p) {//ударна§ волна
			Ro = U_R.ro*((gamma + 1)*P0 + (gamma - 1)*U_R.p) / ((gamma - 1)*P0 + (gamma + 1)*U_R.p);
		}
		if (P0<U_R.p) {//веер волн
			double a_new_R = aR - ((gamma - 1) / 2)*(U_R.u - U0);
			Ro = gamma*P0 / pow(a_new_R, 2);
		}
	}
	else if (U0>0) {
		if (P0>U_L.p) {//ударна§ волна
			Ro = U_L.ro*((gamma + 1)*P0 + (gamma - 1)*U_L.p) / ((gamma - 1)*P0 + (gamma + 1)*U_L.p);
		}
		if (P0<U_L.p) {//веер волн
			double a_new_L = aL + ((gamma - 1) / 2)*(U_L.u - U0);
			Ro = gamma*P0 / pow(a_new_L, 2);
		}
	}

	if (P0>0) { //&& P0<=1){//провер§ем
		U.ro = Ro;
		U.u = U0;
		U.p = P0;
	}
	else {//если конфигураци§ разрывов не вы§вленна
		U.ro = (U_L.ro + U_R.ro) / 2;
		U.u = (U_L.u + U_R.u) / 2;
		U.p = (U_L.p + U_R.p) / 2;
		leftD = U.u;//04.01
		rightD = U.u;
	}

	return U;
}

void Godunov::solver(Array<double> &ro, Array<double> &u, Array<double> &p, int n1) {
	double gamma = 1.4;
	Array<ConservativeVariables> U(n1);
	Array<ConservativeVariables> U_new(n1); 
	Array<Flow> F(n1);
	
	Array<ConservativeVariables> U_drob(n1);
	Array<Flow> F_drob(n1);

	for (int i = 0; i < n1; i++) {
		U[i].U1 = ro[i];
		U[i].U2 = ro[i] * u[i];
		double e = ro[i] * ((1.0 / (1.4 - 1.0))*(p[i] / ro[i]));
		U[i].U3 = (p[i] / (gamma - 1)) + 0.5*ro[i] * pow(u[i], 2);//

	}

	for (int i = 1; i<n1 - 1; i++) {
			struct Flow F_plus, F_minus;
			double E;
			if (p[i - 1] == p[i] && ro[i - 1] == ro[i] && u[i - 1] == u[i]) {

				E = (p[i] / (gamma - 1)) + 0.5*ro[i] * pow(u[i], 2);

				F_minus.F1 = ro[i] * u[i];
				F_minus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_minus.F3 = u[i] * (p[i] + E);
			}
			else {
				Variables U_L = { ro[i - 1], u[i - 1], p[i - 1] };
				Variables U_R = { ro[i], u[i], p[i] };
				Variables U_res;
				U_res = solver_riman(U_L, U_R);

				E = (U_res.p / (gamma - 1)) + 0.5*U_res.ro*pow(U_res.u, 2);

				F_minus.F1 = U_res.ro*U_res.u;
				F_minus.F2 = U_res.p + U_res.ro*pow(U_res.u, 2);
				F_minus.F3 = U_res.u*(U_res.p + E);


			}

			if (p[i] == p[i + 1] && ro[i] == ro[i + 1] && u[i] == u[i + 1]) {
				E = (p[i] / (gamma - 1)) + 0.5*ro[i] * pow(u[i], 2);

				F_plus.F1 = ro[i] * u[i];
				F_plus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_plus.F3 = u[i] * (p[i] + E);
			}
			else {
				Variables U_L = { ro[i], u[i], p[i] };
				Variables U_R = { ro[i + 1], u[i + 1], p[i + 1] };
				Variables U_res;
				U_res = solver_riman(U_L, U_R);

				E = (U_res.p / (gamma - 1)) + 0.5*U_res.ro*pow(U_res.u, 2);

				F_plus.F1 = U_res.ro*U_res.u;
				F_plus.F2 = U_res.p + U_res.ro*pow(U_res.u, 2);
				F_plus.F3 = U_res.u*(U_res.p + E);
			}

			U_new[i].U1 = U[i].U1 - (tau / (h))*(F_plus.F1 - F_minus.F1);
			U_new[i].U2 = U[i].U2 - (tau / (h))*(F_plus.F2 - F_minus.F2);
			U_new[i].U3 = U[i].U3 - (tau / (h))*(F_plus.F3 - F_minus.F3);
		
	}

	for (int i = 1; i< n1-1; i++) {
		ro[i] = U_new[i].U1;
		u[i] = U_new[i].U2 / U_new[i].U1;
		p[i] = (gamma - 1.0)*(U_new[i].U3 - (pow(U_new[i].U2, 2)) / (2.0*U_new[i].U1));
	}
}


void Godunov::setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p,
								   double x0, Variables &left, Variables &right) {

	//проверка что все массивы одной длины
	for (int i = 0; i < p.length(); i++) {
		if ((minX + (h * i)) < x0) {
			p[i] = left.p;
			ro[i] = left.ro;
			u[i] = left.u;
		}
		else {
			p[i] = right.p;
			ro[i] = right.ro;
			u[i] = right.u;
		}
	}
}


void Godunov::setBoundaryConditions(Array<double> &ro, Array<double> &u, Array<double> &p, 
	                                Variables &left, Variables &right) {
	ro.setFirst(left.ro);
	u.setFirst(left.u);
	p.setFirst(left.p);

	ro.setLast(right.ro);
	u.setLast(right.u);
	p.setLast(right.p);
}


void Godunov::solve() {
	double gamma = 1.4;

	const int arrSize = N + 1;

	Array<double> ro(arrSize);
	Array<double> u(arrSize);
	Array<double> p(arrSize);

	setInitialConditions(ro, u, p, x0, bcLeft, bcRight);

	double t = 0;
	double coef = 0.3;
	int k = 0;
	while (t <= maxT) {
		solver(ro, u, p, arrSize);
		//граничные услови§
		setBoundaryConditions(ro, u, p, bcLeft, bcRight);

		//подбор шага по шагу
		t += tau;
		k++;
		tau = min(tau, coef * h / (max(fabs(leftD), fabs(rightD))));
		if (tau == -1) tau = min(tau, coef * h / (u.max()));

		cout << "\n \t T: " << t << "\t tau: " << tau << "\t h: " << h << "\n\n";
	}
	cout << "\n \t T: " << t << "\n\n\n";

	//----------------------Testing----------------------------------------//
	(new SolutionWriter())->write("Godunov", ro, u, p);

	ro = NULL; u = NULL; p = NULL;
}

Godunov::~Godunov()
{
}

