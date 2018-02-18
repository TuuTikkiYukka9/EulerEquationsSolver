#include "Godunov.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

/*Не должен конструктор ничего печатать, начальные и краевые услови§ то же надо передавать 4 вектора по 3 элемента, и еще положение разрыва х0, до звуковое течение*/
Godunov::Godunov(void)
{
	min_x = 0.0;
	max_x = 1.0;

	cout << "Enter the number of partitions for x:";
	cin >> N;
	cout << "Enter max T:";
	cin >> max_T;

	h = (max_x - min_x) / N;
	tau = 0.411*h;
	K = (int)(max_T / tau);

}

Perem Godunov::solver_riman(Perem U_L, Perem U_R) {
	double gam = 1.4;
	struct Perem U;
	double c_L = 1;
	double c_R = 1;
	double Ro;

	double a_L = pow(gam*(U_L.p / U_L.ro), 0.5);
	double a_R = pow(gam*(U_R.p / U_R.ro), 0.5);
	double P0 = (a_R*U_R.ro*U_L.p + a_L*U_L.ro*U_R.p + a_R*U_R.ro*a_L*U_L.ro*(U_L.u - U_R.u)) / (a_R*U_R.ro + a_L*U_L.ro);
	double U0 = (a_R*U_R.ro*U_R.u + a_L*U_L.ro*U_L.u + U_L.p - U_R.p) / (a_R*U_R.ro + a_L*U_L.ro);

	double P_buf = 3;
	double U_buf = 3000;

	double eps = 0.000001;//eps=0.00000000001;
	while ((fabs(P0 - P_buf)) >= eps || (fabs(U0 - U_buf)) >= eps) {
		P_buf = P0;
		U_buf = U0;

		if (P0<U_L.p) {
			c_L = a_L * U_L.ro*(gam - 1)*(1 - (P0 / U_L.p)) / (2 * gam*(1 - pow(P0 / U_L.p, (gam - 1) / (2 * gam))));
			//cout<<"veer voln razr vlevo \n";
			D_L = U_L.u - c_L / U_L.ro;
		}
		else if (P0>U_L.p) {
			c_L = sqrt(U_L.ro*((gam + 1)*P0 + (gam - 1)*U_L.p) / 2.0);
			//cout<<"Udarnaia volna vlevo\n";
			D_L = U_L.u - a_L;
		}


		if (P0<U_R.p) {
			c_R = a_R * U_R.ro*(gam - 1)*(1 - P0 / U_R.p) / (2 * gam*(1 - pow(P0 / U_R.p, (gam - 1) / (2 * gam))));
			//cout<<"veer voln razr vpravo\n";
			D_R = U_R.u + c_R / U_R.ro;
		}
		else if (P0>U_R.p) {
			c_R = sqrt(U_R.ro*((gam + 1)*P0 + (gam - 1)*U_R.p) / 2.0);
			//cout<<"Udarnaia volna vpravo \n";
			D_R = U_R.u + a_R;
		}

		P0 = (c_R * U_L.p + c_L *U_R.p + c_R * c_L * (U_L.u - U_R.u)) / (c_R + c_L);
		U0 = (c_L * U_L.u + c_R * U_R.u + U_L.p - U_R.p) / (c_R + c_L);
	}

	//Ччитаем ро
	Ro = (U_L.ro + U_R.ro) / 2;
	if (U0<0) {
		if (P0>U_R.p) {//ударна§ волна
			Ro = U_R.ro*((gam + 1)*P0 + (gam - 1)*U_R.p) / ((gam - 1)*P0 + (gam + 1)*U_R.p);
		}
		if (P0<U_R.p) {//веер волн
			double a_new_R = a_R - ((gam - 1) / 2)*(U_R.u - U0);
			Ro = gam*P0 / pow(a_new_R, 2);
		}
	}
	else if (U0>0) {
		if (P0>U_L.p) {//ударна§ волна
			Ro = U_L.ro*((gam + 1)*P0 + (gam - 1)*U_L.p) / ((gam - 1)*P0 + (gam + 1)*U_L.p);
		}
		if (P0<U_L.p) {//веер волн
			double a_new_L = a_L + ((gam - 1) / 2)*(U_L.u - U0);
			Ro = gam*P0 / pow(a_new_L, 2);
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
		D_L = U.u;//04.01
		D_R = U.u;
	}

	return U;
}

void Godunov::solver(double *ro, double *u, double *p, int n1) {
	double gam = 1.4;
	struct consv *U;
	struct consv *U_new;
	struct potoc *F;
	F = new struct potoc[n1];
	U = new struct consv[n1];
	U_new = new struct consv[n1];

	struct consv *U_drob;
	struct potoc *F_drob;
	F_drob = new struct potoc[n1];
	U_drob = new struct consv[n1];

	for (int i = 0; i<n1; i++) {

		U[i].U1 = ro[i];
		U[i].U2 = ro[i] * u[i];
		double e = ro[i] * ((1.0 / (1.4 - 1.0))*(p[i] / ro[i]));
		U[i].U3 = (p[i] / (gam - 1)) + 0.5*ro[i] * pow(u[i], 2);//

	}

	for (int i = 1; i<n1 - 1; i++) {
			struct potoc F_plus, F_minus;
			double E;
			if (p[i - 1] == p[i] && ro[i - 1] == ro[i] && u[i - 1] == u[i]) {

				E = (p[i] / (gam - 1)) + 0.5*ro[i] * pow(u[i], 2);

				F_minus.F1 = ro[i] * u[i];
				F_minus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_minus.F3 = u[i] * (p[i] + E);
			}
			else {
				Perem U_L = { ro[i - 1], u[i - 1], p[i - 1] };
				Perem U_R = { ro[i], u[i], p[i] };
				Perem U_res;
				U_res = solver_riman(U_L, U_R);

				E = (U_res.p / (gam - 1)) + 0.5*U_res.ro*pow(U_res.u, 2);

				F_minus.F1 = U_res.ro*U_res.u;
				F_minus.F2 = U_res.p + U_res.ro*pow(U_res.u, 2);
				F_minus.F3 = U_res.u*(U_res.p + E);


			}

			if (p[i] == p[i + 1] && ro[i] == ro[i + 1] && u[i] == u[i + 1]) {
				E = (p[i] / (gam - 1)) + 0.5*ro[i] * pow(u[i], 2);

				F_plus.F1 = ro[i] * u[i];
				F_plus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_plus.F3 = u[i] * (p[i] + E);
			}
			else {
				Perem U_L = { ro[i], u[i], p[i] };
				Perem U_R = { ro[i + 1], u[i + 1], p[i + 1] };
				Perem U_res;
				U_res = solver_riman(U_L, U_R);

				E = (U_res.p / (gam - 1)) + 0.5*U_res.ro*pow(U_res.u, 2);

				F_plus.F1 = U_res.ro*U_res.u;
				F_plus.F2 = U_res.p + U_res.ro*pow(U_res.u, 2);
				F_plus.F3 = U_res.u*(U_res.p + E);
			}


			U_new[i].U1 = U[i].U1 - (tau / (h))*(F_plus.F1 - F_minus.F1);
			U_new[i].U2 = U[i].U2 - (tau / (h))*(F_plus.F2 - F_minus.F2);
			U_new[i].U3 = U[i].U3 - (tau / (h))*(F_plus.F3 - F_minus.F3);
		
	}

	for (int i = 1; i<n1 - 1; i++) {
		ro[i] = U_new[i].U1;
		u[i] = U_new[i].U2 / U_new[i].U1;
		p[i] = (gam - 1.0)*(U_new[i].U3 - (pow(U_new[i].U2, 2)) / (2.0*U_new[i].U1));
	}


	delete[]U;
	delete[]U_new;
	delete[]F;
	delete[]U_drob;
	delete[]F_drob;
	ro = NULL; u = NULL; p = NULL;
	U = NULL; F = NULL;
	U_drob = NULL; F_drob = NULL;
	U_new = NULL;
}

void Godunov::solve() {
	double gam = 1.4;

	const int n1 = N + 1;


	double *ro;
	double *u;
	double *p;
	ro = new double[n1];
	u = new double[n1];
	p = new double[n1];


	for (int i = 0; i < n1; i++) {
		ro[i] = 0;
		p[i] = 0;
		u[i] = 0;
	}


	//начальное усдловие
	for (int i = 0; i < n1; i++) {
		if ((min_x + (h*i))<x0) {// Уут точно х0
			p[i] = kr_l.p;
			ro[i] = kr_l.ro;
			u[i] = kr_l.u;
		}
		else {
			p[i] = kr_r.p;//ro
			ro[i] = kr_r.ro;// ro8u
			u[i] = kr_r.u;//E
		}
	}
	


	double T_vsego = 0;
	int k = 0;
	//for(int k = 0; k < (max_T/tau); k++){
	while (T_vsego <= max_T) {// && k < 100000000){
		solver(ro, u, p, n1);
		//граничные услови§
		ro[0] = kr_l.ro;
		u[0] = kr_l.u;
		p[0] = kr_l.p;
		ro[n1 - 1] = kr_r.ro;
		u[n1 - 1] = kr_r.u;
		p[n1 - 1] = kr_r.p;

		//подбор шага по шагу
		T_vsego += tau;
		k++;
		tau = min(tau, 0.3*h / (max(fabs(D_L), fabs(D_R))));
		//tau=min(tau, 0.3*h/(max(fabs(a_L), fabs(a_R))));
		//tau=min(tau, 0.3*h/(max(u,n1)));
		if (tau == -1) {
			//cout<<"\n \t T: "<<T_vsego<<"\t tau: "<<tau<<"\n\n\n";
			//break;
			tau = min(tau, 0.3*h / (max(u, n1)));
			//tau = 0.3*h;
		}
		cout << "\n \t T: " << T_vsego << "\t tau: " << tau << "\t h: " << h << "\n\n";
	}
	cout << "\n \t T: " << T_vsego << "\n\n\n";

	//----------------------Testing----------------------------------------//
	this->print("Godunov", ro, u, p, n1);


	delete[]ro;
	delete[]u;
	delete[]p;
	ro = NULL; u = NULL; p = NULL;
}

Godunov::~Godunov()
{
}

