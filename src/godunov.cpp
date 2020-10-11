#include "godunov.h"
#include "array.h"


Godunov::Godunov(void) {
	N = 0;
	K = 0;
	maxT = 0;
	h = 0;
	tau = 0;
}


Response Godunov::initСomputationalGrid(const СomputationalGrid &grid, double maxTime) {
	N = grid.numberOfXSplits;
	maxT = maxTime;
	h = (maxX - minX) / N;
	tau = 0.411 * h;
	K = (int)(maxT / tau);
	return Response{ true, "" };
}

double acceleration(const Variables<double> &var) {
	double gamma = 1.4;
	return pow(gamma * (var.p / var.ro), 0.5);
}

Variables<double> averageValue(const Discontinuity<Variables<double>> &u) {
	return { (u.left.ro + u.right.ro) / 2,
			 (u.left.u + u.right.u) / 2,
			 (u.left.p + u.right.p) / 2 };
}


Variables<double> Godunov::solveRiemannProblem(const Discontinuity<Variables<double>> &u) {
	
	double gamma = 1.4;
	Variables<double> riemannPoblemSolution;
	Discontinuity<double> c = { 1, 1 };
	Discontinuity<double> a;
	Variables<double> var;

	a = { acceleration(u.left), acceleration(u.right) };
	var.p = (a.right * u.right.ro * u.left.p 
		        + a.left * u.left.ro * u.right.p 
		        + a.right * u.right.ro * a.left * u.left.ro * (u.left.u - u.right.u)
		     ) / (a.right * u.right.ro + a.left * u.left.ro);
	var.u = (a.right * u.right.ro * u.right.u 
		         + a.left * u.left.ro * u.left.u 
		         + u.left.p - u.right.p
		     ) / (a.right * u.right.ro + a.left * u.left.ro);

	Variables<double> buffer = { 0, 3000 , 3 };
	double eps = 0.000001;
	while ((fabs(var.p - buffer.p)) >= eps || (fabs(var.u - buffer.u)) >= eps) {
		buffer = var;

		if (var.p < u.left.p) {
			c.left = a.left * u.left.ro * (gamma - 1) * (1 - (var.p / u.left.p)) 
				     / (2 * gamma * (1 - pow(var.p / u.left.p, (gamma - 1) / (2 * gamma))));
			//cout<<"veer voln razr vlevo \n";
			leftD = u.left.u - c.left / u.left.ro;
		}
		else if (var.p > u.left.p) {
			c.left = sqrt(u.left.ro * ((gamma + 1) * var.p + (gamma - 1) * u.left.p) / 2.0);
			//cout<<"Udarnaia volna vlevo\n";
			leftD = u.left.u - a.left;
		}

		if (var.p < u.right.p) {
			c.right = a.right * u.right.ro * (gamma - 1) * (1 - var.p / u.right.p) 
				      / (2 * gamma*(1 - pow(var.p / u.right.p, (gamma - 1) / (2 * gamma))));
			//cout<<"veer voln razr vpravo\n";
			rightD = u.right.u + c.right / u.right.ro;
		}
		else if (var.p > u.right.p) {
			c.right = sqrt(u.right.ro * ((gamma + 1) * var.p + (gamma - 1) * u.right.p) / 2.0);
			//cout<<"Udarnaia volna vpravo \n";
			rightD = u.right.u + a.right;
		}

		var.p = (c.right * u.left.p + c.left * u.right.p 
			        + c.right * c.left * (u.left.u - u.right.u)) 
		        / (c.right + c.left);
		var.u = (c.left * u.left.u + c.right * u.right.u + u.left.p - u.right.p) / (c.right + c.left);
	}

	//Считаем ро
	var.ro = averageValue(u).ro;
	if (var.u < 0) {
		if (var.p > u.right.p) {//ударная волна
			var.ro = u.right.ro * ((gamma + 1) * var.p + (gamma - 1) * u.right.p) 
				     / ((gamma - 1) * var.p + (gamma + 1) * u.right.p);
		}
		if (var.p < u.right.p) {//веер волн
			double a_new_R = a.right - ((gamma - 1) / 2) * (u.right.u - var.u);
			var.ro = gamma * var.p / pow(a_new_R, 2);
		}
	}
	else if (var.u > 0) {
		if (var.p > u.left.p) {//ударная волна
			var.ro = u.left.ro * ((gamma + 1) * var.p + (gamma - 1) * u.left.p) 
				     / ((gamma - 1) * var.p + (gamma + 1) * u.left.p);
		}
		if (var.p < u.left.p) {//веер волн
			double a_new_L = a.left + ((gamma - 1) / 2) * (u.left.u - var.u);
			var.ro = gamma * var.p / pow(a_new_L, 2);
		}
	}

	if (var.p > 0) { //проверяем
		riemannPoblemSolution = var;
	}
	else {//Если конфигурация разрывов не выявленна
		riemannPoblemSolution = averageValue(u);
		leftD = riemannPoblemSolution.u;
		rightD = riemannPoblemSolution.u;
	}
	return riemannPoblemSolution;
}


void Godunov::solver(Array<double> &ro, Array<double> &u, Array<double> &p, int n1) {

	if (maxT == 0) return;

	double gamma = 1.4;
	Array<ConservativeVariables> U(n1);
	Array<ConservativeVariables> U_new(n1); 
	Array<Flow> F(n1);
	
	Array<ConservativeVariables> U_drob(n1);
	Array<Flow> F_drob(n1);

	for (int i = 0; i < n1; i++) {
		U[i].U1 = ro[i];
		U[i].U2 = ro[i] * u[i];
		double e = ro[i] * ((1.0 / (1.4 - 1.0)) * (p[i] / ro[i]));
		U[i].U3 = (p[i] / (gamma - 1)) + 0.5 * ro[i] * pow(u[i], 2);//

	}

	for (int i = 1; i<n1 - 1; i++) {
			struct Flow F_plus, F_minus;
			double E;
			if (p[i - 1] == p[i] && ro[i - 1] == ro[i] && u[i - 1] == u[i]) {

				E = (p[i] / (gamma - 1)) + 0.5 * ro[i] * pow(u[i], 2);

				F_minus.F1 = ro[i] * u[i];
				F_minus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_minus.F3 = u[i] * (p[i] + E);
			}
			else {
				Variables<double> U_L = { ro[i - 1], u[i - 1], p[i - 1] };
				Variables<double> U_R = { ro[i], u[i], p[i] };
				Variables<double> U_res;
				U_res = solveRiemannProblem({ U_L, U_R });

				E = (U_res.p / (gamma - 1)) + 0.5 * U_res.ro * pow(U_res.u, 2);

				F_minus.F1 = U_res.ro * U_res.u;
				F_minus.F2 = U_res.p + U_res.ro * pow(U_res.u, 2);
				F_minus.F3 = U_res.u * (U_res.p + E);


			}

			if (p[i] == p[i + 1] && ro[i] == ro[i + 1] && u[i] == u[i + 1]) {
				E = (p[i] / (gamma - 1)) + 0.5 * ro[i] * pow(u[i], 2);

				F_plus.F1 = ro[i] * u[i];
				F_plus.F2 = p[i] + ro[i] * pow(u[i], 2);
				F_plus.F3 = u[i] * (p[i] + E);
			}
			else {
				Variables<double> U_L = { ro[i], u[i], p[i] };
				Variables<double> U_R = { ro[i + 1], u[i + 1], p[i + 1] };
				Variables<double> U_res;
				U_res = solveRiemannProblem({ U_L, U_R });

				E = (U_res.p / (gamma - 1)) + 0.5 * U_res.ro * pow(U_res.u, 2);

				F_plus.F1 = U_res.ro * U_res.u;
				F_plus.F2 = U_res.p + U_res.ro * pow(U_res.u, 2);
				F_plus.F3 = U_res.u * (U_res.p + E);
			}

			U_new[i].U1 = U[i].U1 - (tau / (h))*(F_plus.F1 - F_minus.F1);
			U_new[i].U2 = U[i].U2 - (tau / (h))*(F_plus.F2 - F_minus.F2);
			U_new[i].U3 = U[i].U3 - (tau / (h))*(F_plus.F3 - F_minus.F3);
		
	}

	for (int i = 1; i< n1-1; i++) {
		ro[i] = U_new[i].U1;
		u[i] = U_new[i].U2 / U_new[i].U1;
		p[i] = (gamma - 1.0) * (U_new[i].U3 - (pow(U_new[i].U2, 2)) / (2.0 * U_new[i].U1));
	}
}


void Godunov::setInitialConditions(Array<double> &ro, Array<double> &u, Array<double> &p,
								   double x0, Variables<double> &left, Variables<double> &right) {

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
	                                Variables<double> &left, Variables<double> &right) {
	ro.setFirst(left.ro);
	u.setFirst(left.u);
	p.setFirst(left.p);

	ro.setLast(right.ro);
	u.setLast(right.u);
	p.setLast(right.p);
}


Variables<Array<double>*> Godunov::solve() {
	double gamma = 1.4;

	const int arrSize = N + 1;

	Array<double> *ro = new Array<double>(arrSize);
	Array<double> *u = new Array<double>(arrSize);
	Array<double> *p = new Array<double>(arrSize);

	setInitialConditions(*ro, *u, *p, x0, bcLeft, bcRight);

	double t = 0;
	double coef = 0.3;
	int k = 0;
	while (t <= maxT) {
		solver(*ro, *u, *p, arrSize);
		//граничные услови¤
		setBoundaryConditions(*ro, *u, *p, bcLeft, bcRight);

		//подбор шага по шагу
		t += tau;
		k++;
		tau = min(tau, coef * h / (max(fabs(leftD), fabs(rightD))));
		if (tau == -1) tau = min(tau, coef * h / (u->max()));

		std::cout << "\n \t T: " << t << "\t tau: " << tau << "\t h: " << h << "\n\n";
	}
	std::cout << "\n \t T: " << t << "\n\n\n";

	return { ro, u, p };
}

Godunov::~Godunov() {
}