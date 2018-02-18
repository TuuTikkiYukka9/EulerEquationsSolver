#include "Solver.h"
#include "strConvert.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>




Solver::Solver()
{
}

void Solver::read_file(std::string name_file) {
	using namespace std;
	//std::string name_file = "in.txt";
	ifstream filer(name_file);

	string str;
	if (!filer.is_open()) std::cout << "Can not open file! \n";
	while (getline(filer, str) && str.size() != 0) {
		token(str);

	}
	filer.close();
}

void Solver::token(std::string str) {
	int size = str.size();
	std::string buf="";
	std::vector<std::string> ls;
	for (int i = 0;i < size; i++) {
		if (str[i] == ' ') {
			ls.push_back(buf);
			buf = "";
		}
		else{ buf += str[i]; }
	}
	ls.push_back(buf);


	if (ls.size() == 3) {
		if (ls[0] == "x_0") x0 = to_double(ls[2]);
		if (ls[0] == "l.ro") kr_l.ro = to_double(ls[2]);
		if (ls[0] == "l.u") kr_l.u = to_double(ls[2]);
		if (ls[0] == "l.p") kr_l.p = to_double(ls[2]);
		if (ls[0] == "r.ro") kr_r.ro = to_double(ls[2]);
		if (ls[0] == "r.u") kr_r.u = to_double(ls[2]);
		if (ls[0] == "r.p") kr_r.p = to_double(ls[2]);
	}
}


void Solver::print(std::string name_metod, double *ro, double *u, double *p, int n1) {
	using namespace std;

	std::string name_file_u, name_file_p, name_file_ro;
	
	name_file_u = "u_"+ name_metod +".xls";
	name_file_p = "p_" + name_metod + ".xls";
	name_file_ro = "ro_" + name_metod + ".xls";

	ofstream filew1(name_file_u);
	ofstream filew2(name_file_p);
	ofstream filew3(name_file_ro);
	filew1.setf(ios::fixed);
	filew1 << setprecision(5);
	filew2.setf(ios::fixed);
	filew2 << setprecision(5);
	filew3.setf(ios::fixed);
	filew3 << setprecision(5);

	for (int i = 0; i < n1; i++) {
		filew1 << setw(8) << right << u[i] << "\t";
		filew2 << setw(8) << right << p[i] << "\t";
		filew3 << setw(8) << right << ro[i] << "\t";
	}
	filew1 << '\n';
	filew2 << '\n';
	filew3 << '\n';

	filew1.close();
	filew2.close();
	filew3.close();

	convert(name_file_u);
	convert(name_file_p);
	convert(name_file_ro);

}

void Solver::convert(std::string name_file) {
	using namespace std;

	ifstream filer(name_file);
	string str;
	getline(filer, str);
	filer.close();
	int size = str.size();
	for (int i = 0;i < size;i++) {
		if (str[i] == '.') str[i] = ',';
	}

	ofstream filew(name_file);
	filew << str;
	filew.close();

}

Solver::~Solver()
{
}
