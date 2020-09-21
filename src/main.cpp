#include<iostream>
#include <cmath>
#include <stdlib.h>
#include"Solver.h"
#include "ASUM.h"
#include "Godunov.h"
#include "equationsystemreader.h"
#include "matrix.h"

using namespace std;

int main() {

	std::string name_file = "input\\input.txt";
	IEquationSystemReader *equationSystemReader = new EquationSystemReader();
	EquationSystem eq = equationSystemReader->readFile(name_file);
	int ch = 5;
	while (ch != 3) {
		cout << "Menu: \n";
		cout << "0 - Godunov \n";
		cout << "1 - ASUM \n";
		cout << "2 - input  file name \n";
		cout << "3 - Exit \n";
		cout << "Enter:";
		cin >> ch;
		
		Solver *ptr;
		if (ch == 2) {
			//getline(cin, name_file);
			cin >> name_file;
			cout << "\n file name: " << name_file << "\n";
			eq = equationSystemReader->readFile(name_file);
		}
		else if (ch == 0) {
			ÑomputationalGrid grid;
			double maxT;
			cout << "Enter the number of partitions for x:";
			cin >> grid.numberOfXSplits;
            cout << "Enter max T:";
            cin >> maxT;

			Godunov a;
			ptr = &a;
			ptr->init(eq);
			if ((ptr->initÑomputationalGrid(grid, maxT)).success) {
				ptr->solve();
			}
		}
		else if (ch == 1) {
			ÑomputationalGrid grid;
			double maxT;

			cout << "Enter the number of partitions for x:";
			cin >> grid.numberOfXSplits;
			cout << "Enter the number of partitions for t:";
			cin >> grid.numberOfTimeSplits;
			cout << "Enter max T:";
			cin >> maxT;

			ASUM a;
			ptr = &a;
			ptr->init(eq);
			const Response resp = ptr->initÑomputationalGrid(grid, maxT);
			if (resp.success) {
				ptr->solve();
			}
			else {
				cout << resp.message;
			}
		}
		else if (ch == 3) {
			return 0;
		}
		else cout << "error!";
		cout << "\n\n";
		ch = 5;
	}
	delete equationSystemReader;
	system("pause");
	return 0;
}