#include<iostream>
#include <cmath>
#include <stdlib.h>
#include"Solver.h"
#include "ASUM.h"
#include "Godunov.h"


using namespace std;

int main() {
	std::string name_file = "input\\input.txt";

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
		}
		else if (ch == 0) {
			Godunov a;
			ptr = &a;
			ptr->read_file(name_file);
			ptr->solve();
		}
		else if (ch == 1) {
			ASUM a;
			ptr = &a;
			ptr->read_file(name_file);
			ptr->solve();
		}
		else if (ch == 3) {
			return 0;
		}
		else cout << "error!";
		cout << "\n\n";
		ch = 5;
	}

	system("pause");
	return 0;
}