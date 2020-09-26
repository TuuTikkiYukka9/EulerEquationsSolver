#include "asum.h"
#include "godunov.h"
#include "equationsystemreader.h"
#include "solutionwriter.h"


int main() {

	std::string name_file = "input\\input.txt";
	IEquationSystemReader *equationSystemReader = new EquationSystemReader();
	EquationSystem eq = equationSystemReader->readFile(name_file);
	int ch = 5;
	while (ch != 3) {
		std::cout << "Menu: \n";
		std::cout << "0 - Godunov \n";
		std::cout << "1 - ASUM \n";
		std::cout << "2 - input  file name \n";
		std::cout << "3 - Exit \n";
		std::cout << "Enter:";
		std::cin >> ch;
		
		Solver *ptr;
		
		if (ch == 0) {
			ÑomputationalGrid grid;
			double maxT;
			std::cout << "Enter the number of partitions for x:";
			std::cin >> grid.numberOfXSplits;
			std::cout << "Enter max T:";
			std::cin >> maxT;

			Godunov a;
			Variables<Array<double>*> result;
			ptr = &a;
			ptr->init(eq);
			if ((ptr->initÑomputationalGrid(grid, maxT)).success) {
				result = ptr->solve();
				(new SolutionWriter())->write("Godunov", result);
				delete result.p; delete result.ro, delete result.u;
			}
		}
		else if (ch == 1) {
			ÑomputationalGrid grid;
			double maxT;

			std::cout << "Enter the number of partitions for x:";
			std::cin >> grid.numberOfXSplits;
			std::cout << "Enter the number of partitions for t:";
			std::cin >> grid.numberOfTimeSplits;
			std::cout << "Enter max T:";
			std::cin >> maxT;

			ASUM a;
			Variables<Array<double>*> result;
			ptr = &a;
			ptr->init(eq);
			const Response resp = ptr->initÑomputationalGrid(grid, maxT);
			if (resp.success) {
				result = ptr->solve();
				(new SolutionWriter())->write("AUSM", result);
				delete result.p; delete result.ro, delete result.u;
			}
			else {
				std::cout << resp.message;
			}
		}
		else if (ch == 2) {
			Variables<Array<double>*> result;

			std::cin >> name_file;
			std::cout << "\n file name: " << name_file << "\n";
			eq = equationSystemReader->readFile(name_file);
			
		}
		else if (ch == 3) {
			return 0;
		}
		else {
			std::cout << "The entered item is not in the menu!";
		}
		std::cout << "\n\n";
		ch = 5;
	}
	delete equationSystemReader;
	system("pause");
	return 0;
}