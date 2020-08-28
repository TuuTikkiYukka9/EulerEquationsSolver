#pragma once
#include "solver.h"
#include "iequationsystemreader.h"
#include <string>
#include <vector>


class EquationSystemReader : public IEquationSystemReader {
private:
	std::vector<std::string> split(std::string str, std::string delimiter) {
		std::size_t pos;
		std::vector<std::string> tokens;
		while ((pos = str.find(delimiter)) != std::string::npos) {
			tokens.push_back(str.substr(0, pos));
			str.erase(0, pos + delimiter.length());
		}
		tokens.push_back(str);
		return tokens;
	}

	void parseLine(std::string line, EquationSystem& eq) {
		std::vector<std::string> tokens = split(line, " = ");
		if (tokens.size() >= 2) {
			std::string paramName = tokens[0];
			if (paramName == "x_0") eq.x0 = std::stod(tokens[1]);
			else if (paramName == "l.ro") eq.leftBoundaryCondition.ro = std::stod(tokens[1]);
			else if (paramName == "l.u") eq.leftBoundaryCondition.u = std::stod(tokens[1]);
			else if (paramName == "l.p") eq.leftBoundaryCondition.p = std::stod(tokens[1]);
			else if (paramName == "r.ro") eq.rightBoundaryCondition.ro = std::stod(tokens[1]);
			else if (paramName == "r.u") eq.rightBoundaryCondition.u = std::stod(tokens[1]);
			else if (paramName == "r.p") eq.rightBoundaryCondition.p = std::stod(tokens[1]);
			else if (paramName == "min_x") eq.minX = std::stod(tokens[1]);
			else if (paramName == "max_x") eq.maxX = std::stod(tokens[1]);
		}
	}

public:
	EquationSystem readFile(std::string fileName) {
		using namespace std;
		EquationSystem eq = EquationSystem();
		ifstream filer(fileName);

		string line;
		if (!filer.is_open()) std::cout << "Can not open file! \n";
		while (getline(filer, line) && line.size() != 0) {
			parseLine(line, eq);
		}
		filer.close();
		return eq;
	}
};