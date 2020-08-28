#pragma once
#include "isolutionwriter.h"
#include <string>

class SolutionWriter : public ISolutionWriter {
private:
	std::string convertDoubleToString(double x) {
		std::string str = std::to_string(x);
		for (std::size_t i = 0; i < str.size(); i++) {
			if (str[i] == '.') str[i] = ',';
		}
		return str;
	};

	void writeFile(std::string fileName, double *values, int size) {
		using namespace std;
		
		ofstream fileStream(fileName);
		
		for (int i = 0; i < size; i++) {
			fileStream << setw(8) << right << convertDoubleToString(values[i]) << "\t";
		}
		fileStream << '\n';
		fileStream.close();		
	};

public:
	void write(std::string methodName, double *ro, double *u, double *p, int size) {
		using namespace std;

		std::string uFileName, pFileName, roFilrName;

		uFileName = "u_" + methodName + ".xls";
		pFileName = "p_" + methodName + ".xls";
		roFilrName = "ro_" + methodName + ".xls";

		writeFile(uFileName, u, size);
		writeFile(pFileName, p, size);
		writeFile(roFilrName, ro, size);
	};

};
