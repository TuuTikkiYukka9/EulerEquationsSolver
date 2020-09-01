#pragma once
#include "isolutionwriter.h"
#include <string>
#include "array.h"

class SolutionWriter : public ISolutionWriter {
private:
	std::string convertDoubleToString(double x) {
		std::string str = std::to_string(x);
		for (std::size_t i = 0; i < str.size(); i++) {
			if (str[i] == '.') str[i] = ',';
		}
		return str;
	};

	void writeFile(std::string fileName, Array<double> &values) {
		using namespace std;
		
		ofstream fileStream(fileName);
		
		for (int i = 0; i < values.length(); i++) {
			fileStream << setw(8) << right << convertDoubleToString(values[i]) << "\t";
		}
		fileStream << '\n';
		fileStream.close();		
	};

public:
	void write(std::string methodName, Array<double> &ro, Array<double> &u, Array<double> &p) {
		using namespace std;

		std::string uFileName, pFileName, roFilrName;

		uFileName = "u_" + methodName + ".xls";
		pFileName = "p_" + methodName + ".xls";
		roFilrName = "ro_" + methodName + ".xls";

		writeFile(uFileName, u);
		writeFile(pFileName, p);
		writeFile(roFilrName, ro);
	};

};
