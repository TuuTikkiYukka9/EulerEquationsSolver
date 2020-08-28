#pragma once

#include "solver.h"
#include <string>
/*#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
*/

class IEquationSystemReader {
public:
	EquationSystem virtual readFile(std::string fileName) { return EquationSystem(); }
};
