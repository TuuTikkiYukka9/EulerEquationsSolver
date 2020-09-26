#pragma once

#include "solver.h"
#include <string>

class IEquationSystemReader {
public:
	EquationSystem virtual readFile(std::string fileName) { return EquationSystem(); }
};
