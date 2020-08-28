#pragma once
#include <string>

class ISolutionWriter {
public:
	void virtual write(std::string methodName, double *ro, double *u, double *p, int size) {};
};
