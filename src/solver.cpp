#include "Solver.h"

void Solver::init(EquationSystem &eq) {
	bcLeft = eq.leftBoundaryCondition;
	bcRight = eq.rightBoundaryCondition;
	x0 = eq.x0;
	maxX = eq.maxX;
	minX = eq.minX;
};

