#include "Solver.h"

void Solver::init(EquationSystem &eq) {
	kr_l = eq.leftBoundaryCondition;
	kr_r = eq.rightBoundaryCondition;
	x0 = eq.x0;
	max_x = eq.maxX;
	min_x = eq.minX;
};

