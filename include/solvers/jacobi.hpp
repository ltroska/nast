#ifndef NAST_SOLVERS_JACOBI_HPP_
#define NAST_SOLVERS_JACOBI_HPP_

#include "solvers/solver_base.hpp"

#include <iostream>

namespace nast { namespace solvers {
	
class jacobi : public solver_base {
public:
	virtual void solve(grid::staggered_grid& grid, const parameters::parameters& params) override
	{
		std::cout << "Jacobi" << std::endl;
	}
};
	
} //namespace solvers
} //namespace nast

#endif
