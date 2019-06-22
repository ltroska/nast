#ifndef NAST_SOLVERS_JACOBI_HPP_
#define NAST_SOLVERS_JACOBI_HPP_

#include "solvers/solver_base.hpp"

#include <iostream>

namespace nast { namespace solvers {
	
class jacobi : public solver_base {
public:
	virtual void solve(grid::staggered_grid& grid, const parameters::parameters& params) override
	{
		auto dx_sq = grid.get_dx_sq();
		auto dy_sq = grid.get_dy_sq();
		
		auto old_p = grid.p;
		
		for (auto& id : grid.fluid_cells)
		{
			auto& i = id.first;
			auto& j = id.second;
		
			grid.p(i, j) =
                         ( (old_p(i + 1, j) + old_p(i - 1, j)) * dy_sq
                            + (old_p(i, j + 1) + old_p(i, j - 1)) * dx_sq
                            - dx_sq * dy_sq * grid.rhs(i, j))
                        /
                        (2 * (dx_sq + dy_sq));
		}
	}
};
	
} //namespace solvers
} //namespace nast

#endif
