#ifndef NAST_SOLVERS_SOLVER_BASE_HPP_
#define NAST_SOLVERS_SOLVER_BASE_HPP_

#include "grid/staggered_grid.hpp"
#include "parameters/parameters.hpp"

#include <cmath>

namespace nast { namespace solvers {
	
class solver_base {
public:
	virtual void solve(grid::staggered_grid& grid, const parameters::parameters& params) = 0;
	
	Real compute_residual(grid::staggered_grid& grid)
	{		
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto over_dx_sq = 1./grid.get_dx_sq();
		auto over_dy_sq = 1./grid.get_dy_sq();
		
		Real tmp;
		Real residual = 0;
		
		for (auto& id : grid.fluid_cells)
		{
			auto& i = id.first;
			auto& j = id.second;

			tmp =
				(grid.p(i + 1, j) - 2 * grid.p(i, j) + grid.p(i - 1, j)) * over_dx_sq
				+ (grid.p(i, j + 1) - 2 * grid.p(i, j) + grid.p(i, j - 1)) * over_dy_sq
				- grid.rhs(i, j);

			residual += std::pow(tmp, 2);
        }
             						       
        return residual / ((size_x - 1) * (size_y - 1)) ;
	}
};
	
} //namespace solvers
} //namespace nast

#endif
