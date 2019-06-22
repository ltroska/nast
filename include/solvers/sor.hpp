#ifndef NAST_SOLVERS_SOR_HPP_
#define NAST_SOLVERS_SOR_HPP_

#include "solvers/solver_base.hpp"

namespace nast { namespace solvers {
	
class sor : public solver_base {
public:
	virtual void solve(grid::staggered_grid& grid, const parameters::parameters& params) override
	{	
		auto dx_sq = grid.get_dx_sq();
		auto dy_sq = grid.get_dy_sq();
		
		auto omega = params.omega;
		
		auto part1 = 1. - omega;
        auto part2 = omega / ( 2 * (1/dx_sq + 1/dy_sq));
		
		for (auto& id : grid.fluid_cells)
		{
			auto& i = id.first;
			auto& j = id.second;
			
			grid.p(i, j) =
				part1 * grid.p(i, j)
				+ part2 * (
						(grid.p(i + 1, j) + grid.p(i - 1, j)) / dx_sq
						+ (grid.p(i, j + 1) + grid.p(i, j - 1)) / dy_sq
						- grid.rhs(i, j)
				);		
		}
	}
};
	
} //namespace solvers
} //namespace nast

#endif
