#ifndef NAST_SOLVERS_SOR_HPP_
#define NAST_SOLVERS_SOR_HPP_

#include "solvers/solver_base.hpp"

namespace nast { namespace solvers {
	
class sor : public solver_base {
public:
	virtual void solve(grid::staggered_grid& grid, const parameters::parameters& params) override
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto dx_sq = grid.get_dx_sq();
		auto dy_sq = grid.get_dy_sq();
		
		auto omega = params.omega;
		
		auto part1 = 1. - omega;
        auto part2 = omega / ( 2 * (1/dx_sq + 1/dy_sq));
        
    //    std::cout << "parts " << part1 << " " << part2 << " " << 1./dx_sq << " " << 1./dy_sq << std::endl;
		
		//std::cout << grid.rhs << std::endl;
		
		//std::cout << std::endl;
		//	std::cout << grid.p << std::endl;
		
		for (std::size_t j = 1; j < size_y - 1; ++j)
		{
			for (std::size_t i = 1; i < size_x - 1; ++i)
			{
				auto const& type = grid.cell_type(i, j);
				
				if (type[is_fluid])
				{
					grid.p(i, j) =
						part1 * grid.p(i, j)
						+ part2 * (
								(grid.p(i + 1, j) + grid.p(i - 1, j)) / dx_sq
								+ (grid.p(i, j + 1) + grid.p(i, j - 1)) / dy_sq
								- grid.rhs(i, j)
						);
				}
			}
		}
	}
};
	
} //namespace solvers
} //namespace nast

#endif
