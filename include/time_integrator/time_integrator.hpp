#ifndef NAST_TIME_INTEGRATOR_TIME_INTEGRATOR_HPP_
#define NAST_TIME_INTEGRATOR_TIME_INTEGRATOR_HPP_

#include "solvers/solver_base.hpp"
#include "parameters/parameters.hpp"
#include "grid/staggered_grid.hpp"
#include "grid/boundary_conditions.hpp"

#include "util/derivatives.hpp"

#include <memory>

namespace nast { namespace time_integrator {

class time_integrator {
public:
	Real do_timestep(grid::staggered_grid& grid, std::shared_ptr<solvers::solver_base> solver, const grid::boundary_conditions& bcs, const parameters::parameters& params, Real dt)
	{
		set_boundary_conditions(grid, bcs);

		compute_forces(grid, bcs, params, dt);

		compute_rhs(grid, dt);

		Real residual = 0;
		
		std::size_t iteration = 0;
		
		Real eps_sq = std::pow(params.eps, 2);
				
		do {		
			++iteration;
	
			set_obstacle_pressure(grid, dt);

			solver->solve(grid, params);
			
			residual = solver->compute_residual(grid);	
		} while (residual > eps_sq && iteration < params.max_solver_iterations);		

		if (params.verbose)
			std::cout << "Iterations: " << iteration << " Residual: " << residual << std::endl;			

		auto max_uv = update_velocity(grid, dt);
		
		auto next_dt = params.tau * std::min(params.reynolds / 2 / ( 1 / grid.get_dx_sq() + 1 / grid.get_dy_sq()), std::min(grid.get_dx() / std::abs(max_uv.first), grid.get_dy() / std::abs(max_uv.second)));	
		
		return next_dt;
	}
		
	void set_boundary_conditions(grid::staggered_grid& grid, const grid::boundary_conditions& bcs)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto top_right_fluid_u = grid.u(size_x - 2, size_y - 2);
		auto top_right_fluid_v = grid.v(size_x - 2, size_y - 3);
		
		for (std::size_t j = 0; j < size_y; ++j)
		{
			for (std::size_t i = 0; i < size_x; ++i)
			{
				auto& type = grid.cell_type(i, j);
				
				// ignore boundary cells in the four corners
				if (type[is_boundary] && type.count() > 1)
				{					
					if (type[has_fluid_right])
					{
						switch(bcs.left_type)
						{
							case grid::boundary_condition_type::noslip:
								grid.u(i, j) = 0;
								grid.v(i, j) = 2 * bcs.left.y - grid.v(i + 1, j);
								break;

							case grid::boundary_condition_type::slip:
								grid.u(i, j) = 0;
								grid.v(i, j) = grid.v(i + 1, j);
								break;

							case grid::boundary_condition_type::outstream:
								grid.u(i, j) = grid.u(i + 1, j);
								grid.v(i, j) = grid.v(i + 1, j);
								break;

							case grid::boundary_condition_type::instream:
								grid.u(i, j) = bcs.left.x;
								grid.v(i, j) = 2 * bcs.left.y - grid.v(i + 1, j);
								break;
						}
					}
					else if (type[has_fluid_left])
					{
						switch(bcs.right_type)
						{
							case grid::boundary_condition_type::noslip:
								grid.u(i - 1, j) = 0;
								grid.v(i, j) = 2 * bcs.right.y - grid.v(i - 1, j);
								break;

							case grid::boundary_condition_type::slip:
								grid.u(i - 1, j) = 0;
								grid.v(i, j) = grid.v(i - 1, j);
								break;

							case grid::boundary_condition_type::outstream:
								grid.u(i - 1, j) = grid.u(i - 2, j);
								grid.v(i, j) = grid.v(i - 1, j);
								break;

							case grid::boundary_condition_type::instream:
								grid.u(i - 1, j) = bcs.right.x;
								grid.v(i, j) = 2 * bcs.right.y - grid.v(i - 1, j);
								break;
						}
					}
					if (type[has_fluid_top])
					{
						switch(bcs.bottom_type)
						{
							case grid::boundary_condition_type::noslip:
								grid.u(i, j) = 2 * bcs.bottom.x - grid.u(i, j + 1);
								grid.v(i, j) = 0;
								break;

							case grid::boundary_condition_type::slip:
								grid.u(i, j) = grid.u(i, j + 1);
								grid.v(i, j) = 0;
								break;

							case grid::boundary_condition_type::outstream:
								grid.u(i, j) = grid.u(i, j + 1);
								grid.v(i, j) = grid.v(i, j + 1);
								break;

							case grid::boundary_condition_type::instream:
								grid.u(i, j) = 2 * bcs.bottom.x - grid.u(i, j + 1);
								grid.v(i, j) = bcs.bottom.y;
								break;
						}
					}
					else if (type[has_fluid_bottom])
					{
						
						auto& neighbor_u = (i == size_x - 2 && j == size_y - 1 && false) ? top_right_fluid_u : grid.u(i, j - 1);
						auto& neighbor_v = (i == size_x - 2 && j == size_y - 1 && false) ? top_right_fluid_v : grid.v(i, j - 2);
			
						switch(bcs.top_type)
						{
							case grid::boundary_condition_type::noslip:
								grid.u(i, j) = 2 * bcs.top.x -  neighbor_u;
								grid.v(i, j - 1) = 0;
								break;

							case grid::boundary_condition_type::slip:
								grid.u(i, j) =  neighbor_u;
								grid.v(i, j - 1) = 0;
								break;

							case grid::boundary_condition_type::outstream:
								grid.u(i, j) =  neighbor_u;
								grid.v(i, j - 1) = neighbor_v;
								break;

							case grid::boundary_condition_type::instream:
								grid.u(i, j) = 2 * bcs.top.x -  neighbor_u;
								grid.v(i, j - 1) = bcs.top.y;
								break;
						}
					}
				}
				else if (type[is_obstacle] && type.count() > 1)
				{
                    grid.u(i, j) =
                        - grid.u(i, j - 1) * type[has_fluid_bottom]
                        - grid.u(i, j + 1) * type[has_fluid_top];

                    grid.v(i, j) =
                        - grid.v(i - 1, j) * type[has_fluid_left]
                        - grid.v(i + 1, j) * type[has_fluid_right];
                
				}
			}
		}
	}
	
	void compute_forces(grid::staggered_grid& grid, const grid::boundary_conditions& bcs, const parameters::parameters& params, Real dt)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto dx = grid.get_dx();
		auto dy = grid.get_dy();
		
		auto dx_sq = grid.get_dx_sq();
		auto dy_sq = grid.get_dy_sq();
		
		for (std::size_t j = 0; j < size_y - 1; ++j)
		{
			for (std::size_t i = 0; i < size_x - 1; ++i)
			{
				auto const& type = grid.cell_type(i, j);
				grid.f(i, j) = grid.u(i, j);
				grid.g(i, j) = grid.v(i, j);
				
				if (type[is_fluid])
				{			
					if (type[has_fluid_right])
					{
						/*std::cout << "cell " << i << " " << j << " " << dt << " " << 1./params.reynolds << " " 
						<< util::derivatives::second_derivative_fwd_bkwd_x(grid.u, i, j, dx_sq) << " " << util::derivatives::second_derivative_fwd_bkwd_y(grid.u, i, j, dy_sq)
						<< " " <<util::derivatives::first_derivative_of_square_x(grid.u, i, j, dx, params.alpha) << " " << util::derivatives::first_derivative_of_uv_y(grid.u, grid.v, i, j, dy, params.alpha) << std::endl;
						*/
						grid.f(i, j) +=
							  dt * (
								1. / params.reynolds * (util::derivatives::second_derivative_fwd_bkwd_x(grid.u, i, j, dx_sq)
											+ util::derivatives::second_derivative_fwd_bkwd_y(grid.u, i, j, dy_sq)
										  )

								- util::derivatives::first_derivative_of_square_x(grid.u, i, j, dx, params.alpha)
								- util::derivatives::first_derivative_of_uv_y(grid.u, grid.v, i, j, dy, params.alpha)
								+ bcs.external_forces.x
							);
					}
					
					if (type[has_fluid_top])
					{
						grid.g(i, j) +=
							dt * (
								1. / params.reynolds * (util::derivatives::second_derivative_fwd_bkwd_x(grid.v, i, j, dx_sq)
										   + util::derivatives::second_derivative_fwd_bkwd_y(grid.v, i, j, dy_sq)
										   )

								- util::derivatives::first_derivative_of_uv_x(grid.u, grid.v, i, j, dx, params.alpha)
								- util::derivatives::first_derivative_of_square_y(grid.v, i, j, dy, params.alpha)
								+ bcs.external_forces.y
							);
					}
				}
			}
		}		
	}

	void compute_rhs(grid::staggered_grid& grid, Real dt)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto dx = grid.get_dx();
		auto dy = grid.get_dy();
		
		for (std::size_t j = 1; j < size_y - 1; ++j)
		{
			for (std::size_t i = 1; i < size_x - 1; ++i)
			{
				auto const& type = grid.cell_type(i, j);
				
				if (type[is_fluid])
				{
					//std::cout << i << " " << j << " " 1./dt  << grid.f(i, j) << " " << grid.f(i - 1, j) << " " << 1./ dx
					//					<< " " <<grid.g(i, j) << " " << grid.g(i, j - 1) << " " << 1. / dy << std::endl;
					
					grid.rhs(i, j) =
						1. / dt *   (
										(grid.f(i, j) - grid.f(i - 1, j)) / dx
										+
										(grid.g(i, j) - grid.g(i, j - 1)) / dy
									);
				}	
			}
		}
	}
	
	void set_obstacle_pressure(grid::staggered_grid& grid, Real dt)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto dx_sq = grid.get_dx_sq();
		auto dy_sq = grid.get_dy_sq();
		
		for (std::size_t j = 0; j < size_y; ++j)
		{
			for (std::size_t i = 0; i < size_x ; ++i)
			{
				auto const& type = grid.cell_type(i, j);
				
				if (type[is_boundary])
				{				
					if (type[has_fluid_right]) {
						grid.p(i, j) = grid.p(i + 1, j);
					}
					else if (type[has_fluid_left])
						grid.p(i, j) = grid.p(i - 1, j);
					else if (type[has_fluid_top])
						grid.p(i, j) = grid.p(i, j + 1);
					else if (type[has_fluid_bottom])
						grid.p(i, j) = grid.p(i, j - 1);
				}
				else if (type[is_obstacle] && type.count() > 1)
				{
					grid.p(i, j) =
					(
						dx_sq * (grid.p(i - 1, j) * type.test(has_fluid_left)
						+ grid.p(i + 1, j) * type.test(has_fluid_right))
						+ dy_sq * (grid.p(i, j - 1) * type.test(has_fluid_bottom)
						+ grid.p(i, j + 1) * type.test(has_fluid_top))
					)
					/
					(
						dx_sq * (type.test(has_fluid_left)
						+ type.test(has_fluid_right) )
						+ dy_sq * (type.test(has_fluid_bottom)
						+ type.test(has_fluid_top))
					);
				}

			}
		}
	}

	std::pair<Real, Real> update_velocity(grid::staggered_grid& grid, Real dt)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		
		auto dx = grid.get_dx();
		auto dy = grid.get_dy();
		
		std::pair<Real, Real> max_uv{0, 0};
		
		for (std::size_t j = 1; j < size_y - 1; ++j)
		{
			for (std::size_t i = 1; i < size_x - 1; ++i)
			{
				auto const& type = grid.cell_type(i, j);
				
				if (type[is_fluid])
				{
					if (type[has_fluid_right])
					{
						grid.u(i, j) = grid.f(i, j) - dt / dx *
							(grid.p(i + 1, j) - grid.p(i, j));

						max_uv.first = std::abs(grid.u(i, j)) > max_uv.first ? std::abs(grid.u(i, j)) : max_uv.first;
					}

					if (type[has_fluid_top])
					{
						grid.v(i, j) = grid.g(i, j) - dt / dy *
							(grid.p(i, j + 1) - grid.p(i, j));

						max_uv.second = std::abs(grid.v(i, j)) > max_uv.second ? std::abs(grid.v(i, j)) : max_uv.second;
					}
				}
			}
		}
		
		return max_uv;
	}

};

} //namespace time_integrator
} //namespace nast
#endif
