#ifndef NAST_SIMULATION_SIMULATION_HPP_
#define NAST_SIMULATION_SIMULATION_HPP_

#include "grid/staggered_grid.hpp"
#include "grid/boundary_conditions.hpp"
#include "time_integrator/time_integrator.hpp"
#include "parameters/parameters.hpp"

#include "io/vtk_writer.hpp"

#include "solvers/sor.hpp"
#include "solvers/jacobi.hpp"

namespace nast { namespace simulation {
	
	class simulation
	{
	public:
		simulation() : grid(0, 0)
		{
			set_solver(params.solver);
		}
		
		void set_solver(parameters::solver_type type)
		{
			if (type == parameters::solver_type::Jacobi)
			{
				solver = std::make_shared<solvers::jacobi>();
			}
			else
			{
				solver = std::make_shared<solvers::sor>();
			}
		}
		
		void set_boundary_conditions(grid::boundary_conditions bcs_)
		{
			bcs = bcs_;
		}
		
		void set_parameters(parameters::parameters params_)
		{
			params = params;
			
			set_solver(params.solver);
		}
		
		void set_grid_size(std::size_t size_x, std::size_t size_y)
		{
			grid.resize(size_x, size_y);
		}
		
		void set_length(Real length_x, Real length_y)
		{
			grid.set_length(length_x, length_y);
		}
		
		Real do_timestep(Real dt)
		{					
			if (dt < 1e-12)
			{
				dt = params.initial_dt;
			}
			return integrator.do_timestep(grid, bcs, params, dt);	
		}
		
		void run()
		{
			Real t = 0;
			Real dt = 0;
			
			std::size_t timestep = 0;
			while (true)
			{
				++timestep;
				
				if (timestep == 1)
				{
					dt = params.initial_dt;
				}
					
				if (params.verbose)
				{
					std::cout << "Timestep " << timestep << ": time " << t << " dt " << dt << std::endl;
				}
				
				dt = do_timestep(dt);
			//	write_grid_to_file("test_" + std::to_string(timestep) + ".vtr");
				
				t += dt;
				
				if ( (params.t_end > 0 && t > params.t_end) || (params.max_timesteps > 0 && timestep >= params.max_timesteps) )
				{
					break;
				}
			}					
		}
		
		void write_grid_to_file(std::string filename)
		{
			writer.write_grid(filename, grid);
		}
		
		parameters::parameters params;

		grid::staggered_grid grid;
		grid::boundary_conditions bcs;
		time_integrator::time_integrator integrator;
		std::shared_ptr<solvers::solver_base> solver;
		
	private:
		io::vtk_writer writer;
	};
	
} // namespace simulation
} // namespace nast

#endif
