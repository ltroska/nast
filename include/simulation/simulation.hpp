#ifndef NAST_SIMULATION_SIMULATION_HPP_
#define NAST_SIMULATION_SIMULATION_HPP_

#include "grid/staggered_grid.hpp"
#include "grid/boundary_conditions.hpp"
#include "grid/boundary_data.hpp"
#include "time_integrator/time_integrator.hpp"
#include "time_integrator/particle_tracer.hpp"
#include "grid/particle.hpp"
#include "parameters/parameters.hpp"

#include "io/vtk_writer.hpp"

#include "solvers/sor.hpp"
#include "solvers/jacobi.hpp"

namespace nast { namespace simulation {

	class simulation
	{
	public:
        using index = nast::grid::staggered_grid::index;

		simulation() : grid(0, 0)
		{
		}

		void set_boundary_condition(grid::direction direction, grid::boundary_condition_type type, grid::boundary_data<Real> value = {0, 0})
        {
            bcs.type[direction] = type;
            bcs.value[direction] = value;
        }

        void set_boundary_condition(grid::direction direction, grid::boundary_condition_type type, Real u, Real v)
        {
            set_boundary_condition(direction, type, {u, v});
        }

		void set_boundary_conditions(grid::boundary_conditions bcs_)
		{
			bcs = bcs_;
		}

		void set_parameters(parameters::parameters parameters_)
		{
			parameters = parameters_;
        }

		void set_grid_size(std::size_t size_x, std::size_t size_y)
		{
			grid.resize(size_x, size_y);
		}

		void set_grid_length(Real length_x, Real length_y)
		{
			grid.set_length(length_x, length_y);
		}

		void set_obstacle(std::size_t i_min, std::size_t i_max, std::size_t j_min, std::size_t j_max)
        {
            grid.set_obstacle(i_min, i_max, j_min, j_max);
        }

        void set_fluid(std::size_t i_min, std::size_t i_max, std::size_t j_min, std::size_t j_max)
        {
            grid.set_fluid(i_min, i_max, j_min, j_max);
        }

        void toggle_cell_type(std::size_t i, std::size_t j, std::size_t width)
        {
            grid.toggle_cell_type(i, j, width);
        }

        void sanitize_cell_types()
        {
            grid.sanitize_cell_types();
        }

        void reset_cell_types()
        {
            grid.reset_cell_types();
        }

        std::vector<index> const& get_obstacle_cells() const
        {
            return grid.obstacle_cells;
        }

        void add_particles(grid::particle_distribution distribution, std::size_t amount)
        {
            particle_tracer.add_particles(particles, grid, bcs, distribution, amount);
        }

        void advance_particles(Real dt)
        {
            particle_tracer.advance_particles(particles, grid, bcs, dt);
        }

		Real do_timestep(Real dt)
		{
			if (dt < 1e-12)
			{
				dt = parameters.initial_dt;
			}
			return time_integrator.do_timestep(grid, bcs, parameters, dt);
		}

		void run()
		{
			Real t = 0;
			Real dt = parameters.initial_dt;
            Real old_dt;

			std::size_t timestep = 0;
			while (true)
			{
				++timestep;

                old_dt = dt;

				if (parameters.verbose)
				{
					std::cout << "Timestep " << timestep << ": time " << t << " dt " << dt << std::endl;
				}

				dt = do_timestep(dt);

                advance_particles(old_dt);

				write_grid_to_file("fields_" + std::to_string(timestep) + ".vtr");
                write_particles_to_file("particles_" + std::to_string(timestep) + ".vtu");

                if (timestep % 10 == 0)
                    add_particles(grid::particle_distribution::random, 100);


				t += dt;

				if ( (parameters.t_end > 0 && t > parameters.t_end) || (parameters.max_timesteps > 0 && timestep >= parameters.max_timesteps) )
				{
					break;
				}
			}
		}

		void write_grid_to_file(std::string filename)
		{
			writer.write_grid(filename, grid);
		}

        void write_particles_to_file(std::string filename)
		{
			writer.write_particles(filename, particles);
		}

		parameters::parameters parameters;

		grid::staggered_grid grid;
		grid::boundary_conditions bcs;
		time_integrator::time_integrator time_integrator;
        time_integrator::particle_tracer particle_tracer;

        std::vector<grid::particle> particles;

	private:
		io::vtk_writer writer;
	};

} // namespace simulation
} // namespace nast

#endif
