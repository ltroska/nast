#include <iostream>
#include <random>

#include "simulation/simulation.hpp"
#include "parameters/parameters.hpp"
#include "time_integrator/time_integrator.hpp"
#include "grid/boundary_conditions.hpp"
#include "solvers/sor.hpp"
#include "solvers/jacobi.hpp"
#include "time_integrator/particle_tracer.hpp"

int main(int, char**)
{
	nast::simulation::simulation sim;

    sim.set_boundary_condition(nast::grid::direction::left, nast::grid::boundary_condition_type::instream, {1, 0});
    sim.set_boundary_condition(nast::grid::direction::right, nast::grid::boundary_condition_type::outstream);
    sim.set_boundary_condition(nast::grid::direction::bottom, nast::grid::boundary_condition_type::slip);
    sim.set_boundary_condition(nast::grid::direction::top, nast::grid::boundary_condition_type::slip);

	sim.set_grid_size(100, 40);
	sim.set_grid_length(5, 1);

    sim.parameters.max_timesteps = 1000;
    sim.parameters.t_end = 0;
    sim.parameters.verbose = 1;

	std::cout << sim.parameters << std::endl;
	std::cout << sim.bcs << std::endl;

	sim.set_obstacle(1, 20, 1, 9);
	sim.set_obstacle(25, 35, 11, 18);

    sim.run();

    return 1;
}
