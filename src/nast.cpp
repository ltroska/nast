#include <iostream>

#include "simulation/simulation.hpp"
#include "parameters/parameters.hpp"
#include "time_integrator/time_integrator.hpp"
#include "grid/boundary_conditions.hpp"
#include "solvers/sor.hpp"
#include "solvers/jacobi.hpp"

int main(int argc, char* argv[])
{
	nast::simulation::simulation sim;	
												
	nast::grid::boundary_conditions bcs;
		
	bcs.type[nast::grid::direction::left] = nast::grid::boundary_condition_type::instream;
	bcs.value[nast::grid::direction::left].x = 1;
	
	bcs.type[nast::grid::direction::right] = nast::grid::boundary_condition_type::outstream;

	bcs.type[nast::grid::direction::top] = nast::grid::boundary_condition_type::slip;
	bcs.type[nast::grid::direction::bottom] = nast::grid::boundary_condition_type::slip;
	

	sim.set_boundary_conditions(bcs);
	
	sim.set_grid_size(50, 20);
	sim.set_length(5, 1);
	
	std::cout << sim.params << std::endl;
	std::cout << sim.bcs << std::endl;
	
	sim.grid.set_obstacle(1, 20, 1, 9);
	sim.grid.set_obstacle(25, 35, 11, 18);
	
	sim.params.max_solver_iterations = 1000;
	sim.params.eps = 1e-10;
	sim.params.verbose = true;
	sim.params.t_end = 0;
	sim.params.max_timesteps = 10;
	
	sim.run();
	
    return 1;
}
