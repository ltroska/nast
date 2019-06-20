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
		
	bcs.left_type = nast::grid::boundary_condition_type::instream;
	bcs.left.x = 1;
	
	bcs.right_type = nast::grid::boundary_condition_type::outstream;

	bcs.top_type = nast::grid::boundary_condition_type::slip;
	bcs.bottom_type = nast::grid::boundary_condition_type::slip;
	

	sim.set_boundary_conditions(bcs);
	
	sim.set_grid_size(50, 20);
	sim.set_length(5, 1);
	
	std::cout << sim.params << std::endl;
	std::cout << sim.bcs << std::endl;
	
	sim.grid.toggle_cell_type(1, 20, 1, 9);
	sim.grid.toggle_cell_type(25, 35, 11, 18);
	
	sim.params.verbose = true;
	sim.params.t_end = 0;
	sim.params.max_timesteps = 0;
	
	sim.run();
	
    return 1;
}
