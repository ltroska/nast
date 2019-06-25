#include <iostream>
#include <random>

#include "simulation/simulation.hpp"
#include "parameters/parameters.hpp"
#include "time_integrator/time_integrator.hpp"
#include "grid/boundary_conditions.hpp"
#include "solvers/sor.hpp"
#include "solvers/jacobi.hpp"
#include "time_integrator/particle_tracer.hpp"

int main(int argc, char* argv[])
{
	nast::simulation::simulation sim;	
												
	nast::grid::boundary_conditions bcs;

    nast::io::vtk_writer writer;
		
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
	
    std::vector<nast::grid::particle> particles(10000);

    std::random_device rd;

    std::mt19937 e2(rd());

    std::uniform_real_distribution<> dist_x(0, 5);
    std::uniform_real_distribution<> dist_y(0, 1);

    for (auto& p : particles)
    {
        p.x = dist_x(rd);
        p.y = dist_y(rd);
    }


    writer.write_particles("particles", particles);

    nast::time_integrator::time_integrator integrator;
    nast::time_integrator::particle_tracer tracer;

    Real dt = sim.params.initial_dt;

    Real t = 0;
    Real old_dt;

    std::size_t iteration = 0;

    while (true)
    {
        ++iteration;
        old_dt = dt;
        dt = integrator.do_timestep(sim.grid, bcs, sim.params, dt);

        tracer.advance_particles(particles, sim.grid, bcs, old_dt);

        t += dt;

        writer.write_grid("flow_" + std::to_string(iteration), sim.grid);
        writer.write_particles("particles_" + std::to_string(iteration), particles);

    }

	sim.params.max_solver_iterations = 1000;
	sim.params.eps = 1e-10;
	sim.params.verbose = true;
	sim.params.t_end = 0;
	sim.params.max_timesteps = 10;
	
	sim.run();
	
    return 1;
}
