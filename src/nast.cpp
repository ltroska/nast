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
	
    std::vector<nast::grid::particle> particles;

    particles.reserve(1000);


    writer.write_particles("particles", particles);

    nast::time_integrator::time_integrator integrator;
    nast::time_integrator::particle_tracer tracer;

    Real dt = sim.params.initial_dt;

    Real t = 0;
    Real old_dt;

    std::size_t timestep = 0;

    while (true)
    {
        ++timestep;
        old_dt = dt;
        dt = integrator.do_timestep(sim.grid, bcs, sim.params, dt);

        if (timestep % 10 == 0)
            tracer.add_particles(particles, sim.grid, bcs, nast::time_integrator::particle_distribution::inflow, 100);

        tracer.advance_particles(particles, sim.grid, bcs, old_dt);

        t += dt;

        writer.write_grid("flow_" + std::to_string(timestep), sim.grid);
        writer.write_particles("particles_" + std::to_string(timestep), particles);

        std::cout << "Timestep " << timestep << std::endl;

        if (timestep == 1000)
            break;

    }
	
    return 1;
}
