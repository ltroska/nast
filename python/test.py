import nast

sim = nast.Simulation()

sim.parameters.verbose = 1
sim.parameters.max_timesteps = 1000
sim.parameters.t_end = 0

sim.set_boundary_condition(nast.Direction.left, nast.BoundaryConditionType.instream, (1, 0))
sim.set_boundary_condition(nast.Direction.right, nast.BoundaryConditionType.outstream)
sim.set_boundary_condition(nast.Direction.bottom, nast.BoundaryConditionType.slip)
sim.set_boundary_condition(nast.Direction.top, nast.BoundaryConditionType.slip)

sim.set_grid_size(50, 20)
sim.set_grid_length(5, 2)

sim.toggle_cell_type(10, 10)
sim.toggle_cell_type(9, 10)
#sim.toggle_cell_type(10, 9, 2)
#sim.toggle_cell_type(9, 9)

print(sim.grid.cell_type)

sim.sanitize_cell_types()

print()
print()
print(sim.grid.cell_type)

print(sim.boundary_conditions)


print(sim.parameters)

sim.set_obstacle(1, 20, 1, 9);
sim.set_obstacle(25, 35, 11, 18);

#sim.add_particles(nast.ParticleDistribution.uniform, 100)

#sim.run()

	

