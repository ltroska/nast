import nast

bcs = nast.BoundaryConditions()

bcs.left_type = nast.BoundaryConditionType.instream
bcs.left = (1, 0)

bcs.right_type = nast.BoundaryConditionType.outstream

bcs.bottom_type = nast.BoundaryConditionType.slip
bcs.top_type = nast.BoundaryConditionType.slip	
											
grid = nast.Grid(50, 20)
grid.set_length(5, 1)
	
params = nast.Parameters()

params.verbose = True
params.solver = nast.SolverType.sor
	
grid.set_obstacle(1, 20, 1, 9);
grid.set_obstacle(25, 35, 11, 18);

integrator = nast.TimeIntegrator()

writer = nast.VTKWriter()
		
t = 0
dt = params.initial_dt

timestep = 0
while True:
	timestep += 1

	if params.verbose:
		print("Timestep {}: time {} dt {}".format(timestep, t, dt))
		
	dt = integrator.do_timestep(grid, bcs, params, dt)
	
	writer.write_grid("test_{}.vtr".format(timestep), grid)
	
	t += dt
	
	if (timestep >= 10):
		break
	#if (params.t_end > 0 and t > params.t_end) or (params.max_timesteps > 0 and timestep >= params.max_timesteps):
	#	break

	

