# distutils: language=c++
# cython: language_level=3

cimport parameters.cparameters as cparameters
cimport grid.cgrid as cgrid
cimport timeintegrator.ctimeintegrator as ctimeintegrator
cimport io.cwriter as cwriter

from cython.operator import dereference as deref
from enum import Enum
from functools import partial

cdef class VTKWriter:
	cdef cwriter.vtk_writer _c_vtk_writer
	
	def __cinit__(self):
		self._c_vtk_writer = cwriter.vtk_writer()
		
	def write_grid(self, filename, Grid grid):
		self._c_vtk_writer.write_grid(filename.encode('utf-8'), deref(grid._this_ptr))


cdef class TimeIntegrator:
	cdef ctimeintegrator.time_integrator* _this_ptr
	
	def __cinit__(self):
		self._this_ptr = new ctimeintegrator.time_integrator()
		if self._this_ptr == NULL:
			raise MemoryError()
		
	def __dealloc__(self):
		if self._this_ptr != NULL:
			del self._this_ptr
		
	def do_timestep(self, Grid grid, BoundaryConditions boundary_conditions, Parameters parameters, dt = 0.):
		return deref(self._this_ptr).do_timestep(deref(grid._this_ptr), deref(boundary_conditions._this_ptr), parameters._c_parameters, dt)
	
class BoundaryConditionType(Enum):
	noslip = 0
	slip = 1
	outstream = 2
	instream = 3	
	
class Direction(Enum):
	left = 0
	right = 1
	bottom = 2
	top = 3	
	
class SolverType(Enum):
	jacobi = 0
	sor = 1

cdef class BoundaryConditions:
	cdef cgrid.boundary_conditions* _this_ptr
	
	def __cinit__(self):
		self._this_ptr = new cgrid.boundary_conditions()
		if self._this_ptr == NULL:
			raise MemoryError()
			
	def __dealloc__(self):
		if self._this_ptr != NULL:
			del self._this_ptr
	
	def __repr__(self):
		return deref(self._this_ptr).to_string().decode("UTF_8")
		
	def get_type(self, direction):
		return deref(self._this_ptr).type[direction.value]
		
	def set_type(self, direction, value):
		deref(self._this_ptr).set_type(direction.value, value.value)
		
	def get_value(self, direction):
		return (deref(self._this_ptr).value[direction.value].x, deref(self._this_ptr).value[direction.value].y)
		
	def set_value(self, direction, value):
		u, v = value
		deref(self._this_ptr).set_value(direction.value, u, v)
		
	# hack to avoid code duplication
	left_type = property(partial(lambda d, s: BoundaryConditions.get_type(s, d), Direction.left), partial(lambda d, s, v: BoundaryConditions.set_type(s, d, v), Direction.left))
	right_type = property(partial(lambda d, s: BoundaryConditions.get_type(s, d), Direction.right), partial(lambda d, s, v: BoundaryConditions.set_type(s, d, v), Direction.right))
	bottom_type = property(partial(lambda d, s: BoundaryConditions.get_type(s, d), Direction.bottom), partial(lambda d, s, v: BoundaryConditions.set_type(s, d, v), Direction.bottom))
	top_type = property(partial(lambda d, s: BoundaryConditions.get_type(s, d), Direction.top), partial(lambda d, s, v: BoundaryConditions.set_type(s, d, v), Direction.top))
	
	left = property(partial(lambda d, s: BoundaryConditions.get_value(s, d), Direction.left), partial(lambda d, s, v: BoundaryConditions.set_value(s, d, v), Direction.left))
	right = property(partial(lambda d, s: BoundaryConditions.get_value(s, d), Direction.right), partial(lambda d, s, v: BoundaryConditions.set_value(s, d, v), Direction.right))
	bottom = property(partial(lambda d, s: BoundaryConditions.get_value(s, d), Direction.bottom), partial(lambda d, s, v: BoundaryConditions.set_value(s, d, v), Direction.bottom))
	top = property(partial(lambda d, s: BoundaryConditions.get_value(s, d), Direction.top), partial(lambda d, s, v: BoundaryConditions.set_value(s, d, v), Direction.top))
				
cdef class Parameters:
	cdef cparameters.parameters _c_parameters
	
	
	def __cinit__(self):
		self._c_parameters = cparameters.parameters()

	def __repr__(self):
		return self._c_parameters.to_string().decode("UTF_8")

	def get_reynolds(self):
		return self._c_parameters.reynolds
	
	def set_reynolds(self, value):
		self._c_parameters.reynolds = value	
		
	reynolds = property(get_reynolds, set_reynolds)

	def get_initial_dt(self):
		return self._c_parameters.initial_dt
	
	def set_initial_dt(self, value):
		self._c_parameters.initial_dt = value	
		
	initial_dt = property(get_initial_dt, set_initial_dt)

	def get_t_end(self):
		return self._c_parameters.t_end
	
	def set_t_end(self, value):
		self._c_parameters.t_end = value	
		
	t_end = property(get_t_end, set_t_end)

	def get_max_timesteps(self):
		return self._c_parameters.max_timesteps
	
	def set_max_timesteps(self, value):
		self._c_parameters.max_timesteps = value	
		
	max_timesteps = property(get_max_timesteps, set_max_timesteps)

	def get_tau(self):
		return self._c_parameters.tau
	
	def set_tau(self, value):
		self._c_parameters.tau = value	
		
	tau = property(get_tau, set_tau)

	def get_solver(self):
		return self._c_parameters.solver
	
	def set_solver(self, value):
		self._c_parameters.solver = value.value	
		
	solver = property(get_solver, set_solver)

	def get_max_solver_iterations(self):
		return self._c_parameters.max_solver_iterations
	
	def set_max_solver_iterations(self, value):
		self._c_parameters.max_solver_iterations = value	
		
	max_solver_iterations = property(get_max_solver_iterations, set_max_solver_iterations)

	def get_eps(self):
		return self._c_parameters.eps
	
	def set_eps(self, value):
		self._c_parameters.eps = value	
		
	eps = property(get_eps, set_eps)

	def get_alpha(self):
		return self._c_parameters.alpha
	
	def set_alpha(self, value):
		self._c_parameters.alpha = value	
		
	alpha = property(get_alpha, set_alpha)

	def get_omega(self):
		return self._c_parameters.omega
	
	def set_omega(self, value):
		self._c_parameters.omega = value	
		
	omega = property(get_omega, set_omega)

	def get_verbose(self):
		return self._c_parameters.verbose
	
	def set_verbose(self, value):
		self._c_parameters.verbose = value	
		
	verbose = property(get_verbose, set_verbose)
		
				
cdef class GridDoubleDataView:	
	cdef cgrid.grid_data[double]* _c_grid_data		
		
	def __str__(self):
		return deref(self._c_grid_data).to_string().decode("UTF_8")
	
	def __getitem__(self, idx):
		i, j = idx
		return deref(self._c_grid_data)(i, j)
		
	def __setitem__(self, idx, value):
		i, j = idx
		deref(self._c_grid_data).assign(i, j, value)
		
cdef class Bitset:
	cdef cgrid.bitset* _c_bitset
		
	def count(self):
		return self._c_bitset.count()
		
	def size(self):
		return self._c_bitset.size()
		
	def test(self, pos):
		return self._c_bitset.test(pos)
		
	def any(self):
		return self._c_bitset.any()
		
	def none(self):
		return self._c_bitset.none()
		
	def all(self):
		return self._c_bitset.all()
		
	def set(self):
		self._c_bitset.set()
		
	def set(self, pos, value):
		self._c_bitset.set(pos, value)
		
	def reset(self):
		self._c_bitset.reset()
		
	def reset(self, pos):
		self._c_bitset.reset(pos)
		
	def flip(self):
		self._c_bitset.flip()
		
	def flip(self, pos):
		self._c_biset.flip(pos)
		
	def __str__(self):
		return deref(self._c_bitset).to_string().decode("UTF_8")
			
		
cdef class GridBitsetDataView:	
	cdef cgrid.grid_data[cgrid.bitset]* _c_grid_data			
	def __str__(self):
		return deref(self._c_grid_data).to_string().decode("UTF_8")
	
	def __getitem__(self, idx):
		i, j = idx
		
		bitset = Bitset()
		bitset._c_bitset = &deref(self._c_grid_data)(i, j)
		
		return bitset
		
		
cdef class Grid:
	cdef cgrid.staggered_grid* _this_ptr
	
	def __cinit__(self, size_x = 0, size_y = 0):
		self._this_ptr = new cgrid.staggered_grid(size_x, size_y)
		if self._this_ptr == NULL:
			raise MemoryError()
			
	def __dealloc__(self):
		if self._this_ptr != NULL:
			del self._this_ptr
		
	def resize(self, size_x, size_y):
		deref(self._this_ptr).resize(size_x, size_y)
		
	def set_length(self, length_x, length_y):
		deref(self._this_ptr).set_length(length_x, length_y)
		
	def reset_cell_types(self):
		deref(self._this_ptr).reset_cell_types()
		
	def set_fluid(self, i_min, i_max, j_min, j_max):
		deref(self._this_ptr).set_fluid(i_min, i_max, j_min, j_max)
		
	def set_obstacle(self, i_min, i_max, j_min, j_max):
		deref(self._this_ptr).set_obstacle(i_min, i_max, j_min, j_max)
		
	def get_size_x(self):
		return deref(self._this_ptr).get_size_x()
		
	def get_size_y(self):
		return deref(self._this_ptr).get_size_y()
		
	def get_size(self):
		return deref(self._this_ptr).get_size()
		
	def get_length_x(self):
		return deref(self._this_ptr).get_length_x()
		
	def get_length_y(self):
		return deref(self._this_ptr).get_length_y()
		
	def get_dx(self):
		return deref(self._this_ptr).get_dx()
		
	def get_dy(self):
		return deref(self._this_ptr).get_dy()
		
	def get_dx_sq(self):
		return deref(self._this_ptr).get_dx_sq()
		
	def get_dy_sq(self):
		return deref(self._this_ptr).get_dy_sq()
		
	@property
	def u(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).u
		return data
		
	@property
	def v(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).v
		return data
		
	@property
	def f(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).f
		return data
		
	@property
	def g(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).g
		return data
		
	@property
	def rhs(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).rhs
		return data
		
	@property
	def p(self):
		data = GridDoubleDataView()
		data._c_grid_data = &deref(self._this_ptr).p
		return data
		
	@property
	def cell_type(self):
		data = GridBitsetDataView()
		data._c_grid_data = &deref(self._this_ptr).cell_type
		return data
	
	
