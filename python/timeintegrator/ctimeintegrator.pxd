# distutils: language=c++
# cython: language_level=3

cimport grid.cgrid as cgrid
cimport parameters.cparameters as cparameters


cdef extern from "time_integrator/time_integrator.hpp" namespace "nast::time_integrator":	
	cdef cppclass time_integrator:
		time_integrator() except+
		
		double do_timestep(cgrid.staggered_grid&, const cgrid.boundary_conditions&, const cparameters.parameters&, double)
