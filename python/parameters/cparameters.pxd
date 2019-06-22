# distutils: language=c++
# cython: language_level=3

from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "parameters/parameters.hpp" namespace "nast::parameters":
	enum solver_type:
		Jacobi = 0
		SOR = 1
	
	cdef struct parameters:
		double reynolds
	
		double initial_dt
		double t_end
		int max_timesteps
		
		double tau
		
		solver_type solver
		
		int max_solver_iterations
		double eps
		
		double alpha
		
		double omega
				
		bool verbose
		
		string to_string()
