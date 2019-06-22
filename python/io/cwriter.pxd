# distutils: language=c++
# cython: language_level=3

cimport grid.cgrid as cgrid
from libcpp.string cimport string

cdef extern from "io/vtk_writer.hpp" namespace "nast::io":	
	cdef cppclass vtk_writer:
		vtk_writer() except+
		
		void write_grid(string, cgrid.staggered_grid&)
