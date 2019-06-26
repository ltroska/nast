# distutils: language=c++
# cython: language_level=3

from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair

cimport grid.cgrid as cgrid
cimport parameters.cparameters as cparameters
cimport timeintegrator.ctimeintegrator as ctimeintegrator


cdef extern from "simulation/simulation.hpp" namespace "nast::simulation":
    cdef cppclass simulation:
        simulation() except+

        void set_boundary_condition(cgrid.direction, cgrid.boundary_condition_type, double, double)

        void set_grid_size(double, double)
        void set_grid_length(double, double)

        void set_obstacle(int, int, int, int)
        void set_fluid(int, int, int, int)
        void toggle_cell_type(int, int, int)
        void reset_cell_types()
        void sanitize_cell_types()


        ctypedef pair[unsigned long, unsigned long] index
        vector[index]& get_obstacle_cells()

        void add_particles(cgrid.particle_distribution, int)
        void advance_particles(double)

        double do_timestep(double)
        void run()

        void write_grid_to_file(string)
        void write_particles_to_file(string)

        cgrid.boundary_conditions bcs
        cgrid.staggered_grid grid
        cparameters.parameters parameters
