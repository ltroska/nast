# distutils: language=c++
# cython: language_level=3

from libcpp.string cimport string
from std.cbitset cimport bitset as bitset

cdef extern from "grid/boundary_data.hpp" namespace "nast::grid":

    cdef cppclass boundary_data "nast::grid::boundary_data<double>":
        boundary_data() except+
        boundary_data(double)
        boundary_data(double, double)

        string to_string()

        double x
        double y

cdef extern from "grid/boundary_conditions.hpp" namespace "nast::grid":
    cdef enum boundary_condition_type:
            noslip = 0
            slip = 1
            outstream = 2
            instream = 3

    cdef enum direction:
            left = 0
            right = 1
            bottom = 2
            top = 3

    cdef cppclass boundary_conditions:
        boundary_conditions() except+

        string to_string()

        void set_type(direction, boundary_condition_type)
        void set_value(direction, double, double)

        boundary_data[5] value
        boundary_condition_type type[5]

cdef extern from "grid/grid_data.hpp" namespace "nast::grid":

    cdef cppclass grid_data[T]:

        grid_data() except +

        T& operator()(int, int)

        void assign(int, int, T)

        string to_string()

cdef extern from "grid/particle.hpp" namespace "nast::grid":
    cdef enum particle_distribution:
        uniform = 0
        random = 1
        at_instream = 2

    cdef cppclass particle:
        double x
        double y
        double angle

cdef extern from "grid/staggered_grid.hpp" namespace "nast::grid":

    cdef cppclass staggered_grid:

        staggered_grid() except +
        staggered_grid(int, int) except +

        void set_length(double, double)

        void resize(int, int)

        void set_fluid(int, int, int, int)
        void set_obstacle(int, int, int, int)

        void reset_cell_types()

        int get_size_x()
        int get_size_y()
        int get_size()

        double get_length_x()
        double get_length_y()

        double get_dx()
        double get_dy()

        double get_dx_sq()
        double get_dy_sq()

        grid_data[double] u
        grid_data[double] v
        grid_data[double] f
        grid_data[double] g
        grid_data[double] p
        grid_data[double] rhs

        grid_data[bitset] cell_type

        #grid_data[bitset[7]] cell_type


