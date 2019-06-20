/** All stencils are provided in this header as inline methods
 */
#ifndef NAST_UTIL_DERIVATIVES_HPP_
#define NAST__UTIL_DERIVATIVES_HPP_

#include "grid/grid_data.hpp"

#include <cmath>
#include <cstdlib>

namespace nast { namespace util { namespace derivatives {

typedef grid::grid_data<double> data_type;

inline double second_derivative_fwd_bkwd_x(
    data_type const& grid, std::size_t i, std::size_t j, double dx_sq)
{
    return (grid(i + 1, j) - 2 * grid(i, j) + grid(i - 1, j)) / dx_sq;
}

inline double second_derivative_fwd_bkwd_y(
    data_type const& grid, std::size_t i, std::size_t j, double dy_sq)
{
    return (grid(i, j + 1) - 2 * grid(i, j) + grid(i, j - 1)) / dy_sq;
}

inline double first_derivative_fwd_x(
    data_type const& grid, std::size_t i, std::size_t j, double dx)
{
    return (grid(i + 1, j) - grid(i, j)) / dx;
}

inline double first_derivative_fwd_y(
    data_type const& grid, std::size_t i, std::size_t j, double dy)
{
    return (grid(i, j + 1) - grid(i, j)) / dy;
}

inline double first_derivative_bkwd_x(
    data_type const& grid, std::size_t i, std::size_t j, double dx)
{
    return (grid(i, j) - grid(i - 1, j)) / dx;
}

inline double first_derivative_bkwd_y(
    data_type const& grid, std::size_t i, std::size_t j, double dy)
{
    return (grid(i, j) - grid(i, j - 1)) / dy;
}

inline double first_derivative_of_square_x(
    data_type const& grid, std::size_t i, std::size_t j, double dx,
    double alpha = 0.9)
{
    return 1./dx * (std::pow((grid(i, j) + grid(i + 1, j)) / 2., 2)
                    - std::pow((grid(i - 1, j) + grid(i, j)) / 2., 2))
            + alpha / dx
            * (std::abs(grid(i, j) + grid(i + 1, j)) * (grid(i, j) - grid(i + 1, j)) / 4.
                - std::abs(grid(i - 1, j) + grid(i, j)) * (grid(i - 1, j) - grid(i, j)) / 4.);
}

inline double first_derivative_of_square_y(
    data_type const& grid, std::size_t i, std::size_t j, double dy,
    double alpha = 0.9)
{
    return 1./dy * (std::pow((grid(i, j) + grid(i, j + 1)) / 2., 2)
                        - std::pow((grid(i, j - 1) + grid(i, j)) / 2., 2))
            + alpha / dy
            * (std::abs(grid(i, j) + grid(i, j + 1)) * (grid(i, j) - grid(i, j + 1)) / 4.
                - std::abs(grid(i, j - 1) + grid(i, j)) * (grid(i, j - 1) - grid(i, j)) / 4.);
}

inline double first_derivative_of_uv_x(
    data_type const& u, data_type const& v, std::size_t i, std::size_t j, double dx,
    double alpha = 0.9)
{
    return 1./dx * ((u(i, j) + u(i, j + 1))
                        * (v(i, j) + v(i + 1, j)) / 4.
                        - (u(i - 1, j) + u(i - 1, j + 1)) * (v(i - 1, j ) + v(i, j)) / 4.)
            + alpha / dx * (std::abs(u(i, j) + u(i, j + 1))
                * (v(i, j) - v(i + 1, j)) / 4.
            - std::abs(u(i - 1, j) + u(i - 1, j + 1)) * (v(i - 1, j) - v(i, j)) / 4.);
}

inline double first_derivative_of_uv_y(
    data_type const& u, data_type const& v, std::size_t i, std::size_t j, double dy
    , double alpha = 0.9)
{
    return 1./dy * ((v(i, j) + v(i + 1, j))  * (u(i, j) + u(i, j + 1)) / 4.
            - (v(i, j - 1) + v(i + 1, j - 1)) * (u(i, j - 1) + u(i, j)) / 4.)
            + alpha / dy * (std::abs(v(i, j) + v(i + 1, j))
                * (u(i, j) - u(i, j + 1)) / 4.
            - std::abs(v(i, j - 1) + v(i + 1, j - 1)) * (u(i, j - 1) - u(i, j)) / 4.);
}

}
}
}
#endif
