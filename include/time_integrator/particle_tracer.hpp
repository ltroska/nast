#ifndef NAST_TIME_INTEGRATOR_PARTICLE_TRACER_HPP_
#define NAST_TIME_INTEGRATOR_PARTICLE_TRACER_HPP_

#include "grid/staggered_grid.hpp"

#include <vector>
#include <cmath>

namespace nast { namespace time_integrator {

class particle_tracer
{
public:
    static void advance_particles(std::vector<grid::particle>& particles, grid::staggered_grid const& grid, grid::boundary_conditions const& bcs, Real dt)
    {
        auto dx = grid.get_dx();
        auto dy = grid.get_dy();

        auto length_x = grid.get_length_x();
        auto length_y = grid.get_length_y();

        auto over_dx_dy = 1./(dx * dy);

        auto dx_over_2 = dx/2.;
        auto dy_over_2 = dy/2.;

        Real u, v;

        Real x_1, x_2, y_1, y_2;

        std::size_t i, j;

        for (auto& particle : particles)
        {
            i = std::floor(particle.x / dx) + 1;
            j = std::floor(particle.y/ dy) + 1;

            auto const& type = grid.cell_type(i, j);

            if(!type[is_fluid])
                continue;

            auto& x = particle.x;
            auto& y = particle.y;

            // u interpolation
            i = std::floor(particle.x / dx) + 1;
            j = std::floor( (particle.y + dy_over_2) / dy) + 1;

            x_1 = (i - 1) * dx;
            y_1 = ((j - 1) - 0.5) * dy;

            x_2 = i * dx;
            y_2 = (j - 0.5) * dy;

            auto& u_1 = grid.u(i - 1, j - 1);
            auto& u_2 = grid.u(i, j - 1);
            auto& u_3 = grid.u(i - 1, j);
            auto& u_4 = grid.u(i, j);

            u = over_dx_dy * ((x_2 - x)*(y_2 - y)*u_1 + (x - x_1)*(y_2 - y) * u_2 + (x_2 - x)*(y - y_1)*u_3 + (x - x_1)*(y - y_1)*u_4);


            // v interpolation
            i = std::floor((x + dx_over_2) / dx) + 1;
            j = std::floor(y / dy) + 1;

            x_1 = ((i - 1) - 0.5) * dx;
            y_1 = (j - 1) * dy;

            x_2 = (i - 0.5) * dx;
            y_2 = j * dy;

            auto& v_1 = grid.v(i - 1, j - 1);
            auto& v_2 = grid.v(i, j - 1);
            auto& v_3 = grid.v(i - 1, j);
            auto& v_4 = grid.v(i, j);

            v = over_dx_dy * ((x_2 - x)*(y_2 - y) * v_1 + (x - x_1) * (y_2 - y) * v_2 + (x_2 - x)*(y - y_1) * v_3 + (x - x_1)*(y - y_1)*v_4);

            x += dt * u;
            x = std::max(0., std::min(length_x, x));
            y += dt * v;
            y = std::max(0., std::min(length_y, y));
        }

        particles.erase(std::remove_if(particles.begin(), particles.end(),
                                    [&](const auto& particle)
                                    {
                                        auto& x = particle.x;
                                        auto& y = particle.y;
                                        return x == 0 && bcs.type[grid::direction::left] == grid::boundary_condition_type::outstream
                                            || x == length_x && bcs.type[grid::direction::right] == grid::boundary_condition_type::outstream
                                            || y == 0 && bcs.type[grid::direction::bottom] == grid::boundary_condition_type::outstream
                                            || y == length_y && bcs.type[grid::direction::top] == grid::boundary_condition_type::outstream;
                                    }),
                        particles.end());
    }


};

} // namespace time_integrator
} // namespace nast


#endif
