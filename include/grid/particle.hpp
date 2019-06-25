#ifndef NAST_GRID_PARTICLE_HPP_
#define NAST_GRID_PARTICLE_HPP_

#include "util/defines.hpp"

namespace nast { namespace grid {

struct particle {
    particle(Real x_, Real y_, Real angle_) : x(x_), y(y_), angle(angle_) {}

    Real x, y, angle;
};


} // namespace nast
} // namespace grid

#endif
