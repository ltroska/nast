#ifndef NAST_HPX_GRID_BOUNDARY_DATA_HPP_
#define NAST_HPX_GRID_BOUNDARY_DATA_HPP_

#include "util/defines.hpp"
#include "grid/boundary_data.hpp"

#include <iostream>
#include <string>

namespace nast { namespace grid {
	
namespace {
const std::string bc_to_string[4]  = {"noslip", "slip", "outstream", "instream"};
}
enum boundary_condition_type {
	noslip,
	slip,
	outstream,
	instream
};
	
/// This class represents one set of boundary conditions.
struct boundary_conditions
{
	using data_type = boundary_data<Real>;
	
    boundary_conditions()
    : left(0), right(0), bottom(0), top(0), external_forces(0),
      left_type(boundary_condition_type::noslip), right_type(boundary_condition_type::noslip),
      bottom_type(boundary_condition_type::noslip), top_type(boundary_condition_type::noslip)
    {}

    data_type left, right, bottom, top, external_forces;
    boundary_condition_type left_type, right_type, bottom_type, top_type;
    
    friend std::ostream& operator<<(std::ostream& os, boundary_conditions const& data);
};

std::ostream& operator<<(std::ostream& os, boundary_conditions const& data)
{
	os 		<< "Boundary conditions:"
			<< "\n\t external forces: "
			<< "{" << data.external_forces.x << ", " << data.external_forces.y << "},"
			<< "\n\t left: " << bc_to_string[data.left_type] << " {" << data.left.x << ", " << data.left.y << "}"
			<< "\n\t right: " << bc_to_string[data.right_type] << " {" << data.right.x << ", " << data.right.y << "}"
			<< "\n\t bottom: " << bc_to_string[data.bottom_type] << " {" << data.bottom.x << ", " << data.bottom.y << "}"
			<< "\n\t top: " << bc_to_string[data.top_type] << " {" << data.top.x << ", " << data.top.y << "}\n";
	return os;
}

}
}

#endif
