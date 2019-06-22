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

enum direction {
	left,
	right,
	bottom,
	top,
	external_forces
};
	
/// This class represents one set of boundary conditions.
struct boundary_conditions
{
	using data_type = boundary_data<Real>;

    boundary_conditions()
     
    {
		 type[direction::left] = boundary_condition_type::noslip;
		  type[direction::right]=boundary_condition_type::noslip;
      type[direction::bottom]=boundary_condition_type::noslip;
       type[direction::top]=boundary_condition_type::noslip;
       type[direction::external_forces]=boundary_condition_type::noslip;
       
		 value[direction::left] = 0;
		  value[direction::right]=0;
      value[direction::bottom]=0;
       value[direction::top]=0;
       value[direction::external_forces]=0;
       }

    data_type value[5];
    boundary_condition_type type[5];
    
    void set_type(direction dir, boundary_condition_type t)
    {
		type[dir] = t;
	}
	
    void set_value(direction dir, Real u, Real v)
    {
		value[dir].x = u;
		value[dir].y = v;		
	}
    
    friend std::ostream& operator<<(std::ostream& os, boundary_conditions const& data);
    
	std::string to_string() const
	{
		std::ostringstream stream;
		
		stream << *this;
		
		return stream.str();
	} 
};

std::ostream& operator<<(std::ostream& os, boundary_conditions const& data)
{
	os 		<< "Boundary conditions:"
			<< "\n\t external forces: "
			<< "{" << data.value[direction::external_forces].x << ", " << data.value[direction::external_forces].y << "},"
			<< "\n\t left: " << bc_to_string[data.type[direction::left]] << " {" << data.value[direction::left].x << ", " << data.value[direction::left].y << "}"
			<< "\n\t right: " << bc_to_string[data.type[direction::right]] << " {" << data.value[direction::right].x << ", " << data.value[direction::right].y << "}"
			<< "\n\t bottom: " << bc_to_string[data.type[direction::bottom]] << " {" << data.value[direction::bottom].x << ", " << data.value[direction::bottom].y << "}"
			<< "\n\t top: " << bc_to_string[data.type[direction::top]] << " {" << data.value[direction::top].x << ", " << data.value[direction::top].y << "}\n";
	return os;
}

}
}

#endif
