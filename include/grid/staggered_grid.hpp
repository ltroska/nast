#ifndef NAST_GRID_STAGGERED_GRID_HPP_
#define NAST_GRID_STAGGERED_GRID_HPP_

#include "grid/grid_data.hpp"

#include <bitset>
#include <algorithm>

namespace nast { namespace grid {
	
	class staggered_grid
	{
	public:
		staggered_grid(std::size_t size_x_ = 0, std::size_t size_y_ = 0) 
			: 
			  u(size_x_, size_y_), v(size_x_, size_y_), f(size_x_, size_y_),
			  g(size_x_, size_y_), p(size_x_, size_y_), rhs(size_x_, size_y_),
			  cell_type(size_x_, size_y_),
			  size_x(size_x_), size_y(size_y_), size(size_x_ * size_y_),
			  length_x(1), length_y(1)
		{
			reset_cell_types();
		}
		
		void set_length(Real length_x_, Real length_y_)
		{
			length_x = length_x_;
			length_y = length_y_;
		}
		
		void resize(std::size_t size_x_, std::size_t size_y_)
		{
			size_x = size_x_;
			size_y = size_y_;
			
			u.resize(size_x, size_y);
			v.resize(size_x, size_y);
			f.resize(size_x, size_y);
			g.resize(size_x, size_y);
			p.resize(size_x, size_y);
			rhs.resize(size_x, size_y);
			cell_type.resize(size_x, size_y);
			
			reset_cell_types();
		}
		
		void set_obstacle(std::size_t i_min, std::size_t i_max, std::size_t j_min, std::size_t j_max)
		{
			obstacle_cells.reserve(obstacle_cells.size() + (i_max - i_min + 1) * (j_max - j_min + 1));
			
			for (std::size_t i = i_min; i <= i_max; ++i)
			{
				for (std::size_t j = j_min; j <= j_max; ++j)
				{
					auto& type = cell_type(i, j);
										
					if ( (i == i_min && type[has_fluid_left]) 
						|| (i == i_max && type[has_fluid_right])
						|| (j == j_min && type[has_fluid_bottom])
						|| (j == j_max && type[has_fluid_top]) )
					{
						obstacle_cells.emplace_back(i, j);
					}
					
					fluid_cells.erase(std::remove(fluid_cells.begin(), fluid_cells.end(), std::make_pair(i, j)), fluid_cells.end());
					
										
					type.set(is_obstacle);
					type.set(is_fluid, 0);
					
					cell_type(i - 1, j).set(has_fluid_right, 0);
					cell_type(i + 1, j).set(has_fluid_left, 0);
					cell_type(i, j - 1).set(has_fluid_top, 0);
					cell_type(i, j + 1).set(has_fluid_bottom, 0);					
				}
			}
		}
		
		void set_fluid(std::size_t i_min, std::size_t i_max, std::size_t j_min, std::size_t j_max)
		{
			fluid_cells.reserve(fluid_cells.size() + (i_max - i_min + 1) * (j_max - j_min + 1));
			
			for (std::size_t i = i_min; i <= i_max; ++i)
			{
				for (std::size_t j = j_min; j <= j_max; ++j)
				{
					auto& type = cell_type(i, j);
										
					if ( (i == i_min && type[has_fluid_left]) 
						|| (i == i_max && type[has_fluid_right])
						|| (j == j_min && type[has_fluid_bottom])
						|| (j == j_max && type[has_fluid_top]) )
					{
						obstacle_cells.erase(std::remove(obstacle_cells.begin(), obstacle_cells.end(), std::make_pair(i, j)), obstacle_cells.end());
					}
					
					fluid_cells.emplace_back(i, j);
					
										
					type.set(is_fluid);
					type.set(is_obstacle, 0);
					
					cell_type(i - 1, j).set(has_fluid_right, 1);
					cell_type(i + 1, j).set(has_fluid_left, 1);
					cell_type(i, j - 1).set(has_fluid_top, 1);
					cell_type(i, j + 1).set(has_fluid_bottom, 1);					
				}
			}
		}
								
		void reset_cell_types()
		{			
			obstacle_cells.clear();
			fluid_cells.clear();
			
			fluid_cells.reserve( (size_x - 2) * (size_y - 2) );
			
			for (std::size_t i = 0; i < size_x; ++i)
			{
				for (std::size_t j = 0; j < size_y; ++j)
				{
					auto& type = cell_type(i, j);
					
					type.reset();
					
					if (i == 0 || i == size_x - 1 || j == 0 || j == size_y - 1)
					{
						type.set(is_boundary);
						
						u(i, j) = 0;
						v(i, j) = 0;
						f(i, j) = 0;
						g(i, j) = 0;
						p(i, j) = 0;
						
						if (i == 0 && j != 0 && j != size_y - 1) type.set(has_fluid_right);
						if (i == size_x -1 && j != 0 && j != size_y - 1) type.set(has_fluid_left);
						if (j == 0 && i != 0 && i != size_x - 1) type.set(has_fluid_top);
						if (j == size_y - 1 && i != 0 && i != size_x - 1) type.set(has_fluid_bottom);
					}
					else
					{
						type.set(is_fluid);
						
						if (i - 1  != 0) type.set(has_fluid_left);
						if (i + 1 != size_x -1) type.set(has_fluid_right);
						if (j - 1 != 0) type.set(has_fluid_bottom);
						if (j + 1 != size_y - 1) type.set(has_fluid_top);
						
						fluid_cells.emplace_back(i, j);
					}					
				}
			}
		}
		
		std::size_t get_size_x() const { return size_x; }
		std::size_t get_size_y() const { return size_y; }
		std::size_t get_size() const { return size; }
		
		Real get_length_x() const { return length_x; }
		Real get_length_y() const { return length_y; }
		
		Real get_dx() const { return length_x/(size_x - 2); }
		Real get_dy() const { return length_y/(size_y - 2); }
		
		Real get_dx_sq() const { auto dx = get_dx(); return dx * dx; }
		Real get_dy_sq() const { auto dy = get_dy(); return dy * dy; }
		
		grid_data<Real> u;
		grid_data<Real> v;
		grid_data<Real> f;
		grid_data<Real> g;
		grid_data<Real> p;
		grid_data<Real> rhs;
		
		grid_data<std::bitset<NUM_BITS> > cell_type;
		
		std::vector<std::pair<std::size_t, std::size_t> > fluid_cells;
		std::vector<std::pair<std::size_t, std::size_t> > obstacle_cells;
		
	private:
		std::size_t size_x, size_y, size;
		Real length_x, length_y;
		
	};
	
}//namespace grid
}//namespace nast


#endif
