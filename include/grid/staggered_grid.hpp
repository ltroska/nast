#ifndef NAST_GRID_STAGGERED_GRID_HPP_
#define NAST_GRID_STAGGERED_GRID_HPP_

#include <bitset>

#include "grid/grid_data.hpp"

namespace nast { namespace grid {
	
	class staggered_grid
	{
	public:
		staggered_grid(std::size_t size_x_, std::size_t size_y_) 
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
				
		void toggle_cell_type(std::size_t i_min, std::size_t i_max, std::size_t j_min, std::size_t j_max)
		{
			for (std::size_t i = i_min; i <= i_max; ++i)
			{
				for (std::size_t j = j_min; j <= j_max; ++j)
				{
					auto& type = cell_type(i, j);
										
					type.flip(is_obstacle);
					type.flip(is_fluid);
					
					cell_type(i - 1, j).flip(has_fluid_right);
					cell_type(i + 1, j).flip(has_fluid_left);
					cell_type(i, j - 1).flip(has_fluid_top);
					cell_type(i, j + 1).flip(has_fluid_bottom);					
				}
			}
		}
				
		void reset_cell_types()
		{			
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
		
	private:
		std::size_t size_x, size_y, size;
		Real length_x, length_y;
		
	};
	
}//namespace grid
}//namespace nast


#endif
