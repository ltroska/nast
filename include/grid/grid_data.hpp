#ifndef NAST_GRID_GRID_DATA_HPP_
#define NAST_GRID_GRID_DATA_HPP_

#include "util/defines.hpp"

#include <vector>
#include <iostream>
#include <sstream>

namespace nast { namespace grid {

/// This class represents one variable of a grid.
template<typename T = Real>
struct grid_data
{
public:

    grid_data()
    : size_x_(0),
      size_y_(0),
      size_(0)
    {}

    grid_data(std::size_t size_x, std::size_t size_y, T val = T())
    : data_(size_x * size_y, val),
      size_x_(size_x),
      size_y_(size_y),
      size_(size_x * size_y)
    {}

    void resize(std::size_t size_x, std::size_t size_y, T val = T())
    {
        size_x_ = size_x;
        size_y_ = size_y;
        size_ = size_x * size_y;
        data_.resize(size_x_ * size_y_, val);
    }

    void clear(T val = T())
    {
        for (auto& elem : data_)
            elem = val;
    }

    inline T operator[](std::size_t idx) const 
    { 
		return data_[idx];
	}
		
    inline T& operator[](std::size_t idx)
    { 
		return data_[idx];
	}

    inline T& operator()(std::size_t idx, std::size_t idy)
    {
		return data_[to_flat_index(idx, idy)];
	}

    inline T const& operator()(std::size_t idx, std::size_t idy) const
    {
		return data_[to_flat_index(idx, idy)];
	}
	
	void assign(std::size_t idx, std::size_t idy, T value)
	{
		(*this)(idx, idy) = value;
	}

    typename std::vector<T>::iterator begin()
    { 
		return data_.begin();
	}
	
    typename std::vector<T>::iterator end()
    { 
		return data_.end(); 
	}

    std::vector<T> data_;
    std::size_t size_x_;
    std::size_t size_y_;
    std::size_t size_;
    
    template<typename U>
    friend std::ostream& operator<<(std::ostream &os, const grid_data<U>& data);
    
	std::string to_string() const
	{
		std::ostringstream stream;
		
		stream << *this;
		
		return stream.str();
	} 
    
protected:
	std::size_t to_flat_index(std::size_t idx, std::size_t idy) const
	{
		return idx + idy * size_x_;
	}
   
};

	
template<typename T>
std::ostream& operator<<(std::ostream &os, const grid_data<T>& data)
{
	if (data.size_y_ == 0)
	{
		os << "empty\n";
	}
	else
	{	
		for (std::size_t j = data.size_y_ - 1; j >= 0; --j)
		{
			for (std::size_t i = 0; i < data.size_x_; ++i)
			{
				os << data(i, j) << " ";
			}		
			if (j == 0)
			{
				break;
			}
			
			os << "\n";
		}
	}	
	return os;
}

}//namespace grid
}//namespace nast

#endif
