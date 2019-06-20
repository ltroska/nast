#ifndef NAST_GRID_BOUNDARY_DATA_HPP_
#define NAST_GRID_BOUNDARY_DATA_HPP_

namespace nast { namespace grid {

template<typename T>
struct boundary_data {

    boundary_data() : x(), y() {}

    boundary_data(T value) : x(value), y(value) {}

    boundary_data(T x_value, T y_value) : x(x_value), y(y_value) {}

    boundary_data(boundary_data<T> const& other) : x(other.x), y(other.y) {}

    boundary_data(boundary_data<T>&& other) : x(std::move(other.x)), y(std::move(other.y)) {}

    boundary_data<T>& operator=(boundary_data<T> const& other)
    {
        x = other.x;
        y = other.y;

        return *this;
    }

    boundary_data<T>& operator=(boundary_data<T>&& other)
    {
        x = std::move(other.x);
        y = std::move(other.y);

        return *this;
    }

    T x, y;

    friend std::ostream& operator<<(std::ostream& os, boundary_data<T> const& data)
    {
        os << "{" << data.x << "," << data.y << "}";
        return os;
    }

};

} //namespace grid
} //namespace nast

#endif
