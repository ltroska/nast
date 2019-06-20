#ifndef NAST_IO_VTK_WRITER_HPP_
#define NAST_IO_VTK_WRITER_HPP_

#include "grid/staggered_grid.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>

namespace nast { namespace io {

namespace {
bool ends_with(std::string str, std::string suffix)
{
  if (str.length() < suffix.length())
    return false;

  return str.substr(str.length() - suffix.length()) == suffix;
}
} //anonymous namespace

struct vtk_writer
{
	void write_grid(std::string filename, grid::staggered_grid grid)
	{
		auto size_x = grid.get_size_x();
		auto size_y = grid.get_size_y();
		auto cells_x = size_x - 2;
		auto cells_y = size_y - 2;

		Real dx = grid.get_dx();
		Real dy = grid.get_dy();

		std::vector<double> strom(cells_x + 2, 0);

		if (!ends_with(filename, ".vtr")) 
		{
			filename.append (".vtr");
		}

		std::filebuf fb;
		fb.open (const_cast < char *>(filename.c_str ()), std::ios::out);
		std::ostream os (&fb);

		std::string coordinate_x;
		std::string coordinate_y;
		std::string coordinate_z;

		for (std::size_t x = 0; x <= cells_x; ++x)
			coordinate_x += std::to_string(dx * x) + " ";

		coordinate_x += "\n";

		for (std::size_t y = 0; y <= cells_y; ++y)
			coordinate_y += std::to_string(dy * y) + " ";

		coordinate_y += "\n";


		coordinate_z += "0\n";
				
		std::stringstream p_stream;
		p_stream << std::setprecision(std::numeric_limits<double>::digits10);

		std::stringstream obstacle_stream;

		std::stringstream uv_stream;
		uv_stream << std::setprecision(std::numeric_limits<double>::digits10);

		std::stringstream vorticity_stream;
		vorticity_stream
			<< std::setprecision(std::numeric_limits<double>::digits10);

		std::stringstream strom_stream;
		strom_stream << std::setprecision(std::numeric_limits<double>::digits10);

		std::stringstream heat_stream;
		heat_stream << std::setprecision(std::numeric_limits<double>::digits10);

		std::stringstream temp_stream;
		temp_stream << std::setprecision(std::numeric_limits<double>::digits10);

		for (std::size_t j = 1; j < size_y - 1; ++j) 
		{
			for (std::size_t i = 1; i < size_x - 1; ++i)
			{
				//uv_stream << i << " " << j << "\n";
				
				auto type = grid.cell_type(i, j);

				if (type[is_fluid])
				{
					obstacle_stream << "0\n";
					
					p_stream << grid.p(i, j)  << "\n";

					uv_stream << (grid.u(i, j) + grid.u(i - 1, j)) / 2. << " ";

					uv_stream << (grid.v(i, j) + grid.v(i, j - 1)) / 2. << " 0\n";
				}
				else 
				{
					obstacle_stream << "1\n";
					
					p_stream << "0\n";
					
					uv_stream << "0 0 0\n";
				}
			}
		}
				
		os  << std::setprecision(std::numeric_limits<double>::digits10)
			<< "<?xml version=\"1.0\"?>" << std::endl
			<< "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
				<< std::endl
			<< "<RectilinearGrid WholeExtent=\""
				<< 0 << " " << cells_x << " " << 0 << " " << cells_y << " " << 0 << " " << 0 << "\">" << std::endl
			<< "<Piece Extent=\""
				<< 0 << " " << cells_x << " " << 0 << " " << cells_y << " " << 0 << " " << 0 << "\">" << std::endl
			
			<< "<PointData>" << std::endl
			<< "</PointData>" << std::endl
			
			<< "<CellData>" << std::endl
			<< "<DataArray type=\"Float32\" Name=\"pressure\">" << std::endl
			<< p_stream.str() << std::endl
			<< "</DataArray>" << std::endl
			
			<< "<DataArray type=\"Int32\" Name=\"obstacle\">" << std::endl
			<< obstacle_stream.str() << std::endl
			<< "</DataArray>" << std::endl
			
			<< "<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\">" << std::endl
			<< uv_stream.str() << std::endl
			<< "</DataArray>" << std::endl
			<< "</CellData>" << std::endl
			
			<< "<Coordinates>" << std::endl			
			<< "<DataArray type=\"Float32\" Name=\"X_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
				<< std::endl
			<< coordinate_x << std::endl
			<< "</DataArray>" << std::endl
			<< "<DataArray type=\"Float32\" Name=\"Y_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
				<< std::endl
			<< coordinate_y << std::endl
			<< "</DataArray>" << std::endl
			<< "<DataArray type=\"Float32\" Name=\"Z_COORDINATES\"  NumberOfComponents=\"1\" format=\"ascii\">"
				<< std::endl
			<< coordinate_z << std::endl
			<< "</DataArray>" << std::endl
			<< "</Coordinates>" << std::endl
			
			<< "</Piece>" << std::endl
			<< "</RectilinearGrid>" << std::endl
			<< "</VTKFile>" << std::endl;

		fb.close();
		
	}
};

} //namespace io
} //namespace nast

#endif
