//===========================================================================
//
// File: CpGrid.cpp
//
// Created: Thu Jun  4 12:55:28 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#include "../CpGrid.hpp"
#include "CpGridData.hpp"
#include "CpGridData.cpp"
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <fstream>
#include <iostream>

namespace Dune
{

    CpGrid::CpGrid()
    : use_unique_boundary_ids_(false),
      data_( new CpGridData(*this) ), current_view_data_(data_)
   {}

    CpGrid::~CpGrid()
    {
        delete data_;
    }

    /// Initialize the grid.
    void CpGrid::init(const Opm::parameter::ParameterGroup& param)
    {
	std::string fileformat = param.get<std::string>("fileformat");
	if (fileformat == "sintef_legacy") {
	    std::string grid_prefix = param.get<std::string>("grid_prefix");
	    readSintefLegacyFormat(grid_prefix);
	} else if (fileformat == "eclipse") {
	    std::string filename = param.get<std::string>("filename");
	    double z_tolerance = param.getDefault<double>("z_tolerance", 0.0);
	    bool periodic_extension = param.getDefault<bool>("periodic_extension", false);
	    bool turn_normals = param.getDefault<bool>("turn_normals", false);
	    readEclipseFormat(filename, z_tolerance, periodic_extension, turn_normals);
	} else if (fileformat == "cartesian") {
	    array<int, 3> dims = {{ param.getDefault<int>("nx", 1),
				    param.getDefault<int>("ny", 1),
				    param.getDefault<int>("nz", 1) }};
	    array<double, 3> cellsz = {{ param.getDefault<double>("dx", 1.0),
					 param.getDefault<double>("dy", 1.0),
					 param.getDefault<double>("dz", 1.0) }};
	    createCartesian(dims, cellsz);
	} else {
	    OPM_THROW(std::runtime_error, "Unknown file format string: " << fileformat);
	}
    }




    void CpGrid::createCartesian(const array<int, 3>& dims,
				 const array<double, 3>& cellsize)
    {
	// Make the grdecl format arrays.
	// Pillar coords.
	std::vector<double> coord;
	coord.reserve(6*(dims[0] + 1)*(dims[1] + 1));
	double bot = 0.0;
	double top = dims[2]*cellsize[2];
	// i runs fastest for the pillars.
	for (int j = 0; j < dims[1] + 1; ++j) {
	    double y = j*cellsize[1];
	    for (int i = 0; i < dims[0] + 1; ++i) {
		double x = i*cellsize[0];
		double pillar[6] = { x, y, bot, x, y, top };
		coord.insert(coord.end(), pillar, pillar + 6);
	    }
	}
	std::vector<double> zcorn(8*dims[0]*dims[1]*dims[2]);
	const int num_per_layer = 4*dims[0]*dims[1];
	double* offset = &zcorn[0];
	for (int k = 0; k < dims[2]; ++k) {
	    double zlow = k*cellsize[2];
	    std::fill(offset, offset + num_per_layer, zlow);
	    offset += num_per_layer;
	    double zhigh = (k+1)*cellsize[2];
	    std::fill(offset, offset + num_per_layer, zhigh);
	    offset += num_per_layer;
	}
	std::vector<int> actnum(dims[0]*dims[1]*dims[2], 1);

	// Process them.
	grdecl g;
	g.dims[0] = dims[0];
	g.dims[1] = dims[1];
	g.dims[2] = dims[2];
	g.coord = &coord[0];
	g.zcorn = &zcorn[0];
	g.actnum = &actnum[0];
	current_view_data->processEclipseFormat(g, 0.0, false, false);
    }
} // namespace Dune
