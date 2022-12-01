
/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2022 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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

#define BOOST_TEST_MODULE GeometryTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>

#include <sstream>
#include <iostream>

void check_refinedPatch_grid(const std::array<int,3>& cells_per_dim,
                        const std::array<int,3>& start_ijk,
                        const std::array<int,3>& end_ijk,
                             Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<3, 3>,0>& refined_cells,
                             Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<2,3>,1>& refined_faces,
                             Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>>& refined_corners)
{
    const std::array<int,3> patch_dim = {end_ijk[0]-start_ijk[0], end_ijk[1]-start_ijk[1], end_ijk[2]-start_ijk[2]};
    if ((patch_dim[0] == 0) || (patch_dim[1] == 0) || (patch_dim[2] == 0)) {
                    OPM_THROW(std::logic_error, "Empty patch. Cannot convert patch into cell.");
    }
    // Check amount of refined faces.
    int count_faces = (cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*((cells_per_dim[2]*patch_dim[2])+1)) // 'bottom/top faces'
        +  (((cells_per_dim[0]*patch_dim[0])+1)*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2]) // 'front/back faces'
        + (cells_per_dim[0]*patch_dim[0]*((cells_per_dim[1]*patch_dim[1]) +1)*cells_per_dim[2]*patch_dim[2]);  // 'left/right faces'
    BOOST_CHECK_EQUAL(refined_faces.size(), count_faces);
    // Check amount of refined corners.
    int count_corners = ((cells_per_dim[0]*patch_dim[0])+1)*((cells_per_dim[1]*patch_dim[1])+1)*((cells_per_dim[2]*patch_dim[2])+1);
    BOOST_CHECK_EQUAL(refined_corners.size(), count_corners);

    int count_cells = cells_per_dim[0]*patch_dim[0]*cells_per_dim[1]*patch_dim[1]*cells_per_dim[2]*patch_dim[2];
    BOOST_CHECK_EQUAL(refined_cells.size(), count_cells);
}


void refinePatch_and_check(const std::array<int, 3>& cells_per_dim,
                           const std::array<int,3>& start_ijk,
                           const std::array<int,3>& end_ijk)
{
    using Dune::cpgrid::DefaultGeometryPolicy;
    Dune::CpGrid refined_grid;
    auto& child_view_data = *refined_grid.current_view_data_;
    Dune::cpgrid::OrientedEntityTable<0, 1>& cell_to_face = child_view_data.cell_to_face_;
    Opm::SparseTable<int>& face_to_point = child_view_data.face_to_point_;
    Dune::cpgrid::DefaultGeometryPolicy& geometries = child_view_data.geometry_;
    std::vector<std::array<int, 8>>& cell_to_point = child_view_data.cell_to_point_;
    Dune::cpgrid::OrientedEntityTable<1,0>& face_to_cell = child_view_data.face_to_cell_;
    Dune::cpgrid::EntityVariable<enum face_tag, 1>& face_tags = child_view_data.face_tag_;
    Dune::cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>, 1>& face_normals = child_view_data.face_normals_;

     // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};   
    std::array<int, 3> start_ijk_A = {1,0,1};
    std::array<int, 3> end_ijk_A = {3,2,3};  // then patch_dim = {3-1, 2-0, 3-1} ={0,0,0}
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    
    // Call refinedBlockPatch()
    coarse_grid.current_view_data_->refineBlockPatch(cells_per_dim_patch, start_ijk_A, end_ijk_A);
    // Create a pointer pointing at the CpGridData object coarse_grid.current_view_data_.
    std::shared_ptr<Dune::cpgrid::CpGridData> coarse_grid_ptr =  std::make_shared<Dune::cpgrid::CpGridData>();
    *coarse_grid.current_view_data_ = *coarse_grid_ptr;
    // Create a vector of shared pointers of CpGridData type.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> data;
    // Add coarse_grid_ptr to data.
    data.push_back(coarse_grid_ptr);
    
    // Call getLeafView2Levels()
    coarse_grid.getLeafView2Levels(data, cells_per_dim_patch, start_ijk_A, end_ijk_A);
    
    // Call addLevel()
    const int level_to_refine = 0;
    std::vector<std::array<int,2>> future_leaf_corners;
    std::vector<std::array<int,2>> future_leaf_faces;
    std::vector<std::array<int,2>> future_leaf_cells;
    coarse_grid.addLevel(data, level_to_refine, cells_per_dim_patch, start_ijk_A, end_ijk_A,
    future_leaf_corners, future_leaf_faces, future_leaf_cells);

    Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<3, 3>,0> refined_cells;
    Dune::cpgrid::EntityVariable<Dune::cpgrid::Geometry<2,3>,1> refined_faces;
    Dune::cpgrid::EntityVariableBase<Dune::cpgrid::Geometry<0,3>> refined_corners;
    check_refinedPatch_grid(cells_per_dim_patch, start_ijk_A, end_ijk_A, refined_cells, refined_faces, refined_corners);
                            //geometries.template geomVector<0>(), geometries.template geomVector<1>(),
                            //geometries.template geomVector<3>());

    /*  cpgrid::OrientedEntityTable<1,0> face_to_cell_computed;
    cell_to_face.makeInverseRelation(face_to_cell_computed);
    BOOST_CHECK(face_to_cell_computed == face_to_cell); */
} 

/*BOOST_AUTO_TEST_CASE(refine_patch)
{
     // Create a grid
    Dune::CpGrid coarse_grid;
    std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    std::array<int, 3> grid_dim = {4,3,3};
    std::array<int, 3> cells_per_dim_patch = {2,2,2};   
    std::array<int, 3> start_ijk = {1,0,1};
    std::array<int, 3> end_ijk = {3,2,3};  // then patch_dim = {3-1, 2-0, 3-1} ={0,0,0}
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    // Call refinedBlockPatch()
    coarse_grid.current_view_data_->refineBlockPatch(cells_per_dim_patch, start_ijk, end_ijk);
    // Create a pointer pointing at the CpGridData object coarse_grid.current_view_data_.
    std::shared_ptr<Dune::cpgrid::CpGridData> coarse_grid_ptr =  std::make_shared<Dune::cpgrid::CpGridData>();
    *coarse_grid.current_view_data_ = *coarse_grid_ptr;
    // Create a vector of shared pointers of CpGridData type.
    std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> data;
    // Add coarse_grid_ptr to data.
    data.push_back(coarse_grid_ptr);
    // Call getLeafView2Levels()
    coarse_grid.getLeafView2Levels(data, cells_per_dim, start_ijk, end_ijk);
    // Call addLevel()
    const int level_to_refine = 0;
    std::vector<std::array<int,2>> future_leaf_corners;
    std::vector<std::array<int,2>> future_leaf_faces;
    std::vector<std::array<int,2>> future_leaf_cells;
    coarse_grid.addLevel(data, level_to_refine, cells_per_dim_pacth, start_ijk, end_ijk,
    future_leaf_corners, future_leaf_faces, future_leaf_cells);

    refinePatch_and_check(cells_per_dim_patch, start_ijk, end_ijk);
    }*/
