//===========================================================================
//
// File: LookUpData.hh
//
// Created: Tue May 23 14:44:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
Copyright 2023 Equinor ASA.

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

#include <dune/grid/common/mcmgmapper.hh>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp> // "CartesianIndexMapper not specialized for given grid"

namespace Dune
{
template <typename Grid, typename GridView>
class LookUpData
{
public:
    // Constructor taking a Grid object
    LookUpData(const Grid& grid) :
        gridView_(grid.leafGridView()),
        elemMapper_(gridView_, Dune::mcmgElementLayout()),
        cartMapper_(grid)
    {
    }

    // Constructor taking a GridView, ElementMapper, CartesianMapper
    //template<typename GridView>
    LookUpData(const  GridView& gridView,
               const  Dune::MultipleCodimMultipleGeomTypeMapper<GridView>& elemMapper,
               const  Dune::CartesianIndexMapper<Grid>& cartMapper) :
        gridView_(gridView),
        elemMapper_(gridView_, Dune::mcmgElementLayout()),
        cartMapper_(gridView_.grid())
    {
    }

     // Constructor taking a GridView
    //template <typename GridView>
    LookUpData(const GridView& gridView) :
        gridView_(gridView)
    {
    }

    template<typename EntityType, typename FeatureType>
    FeatureType operator()(const EntityType& elem, const std::vector<FeatureType>& feature_vec) const
    {
        assert(0 <= elemMapper_.index(elem) && static_cast<int>(feature_vec.size()) > elemMapper_.index(elem));
        // Assuming feature is given for gridView_
        return feature_vec[elemMapper_.index(elem)];
    }


    template<typename FeatureType>
    FeatureType operator()(const int& elemIdx, const std::vector<FeatureType> feature_vec) const
    {
        const int& cartIdx = cartMapper_.cartesianIndex(elemIdx); 
        assert(0 <= cartIdx && static_cast<int>(feature_vec.size()) > cartIdx);
        return feature_vec[cartIdx]; 
    }

   
    // getOriginIdx() For general grids: retunrs a copy of the same index.
    //                For CpGrid: returns index of origin cell (parent cell or equivalent cell when no father) in level 0
    int getOriginIndex(const int& elemIdx) // elemIdx is supposed to be an index of a leafview cell
    {
        if (std::is_same<Grid,Dune::CpGrid>::value) {
            const Dune::cpgrid::Entity<0>& elem = Dune::cpgrid::Entity<0>(gridView_, elemIdx, true);
            return elem.getOrigin().index();
        }
        else{
            return elemIdx;
        }
    }
    
protected:
    GridView gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    Dune::CartesianIndexMapper<Grid> cartMapper_; 


}; // end LookUpData class
}
// end namespace Dune
