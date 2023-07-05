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

namespace Dune
{
template <typename Grid, typename GridView>
class LookUpData
{
public:
    // Constructor taking a Grid object
    LookUpData(const GridView& gridView,
               const CartesianIndexMapper<Grid>& mapper) :
        gridView_(gridView),
        elemMapper_(gridView_, Dune::mcmgElementLayout()),
        cartMapper_(&mapper)
    {
    }

    // Constructor taking a GridView
    LookUpData(const  GridView& gridView) :
        gridView_(gridView),
        elemMapper_(gridView_, Dune::mcmgElementLayout()),
        cartMapper_()
    {
    }

    // operator()(Entity, Vector) Call operator taking an EntityObject and a FeatureVector.
    //                            Return feature of the entity, via (ACTIVE) INDEX
    //                            For general grids, the feature vector is given for the gridView_.
    //                            [For CpGrid, the feature vector is given for level 0.]
    template<typename EntityType, typename FeatureType>
    FeatureType operator()(const EntityType& elem, const std::vector<FeatureType>& feature_vec) const;

    // operator()(EntityIndex, Vector) Call operator taking an EntityObject and a FeatureVector.
    //                                 Return feature of the entity, via CARTESIAN INDEX
    //                                 For general grids, the feature vector is given for the gridView_.
    //                                 [For CpGrid, the feature vector is given for level 0.]
    template<typename FeatureType>
    FeatureType operator()(const int& elemIdx, const std::vector<FeatureType> feature_vec) const;

    // getOriginIdx() For general grids: retunrs a copy of the same index.
    //                [For CpGrid: returns index of origin cell (parent cell or equivalent cell when no father) in level 0]
    int getOriginIndex(const int& elemIdx) const; // elemIdx is supposed to be an index of a leafview cell

protected:
    const GridView& gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView> elemMapper_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;

}; // end LookUpData class
}
// end namespace Dune

template<typename Grid, typename GridView>
template<typename EntityType, typename FeatureType>
FeatureType Dune::LookUpData<Grid,GridView>::operator()(const EntityType& elem, const std::vector<FeatureType>& feature_vec) const
{
    assert(0 <= elemMapper_.index(elem) && static_cast<int>(feature_vec.size()) > elemMapper_.index(elem));
    // Assuming feature is given for gridView_
    return feature_vec[elemMapper_.index(elem)];
}

template<typename Grid, typename GridView>
template<typename FeatureType>
FeatureType Dune::LookUpData<Grid,GridView>::operator()(const int& elemIdx, const std::vector<FeatureType> feature_vec) const
{
    assert(cartMapper_);
    const int& cartIdx = cartMapper_.cartesianIndex(elemIdx);
    assert(0 <= cartIdx && static_cast<int>(feature_vec.size()) > cartIdx);
    return feature_vec[cartIdx];
}

template<typename Grid, typename GridView>
int Dune::LookUpData<Grid,GridView>::getOriginIndex(const int& elemIdx) const
{
    return elemIdx;
}
