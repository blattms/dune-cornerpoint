//===========================================================================
//
// File: LookUpDataCpGrid.hh
//
// Created: Tue May 25 11:45:00 2023
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

namespace Dune
{
/// Specialization for CpGrid
template<typename GridView>
class LookUpData<Dune::CpGrid, GridView>
{
public:
    // Constructor taking a CpGrid object
    LookUpData(const Dune::CpGrid&){
    }

    // operator()(Entity, Vector) Call operator taking an Entity and a FeatureVector.
    //                            Return feature of the entity, via CARTESIAN INDEX
    //                            For CpGrid, the feature vector is given for level 0.
    //                            [For general grids, the feature vector is given for the gridView_.]
    template<typename FeatureType>
    FeatureType operator()(const Dune::cpgrid::Entity<0>& elem, const std::vector<FeatureType>& feature_vec) const;

    // getOriginIdx() For CpGrid: returns index of origin cell (parent cell or equivalent cell when no father) in level 0
    //                [For general grids: retunrs a copy of the same index.]
    int getOriginIndex(const int& elemIdx) // elemIdx is supposed to be an index of a leafview cell

template<typename FeatureType>
FeatureType Dune::LookUpData<Dune::CpGrid, GridView>::operator()(const Dune::cpgrid::Entity<0>& elem,
                                                                 const std::vector<FeatureType>& feature_vec) const
    {
    // elem.getOrigin() Get entity in level 0 (either parent cell, or equivalent cell, or 'itself' if grid_ = level 0)
    return feature_vec[elem.getOrigin().index()];
    }

}; // end LookUpData<CpGrid> class
int  Dune::LookUpData<Dune::CpGrid, GridView>::getOriginIndex(const int& elemIdx) const
{
    const Dune::cpgrid::Entity<0>& elem = Dune::cpgrid::Entity<0>(this->gridView_, elemIdx, true);
    return elem.getOrigin().index();
}
// end namespace Dune
