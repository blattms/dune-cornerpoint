//===========================================================================
//
// File: Entity.hpp
//
// Created: Fri May 29 20:26:48 2009
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

#ifndef OPM_ENTITY_HEADER
#define OPM_ENTITY_HEADER

#include <dune/geometry/type.hh>
#include <dune/grid/common/gridenums.hh>
#include "EntityRep.hpp"

namespace Dune
{
    namespace cpgrid
    {


	/// @brief
	/// @todo Doc me!
	/// @tparam
	template <int codim>
	class Entity : public EntityRep<codim>
	{
	public:
	/// @brief
	/// @todo Doc me!
	    enum { codimension = codim };
	    enum { dimension = 3 };
	    enum { mydimension = dimension - codimension };
	    enum { dimensionworld = 3 };


	    typedef EntityPointer<codim> EntityPointer;

	    /// @brief
	    /// @todo Doc me!
	    /// @tparam
	    template <int cd>
	    struct Codim
	    {
		typedef cpgrid::EntityPointer<codim> EntityPointer;
	    };

	    typedef cpgrod::Geometry<3-codim,3> Geometry;
	    typedef Geometry LocalGeometry;

	    typedef cpgrid::LeafIntersectionIterator LeafIntersectionIterator;
	    typedef cpgrid::LevelIntersectionIterator LevelIntersectionIterator;
	    typedef cpgrid::HierarchicIterator HierarchicIterator;

	    typedef double ctype;

	    /// Constructor taking a grid and an integer entity representation.
	    /// This constructor should probably be removed, since it exposes
	    /// details of the implementation of \see EntityRep, see comment in
	    /// EntityRep<>::EntityRep(int).
// 	    Entity(const CpGrid& grid, int entityrep)
// 		: EntityRep<codim>(entityrep), pgrid_(&grid)
// 	    {
// 	    }

	    /// Constructor taking a grid and an entity representation.
	    Entity(const Cpgrid& grid, EntityRep<codim> entityrep)
		: EntityRep<codim>(entityrep), pgrid_(&grid)
	    {
	    }

	    /// Constructor taking a grid, entity index and orientation.
	    Entity(const CpGrid& grid, int index, bool orientation)
		: EntityRep<codim>(index, orientation), pgrid_(&grid)
	    {
	    }

	    /// Equality.
	    bool operator==(const Entity& other) const
	    {
		return EntityRep<codim>::operator==(other)  &&  pgrid_ == other.pgrid_;
	    }

	    /// Inequality.
	    bool operator!=(const Entity& other) const
	    {
		return !operator==(other);
	    }

	    /// Returns the geometry of the entity (does not depend on its orientation).
	    const Geometry& geometry() const
	    {
		return (*pgrid_).template geomVector<codim>()[*this];
	    }

	    /// We do not support refinement, so level() is always 0.
	    int level() const
	    {
		return 0;
	    }

	    /// The entity is always on the leaf grid, since we have no refinement.
	    bool isLeaf() const
	    {
		return true;
	    }

	    /// Refinement is not defined for CpGrid.
	    bool isRegular() const
	    {
		return true;
	    }

	    /// For now, the grid is serial and the only partitionType() is InteriorEntity.
	    PartitionType partitionType() const
	    {
		return InteriorEntity;
	    }

	    /// Using the cube type for all entities now (cells and vertices).
	    GeometryType type() const
	    {
		GeometryType t;
		t.makeCube(3 - codim);
		return t;
	    }

	    /// The count of subentities of codimension cc
	    template <int cc>
	    int count() const
	    {
		static_assert(codim == 0, "");
 		if (cc == 0) {
 		    return 1;
 		} else if (cc == 3) {
		    return 8;
		} else {
		    return 0;
		}
// 		if (cc == 0) {
// 		    return 1;
// 		} else if (cc == 1) {
// 		    return pgrid_->cell_to_face_[*this].size();
// 		} else {
// 		    return pgrid_->cell_to_point_[*this].size();
// 		}
	    }

	    /// Obtain subentity.
	    template <int cc>
	    typename Codim<cc>::EntityPointer subEntity(int i) const
	    {
		static_assert(codim == 0, "");
		if (cc == 0) {
		    assert(i == 0);
		    typename Codim<cc>::EntityPointer se(*pgrid_, EntityRep<codim>::index(), EntityRep<codim>::orientation());
		    return se;
		} else if (cc == 3) {
		    assert(i >= 0 && i < 8);
		    int corner_index = pgrid_->cell_to_point_[EntityRep<codim>::index()][i];
		    typename Codim<cc>::EntityPointer se(*pgrid_, corner_index, true);
		    return se;
		} else {
		    OPM_THROW(std::runtime_error, "No subentity exists of codimension " << cc);
		}
// 		int index = 0;
// 		if (cc == 1) {
// 		    index = pgrid_->cell_to_face_[*this][i].index();
// 		} else if (cc == 3) {
// 		    index = pgrid_->cell_to_point_[*this][i].index();
// 		}
// 		typename Codim<cc>::EntityPointer se(*pgrid_, index); <--- buggy anyway?
// 		return se;
	    }

	    /// Start iterator for the cell-cell intersections of this entity.
	    LevelIntersectionIterator ilevelbegin() const
	    {
		static_assert(codim == 0, "");
		return LevelIntersectionIterator(*pgrid_, *this, false);
	    }

	    /// End iterator for the cell-cell intersections of this entity.
	    LevelIntersectionIterator ilevelend() const
	    {
		static_assert(codim == 0, "");
		return LevelIntersectionIterator(*pgrid_, *this, true);
	    }

	    /// Start iterator for the cell-cell intersections of this entity.
	    LeafIntersectionIterator ileafbegin() const
	    {
		static_assert(codim == 0, "");
		return LeafIntersectionIterator(*pgrid_, *this, false);
	    }

	    /// End iterator for the cell-cell intersections of this entity.
	    typename LeafIntersectionIterator ileafend() const
	    {
		static_assert(codim == 0, "");
		return LeafIntersectionIterator(*pgrid_, *this, true);
	    }

	    /// Dummy first child iterator.
	    typename HierarchicIterator hbegin(int) const
	    {
		return HierarchicIterator(*pgrid_);
	    }

	    /// Dummy beyond last child iterator.
	    typename HierarchicIterator hend(int) const
	    {
		return HierarchicIterator(*pgrid_);
	    }

	    /// \brief Returns true, if the entity has been created during the last call to adapt(). Dummy.
	    bool isNew() const
	    {
		return false;
	    }
  
	    /// \brief Returns true, if entity might disappear during the next call to adapt(). Dummy.
	    bool mightVanish() const
	    {
		return false;
	    }

	    /// No hierarchy in this grid.
	    bool hasFather() const
	    {
		return false;
	    }


	    /// Dummy, returning this.
	    EntityPointer father() const
	    {
		return EntityPointer(*this);
	    }


	    /// Dummy, returning default geometry.
	    LocalGeometry geometryInFather() const
	    {
		return LocalGeometry();
	    }

	    /// Returns true if any of my intersections are on the boundary.
	    /// Implementation note:
	    /// This is a slow, computed, function. Could be speeded
	    /// up by putting boundary info in the CpGrid class.
	    bool hasBoundaryIntersections() const
	    {
		// Copied implementation from EntityDefaultImplementation,
		// except for not checking LevelIntersectionIterators.
		typedef LeafIntersectionIterator IntersectionIterator; 
		IntersectionIterator end = ileafend();
		for (IntersectionIterator it = ileafbegin(); it != end; ++it) {
		    if (it->boundary()) return true;
		}
		return false;
	    }

            // Mimic Dune entity wrapper
            
            const Entity& impl() const
            {
                return *this;
            }
            
            Entity& impl()
            {
                return *this;
            }
	protected:
	    const CpGrid* pgrid_;

	    bool valid() const
	    {
		return EntityRep<codim>::index() < pgrid_->size(codim);
	    }
	};




	/// \brief Class representing a pointer to an entity.
	/// Implementation note:
	/// Since our entities are quite lightweight, we have chosen
	/// to implement EntityPointer by inheritance from
	/// Entity. Thus all dereferencing operators return the object
	/// itself as an Entity.
	template <int codim>
	class EntityPointer : public cpgrid::Entity<codim>
	{
	public:
	    typedef cpgrid::Entity<codim> Entity;

	    /// Construction from entity.
	    explicit EntityPointer(const Entity& e)
		: Entity(e)
	    {
	    }

	    /// Constructor taking a grid and an entity representation.
	    EntityPointer(const CpGrid& grid, EntityRep<codim> entityrep)
		: Entity(grid, entityrep)
	    {
	    }

	    /// Constructor taking a grid, entity index and orientation.
	    EntityPointer(const CpGrid& grid, int index, bool orientation)
		: Entity(grid, index, orientation)
	    {
	    }

// 	    /// Member by pointer operator.
// 	    Entity* operator->()
// 	    {
// 		assert(Entity::valid());
// 		return this;
// 	    }

// 	    /// Const member by pointer operator.
// 	    const Entity* operator->() const
// 	    {
// 		assert(Entity::valid());
// 		return this;
// 	    }

// 	    /// Dereferencing operator.
// 	    Entity& operator*()
// 	    {
// 		assert(Entity::valid());
// 		return *this;
// 	    }

// 	    /// Const dereferencing operator.
// 	    const Entity& operator*() const
// 	    {
// 		assert(Entity::valid());
// 		return *this;
// 	    }

	    /// Const member by pointer operator.
	    Entity* operator->() const
	    {
		assert(Entity::valid());
		return const_cast<EntityPointer*>(this); // const_cast-hack added because of error in vtkwriter.hh
	    }

	    /// Const dereferencing operator.
	    Entity& operator*() const
	    {
		assert(Entity::valid());
		return const_cast<EntityPointer&>(*this); // const_cast-hack added because of error in vtkwriter.hh
	    }


	    /// Minimizes memory usage.
	    /// Nothing to do in our case.
	    void compactify()
	    {
	    }
	};




    } // namespace cpgrid
} // namespace Dune


#endif // OPM_ENTITY_HEADER
