//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri May 29 20:26:36 2009
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

#ifndef OPM_CPGRID_HEADER
#define OPM_CPGRID_HEADER

#include <string>
#include <map>
#include <array>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#include "cpgrid/Entity.hpp"
#include "cpgrid/Geometry.hpp"
#include "cpgrid/Intersection.hpp"
#include "cpgrid/Iterators.hpp"
#include "cpgrid/Indexsets.hpp"
#include "cpgrid/DefaultGeometryPolicy.hpp"
#include <opm/core/grid/cpgpreprocess/preprocess.h>

#include <iostream>


namespace Opm {
    class EclipseGridParser;
    namespace parameter {
	class ParameterGroup;
    }
}

namespace Dune
{

    class CpGrid;
    
    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridTraits
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridTraits
    {
	/// \brief The type that implements the grid.
	typedef CpGrid Grid;

	/// \brief The type of the intersection at the leafs of the grid.
	typedef cpgrid::Intersection LeafIntersection;
	/// \brief The type of the intersection at the levels of the grid.
	typedef cpgrid::Intersection LevelIntersection;
	/// \brief The type of the intersection iterator at the leafs of the grid.
	typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
	/// \brief The type of the intersection iterator at the levels of the grid.
	typedef cpgrid::IntersectionIterator LevelIntersectionIterator;

	/// \brief The type of the  hierarchic iterator.
	typedef cpgrid::HierarchicIterator HierarchicIterator;

	/// \brief Traits associated with a specific codim.
	/// \tparam cd The codimension.
	template <int cd>
	struct Codim
	{
	    /// \brief The type of the geometry associated with the entity.
	    /// IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
 	    typedef cpgrid::Geometry<3-cd, 3> Geometry;
	    //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> Geometry;
	    /// \brief The type of the local geometry associated with the entity.
 	    typedef cpgrid::Geometry<3-cd, 3> LocalGeometry;
	    //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> LocalGeometry;
	    /// \brief The type of the entity.
	    typedef cpgrid::Entity<cd> Entity;

	    /// \brief The type of the iterator over all level entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition> LevelIterator;

	    /// \brief The type of the iterator over all leaf entities of this codim.
	    typedef cpgrid::Iterator<cd, All_Partition> LeafIterator;

	    /// \brief The type of the entity pointer for entities of this codim.
	    typedef cpgrid::EntityPointer<cd> EntityPointer;

	    /// \brief The type of the entity pointer for entities of this codim.
	    typedef cpgrid::EntityPointer<cd> EntitySeed;

	    /// \brief Traits associated with a specific grid partition type.
	    /// \tparam pitype The type of the grid partition.
	    template <PartitionIteratorType pitype>
	    struct Partition
	    {
		/// \brief The type of the iterator over the level entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype> LevelIterator;
		/// \brief The type of the iterator over the leaf entities of this codim on this partition.
		typedef cpgrid::Iterator<cd, pitype> LeafIterator;
	    };
	};

	/// \brief Traits associated with a specific grid partition type.
	/// \tparam pitype The type of the grid partition.
	template <PartitionIteratorType pitype>
	struct Partition
	{
	    /// \brief The type of the level grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LevelGridView;
	    typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid, pitype> > LevelGridView;

	    /// \brief The type of the leaf grid view associated with this partition type.
// 	    typedef cpgrid::GridView<pitype> LeafGridView;
	    typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid, pitype> > LeafGridView;
	};

	/// \brief The type of the level index set.
	typedef cpgrid::IndexSet LevelIndexSet;
	/// \brief The type of the leaf index set.
	typedef cpgrid::IndexSet LeafIndexSet;
	/// \brief The type of the global id set.
	typedef cpgrid::IdSet GlobalIdSet;
	/// \brief The type of the local id set.
	typedef cpgrid::IdSet LocalIdSet;

	/// \brief The type of the collective communication.

    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    };



    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridFamily
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridFamily
    {
	typedef CpGridTraits Traits;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGrid
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> provides \ref Dune::Grid </em>]
    class CpGrid
	: public GridDefaultImplementation<3, 3, double, CpGridFamily >
    {
    public:

	// --- Typedefs ---


	/// Family typedef, why is this not defined by Grid<>?
	typedef CpGridFamily GridFamily;


	// --- Methods ---


	/// Default constructor
	CpGrid()
	    : ccobj_(Dune::MPIHelper::getCommunicator()),
          index_set_(*this),
          use_unique_boundary_ids_(false),
          idSet_( *this )
	{
	}


        /// Initialize the grid from parameters.
	void init(const Opm::parameter::ParameterGroup& param);


	/// Read the Sintef legacy grid format ('topogeom').
	/// \param grid_prefix the grid name, such that topology is
	/// found in <grid_prefix>-topo.dat etc.
	void readSintefLegacyFormat(const std::string& grid_prefix);

	/// Write the Sintef legacy grid format ('topogeom').
	/// \param grid_prefix the grid name, such that topology will be
	/// found in <grid_prefix>-topo.dat etc.
        void writeSintefLegacyFormat(const std::string& grid_prefix) const;

	/// Read the Eclipse grid format ('grdecl').
	/// \param filename the name of the file to read.
	/// \param z_tolerance points along a pillar that are closer together in z
	///        coordinate than this parameter, will be replaced by a single point.
	/// \param periodic_extension if true, the grid will be (possibly) refined, so that
	///        intersections/faces along i and j boundaries will match those on the other
	///        side. That is, i- faces will match i+ faces etc.
	void readEclipseFormat(const std::string& filename, double z_tolerance, bool periodic_extension, bool turn_normals = false);

	/// Read the Eclipse grid format ('grdecl').
	/// \param input_data the data contained in a parser object.
	/// \param z_tolerance points along a pillar that are closer together in z
	///        coordinate than this parameter, will be replaced by a single point.
	/// \param periodic_extension if true, the grid will be (possibly) refined, so that
	///        intersections/faces along i and j boundaries will match those on the other
	///        side. That is, i- faces will match i+ faces etc.
	/// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
	/// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
	void processEclipseFormat(const Opm::EclipseGridParser& input_parser, double z_tolerance, bool periodic_extension, bool turn_normals = false, bool clip_z = false);

	/// Read the Eclipse grid format ('grdecl').
	/// \param input_data the data in grdecl format, declared in preprocess.h.
	/// \param z_tolerance points along a pillar that are closer together in z
	///        coordinate than this parameter, will be replaced by a single point.
	/// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
	void processEclipseFormat(const grdecl& input_data, double z_tolerance, bool remove_ij_boundary, bool turn_normals = false);

	/// Create a cartesian grid.
	/// \param dims the number of cells in each cartesian direction.
	/// \param cellsize the size of each cell in each dimension.
	void createCartesian(const array<int, 3>& dims,
			     const array<double, 3>& cellsize);

	/// The logical cartesian size of the grid.
	/// This function is not part of the Dune grid interface,
	/// and should be used with caution.
        const std::array<int, 3>& logicalCartesianSize() const
        {
            return logical_cartesian_size_;
        }

	/// Retrieve mapping from internal ("compressed") active grid
	/// cells to external ("uncompressed") cells.  Specifically,
	/// @code globalCell()[i] @endcode is the linearized Cartesian
	/// index of grid cell @code i @endcode.  This method should
	/// only be used by classes which really need it, such as
	/// those dealing with permeability fields from the input deck
	/// from whence the current CpGrid was constructed.
        const std::vector<int>& globalCell() const
        {
            return global_cell_;
        }

        /// @brief
        ///    Extract Cartesian index triplet (i,j,k) of an active cell.
        ///
        /// @param [in] c
        ///    Active cell index.
        ///
        /// @param [out] ijk  Cartesian index triplet
        void getIJK(const int c, std::array<int,3>& ijk) const
        {
            int gc = global_cell_[c];
            ijk[0] = gc % logical_cartesian_size_[0];  gc /= logical_cartesian_size_[0];
            ijk[1] = gc % logical_cartesian_size_[1];
            ijk[2] = gc / logical_cartesian_size_[1];
        }

	/// Is the grid currently using unique boundary ids?
	/// \return true if each boundary intersection has a unique id
	///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
	bool uniqueBoundaryIds() const
	{
	    return use_unique_boundary_ids_;
	}

	/// Set whether we want to have unique boundary ids.
	/// \param uids if true, each boundary intersection will have a unique boundary id.
	void setUniqueBoundaryIds(bool uids)
	{
	    use_unique_boundary_ids_ = uids;
	    if (use_unique_boundary_ids_ && unique_boundary_ids_.empty()) {
		computeUniqueBoundaryIds();
	    }
	}

	// --- Dune interface below ---


        /// Return grid name. It's the same as the class name.
	/// What did you expect, something funny?
        std::string name() const
        {
            return "CpGrid";
        }


        /// Return maximum level defined in this grid. Levels are numbered
        /// 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const
	{
            return 0;
        }


        /// Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim, All_Partition>(*this, 0, true);
        }


        /// one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,All_Partition>(*this, size(codim), true );
        }


        /// Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*this, 0, true );
        }


        /// one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*this, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim, All_Partition>(*this, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim, All_Partition>(*this, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const
	{
            return cpgrid::Iterator<codim, PiType>(*this, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const
	{
            return cpgrid::Iterator<codim, PiType>(*this, size(codim), true);
        }


        /// \brief Number of grid entities per level and codim
        int size (int level, int codim) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(codim);
        }


        /// number of leaf entities per codim in this process
        int size (int codim) const
	{
	    switch (codim) {
	    case 0: return cell_to_face_.size();
	    case 1: return 0;
	    case 2: return 0;
	    case 3: return geomVector<3>().size();
	    default: return 0;
	    }
        }


        /// number of entities per level and geometry type in this process
        int size (int level, GeometryType type) const
	{
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(type);
        }


        /// number of leaf entities per geometry type in this process
        int size (GeometryType type) const
        {
	    if (type.isCube()) {
		return size(3 - type.dim());
	    } else {
		return 0;
	    }
	}


        /// \brief Access to the GlobalIdSet
        const Traits::GlobalIdSet& globalIdSet() const
	{
            return idSet_;
        }


        /// \brief Access to the LocalIdSet
        const Traits::LocalIdSet& localIdSet() const
	{
            return idSet_;
        }


        /// \brief Access to the LevelIndexSets
        const Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return index_set_;
        }


        /// \brief Access to the LeafIndexSet
        const Traits::LeafIndexSet& leafIndexSet() const
        {
            return index_set_;
        }


        /// global refinement
        void globalRefine (int refCount)
        {
            std::cout << "Warning: Global refinement not implemented, yet." << std::endl;
        }

        const std::vector< Dune :: GeometryType >& geomTypes( const int codim ) const 
        {
          return leafIndexSet().geomTypes( codim );
        }

	/*  No refinement implemented. GridDefaultImplementation's methods will be used.

        /// \brief Mark entity for refinement
	///
	/// This only works for entities of codim 0.
	/// The parameter is currently ignored
	///
	/// \return <ul>
	/// <li> true, if marking was succesfull </li>
	/// <li> false, if marking was not possible </li>
	/// </ul>

        bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
        {
            return hostgrid_->mark(refCount, getHostEntity<0>(*e));
        }

        /// \brief Return refinement mark for entity
	///
	/// \return refinement mark (1,0,-1)

        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            return hostgrid_->getMark(getHostEntity<0>(*e));
        }

        /// \todo Please doc me !
        bool preAdapt() {
            return hostgrid_->preAdapt();
        }


        /// Triggers the grid refinement process
        bool adapt()
        {
            return hostgrid_->adapt();
        }

        /// \brief Clean up refinement markers
        void postAdapt() {
            return hostgrid_->postAdapt();
        }

	end of refinement section */


	/* No parallelism implemented.  GridDefaultImplementation's methods will be used.

        /// \brief Size of the overlap on the leaf level
        unsigned int overlapSize(int codim) const {
            return hostgrid_->overlapSize(codim);
        }


        /// \brief Size of the ghost cell layer on the leaf level
        unsigned int ghostSize(int codim) const {
            return hostgrid_->ghostSize(codim);
        }


        /// \brief Size of the overlap on a given level
        unsigned int overlapSize(int level, int codim) const {
            return hostgrid_->overlapSize(level,codim);
        }


        /// \brief Size of the ghost cell layer on a given level
        unsigned int ghostSize(int level, int codim) const {
            return hostgrid_->ghostSize(level,codim);
        }


        /// \brief Distributes this grid over the available nodes in a distributed machine
	///
	/// \param minlevel The coarsest grid level that gets distributed
	/// \param maxlevel does currently get ignored
        void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
            DUNE_THROW(NotImplemented, "CpGrid::loadBalance()");
        }

        /// \brief The communication interface
	///  @param T: array class holding data associated with the entities
	///  @param P: type used to gather/scatter data in and out of the message buffer
	///  @param codim: communicate entites of given codim
	///  @param if: one of the predifined interface types, throws error if it is not implemented
	///  @param level: communicate for entities on the given level
	///
	///  Implements a generic communication function sending an object of type P for each entity
	///  in the intersection of two processors. P has two methods gather and scatter that implement
	///  the protocol. Therefore P is called the "protocol class".
        template<class T, template<class> class P, int codim>
        void communicate (T& t, InterfaceType iftype, CommunicationDirection dir, int level);

        /// The new communication interface.
	/// communicate objects for all codims on a given level
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int level) const
        {}

        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {}

	end of parallel section */

        /// dummy collective communication
        const CollectiveCommunication& comm () const
        {
            return ccobj_;
        }

        // ------------ End of Dune interface, start of simplified interface --------------

        // enum { dimension = 3 }; // already defined

        typedef Dune::FieldVector<double, 3> Vector;

        // Topology
        int numCells() const
        {
            return cell_to_face_.size();
        }
        int numFaces() const
        {
            return face_to_cell_.size();
        }
        int numVertices() const
        {
            return geomVector<3>().size();
        }

        int numCellFaces(int cell) const
        {
            return cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
        }
        int cellFace(int cell, int local_index) const
        {
            return cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
        }
        int faceCell(int face, int local_index) const
        {
            cpgrid::OrientedEntityTable<1,0>::row_type r
                = face_to_cell_[cpgrid::EntityRep<1>(face, true)];
            bool a = (local_index == 0);
            bool b = r[0].orientation();
            bool use_first = a ? b : !b;
            if (r.size() == 2) {
                assert(r[0].orientation() != r[1].orientation());
                return use_first ? r[0].index() : r[1].index();
            } else {
                return use_first ? r[0].index() : -1;
            }
        }
        int numFaceVertices(int face) const
        {
            return face_to_point_[face].size();
        }
        int faceVertex(int face, int local_index) const
        {
            return face_to_point_[face][local_index];
        }

        // Geometry
        const Vector& vertexPosition(int vertex) const
        {
            return geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
        }
        double faceArea(int face) const
        {
            return geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
        }
        const Vector& faceCentroid(int face) const
        {
            return geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
        }
        const Vector& faceNormal(int face) const
        {
            return face_normals_.get(face);
        }
        double cellVolume(int cell) const
        {
            return geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
        }
        const Vector& cellCentroid(int cell) const
        {
            return geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
        }

        // Extra
        int boundaryId(int face) const
        {
            int ret = 0;
            cpgrid::EntityRep<1> f(face, true);
            if (face_to_cell_[f].size() == 1) {
                if (uniqueBoundaryIds()) {
                    // Use the unique boundary ids.
                    ret = unique_boundary_ids_[f];
                } else {
                    // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
                    const bool normal_is_in = !(face_to_cell_[f][0].orientation());
                    enum face_tag tag = face_tag_[f];
                    switch (tag) {
                    case LEFT:
                        //                   LEFT : RIGHT
                        ret = normal_is_in ? 1    : 2; // min(I) : max(I)
                        break;
                    case BACK:
                        //                   BACK : FRONT
                        ret = normal_is_in ? 3    : 4; // min(J) : max(J)
                        break;
                    case TOP:
                        // Note: TOP at min(K) as 'z' measures *depth*.
                        //                   TOP  : BOTTOM
                        ret = normal_is_in ? 5    : 6; // min(K) : max(K)
                        break;
                    }
                }
            }
            return ret;
        }


        // ------------ End of simplified interface --------------

    private:

	// --------- Friends ---------

	template <int cd, class GridType>
	friend class cpgrid::Entity;

	template <class GridType>
	friend class cpgrid::Intersection;

	template <class GridType>
	friend class cpgrid::IndexSet;


	// --------- Data members ---------

	// Dune-ish stuff
        CollectiveCommunication ccobj_;
	cpgrid::IndexSet index_set_;

	// Representing the topology
	cpgrid::OrientedEntityTable<0, 1> cell_to_face_;
	cpgrid::OrientedEntityTable<1, 0> face_to_cell_;
	Opm::SparseTable<int> face_to_point_;
	std::vector< array<int,8> > cell_to_point_;
        // Size of Cartesian bounding box.
	std::array<int, 3> logical_cartesian_size_;
        std::vector<int>                  global_cell_;
        cpgrid::EntityVariable<enum face_tag, 1> face_tag_; // {LEFT, BACK, TOP}
        
	// Representing geometry
	typedef cpgrid::DefaultGeometryPolicy Geom;
	Geom geometry_;
	typedef FieldVector<ctype, 3> PointType;
	cpgrid::SignedEntityVariable<PointType, 1> face_normals_;
	std::vector<PointType> allcorners_; // Yes, this is already stored in the point geometries. \TODO Improve.
	// Boundary information (optional).
	bool use_unique_boundary_ids_;
	cpgrid::EntityVariable<int, 1> unique_boundary_ids_;

  const Traits :: LocalIdSet idSet_;

	// --------- Methods ---------

	// Return the geometry vector corresponding to the given codim.
	template <int codim>
	const cpgrid::EntityVariable< cpgrid::Geometry<3 - codim, 3>, codim>& geomVector() const
	{
	    return geometry_.geomVector<codim>();
	}

	// Make unique boundary ids for all intersections.
	void computeUniqueBoundaryIds();

    }; // end Class CpGrid



    namespace Capabilities
    {
	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 0>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
	template <>
	struct hasEntity<CpGrid, 3>
	{
	    static const bool v = true;
	};

	/// \todo Please doc me !
	template <>
	struct isParallel<CpGrid>
	{
	    static const bool v = false;
	};

	/// \todo Please doc me !
	template <>
	struct hasBackupRestoreFacilities<CpGrid>
	{
	    static const bool v = false;
	};

    }

} // namespace Dune

#include <dune/grid/cpgrid/PersistentContainer.hpp>
#endif // OPM_CPGRID_HEADER
