//===========================================================================
//
// File: Geometry.hpp
//
// Created: Fri May 29 23:29:24 2009
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
  Copyright 2009, 2010, 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2011, 2022 Equinor ASA.

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

#ifndef OPM_GEOMETRY_HEADER
#define OPM_GEOMETRY_HEADER

#include <cmath>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/geometry/type.hh>

#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/common/Volumes.hpp>
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/grid/utility/ErrorMacros.hpp>

namespace Dune
{
    namespace cpgrid
    {

        /// This class encapsulates geometry for vertices,
        /// intersections, and cells. The main template is empty,
        /// the actual dim == 3 (cell), dim == 2 (intersection),
        /// and dim == 0 (vertex) cases have specializations.
        /// For vertices and cells we use the cube type, and provide
        /// constant (vertex) or trilinear (cell) mappings.
        /// For intersections, we use the singular geometry type
        /// (None), and provide no mappings.
        template <int mydim, int cdim>
        class Geometry
        {
        };




        /// Specialization for 0 dimensional geometries, i.e. vertices.
        template <int cdim> // GridImp arg never used
        class Geometry<0, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 0};
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;


            /// @brief Construct from vertex position
            /// @param pos the position of the vertex
            Geometry(const GlobalCoordinate& pos)
                : pos_(pos)
            {
            }

            /// @brief Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0)
            {
            }

            /// Returns the position of the vertex.
            const GlobalCoordinate& global(const LocalCoordinate&) const
            {
                return pos_;
            }

            /// Meaningless for the vertex geometry.
            LocalCoordinate local(const GlobalCoordinate&) const
            {
                // return 0 to make the geometry check happy.
                return LocalCoordinate(0.0);
            }

            /// Returns 1 for the vertex geometry.
            double integrationElement(const LocalCoordinate&) const
            {
                return volume();
            }

            /// Using the cube type for vertices.
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(mydimension);
            }

            /// A vertex is defined by a single corner.
            int corners() const
            {
                return 1;
            }

            /// Returns the single corner: the vertex itself.
            GlobalCoordinate corner(int cor) const
            {
                static_cast<void>(cor);
                assert(cor == 0);
                return pos_;
            }

            /// Volume of vertex is arbitrarily set to 1.
            ctype volume() const
            {
                return 1.0;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// This method is meaningless for singular geometries.
            FieldMatrix<ctype, mydimension, coorddimension>
            jacobianTransposed(const LocalCoordinate& /* local */) const
            {

                // Meaningless to call jacobianTransposed() on singular geometries. But we need to make DUNE happy.
                return FieldMatrix<ctype, mydimension, coorddimension>();
            }

            /// This method is meaningless for singular geometries.
            FieldMatrix<ctype, coorddimension, mydimension>
            jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
            {
                // Meaningless to call jacobianInverseTransposed() on singular geometries. But we need to make DUNE happy.
                return FieldMatrix<ctype, coorddimension, mydimension>();
            }

            /// The mapping implemented by this geometry is constant, therefore affine.
            bool affine() const
            {
                return true;
            }

        private:
            GlobalCoordinate pos_;
        };  // class Geometry<0,cdim>




        /// Specialization for 3-dimensional geometries, i.e. cells.
        template <int cdim>
        class Geometry<3, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 3 };
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;

            typedef Dune::Impl::FieldMatrixHelper< double >  MatrixHelperType;

            /// @brief Construct from centroid, volume (1- and 0-moments) and
            ///        corners.
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            /// @param allcorners array of all corner positions in the grid
            /// @param corner_indices array of 8 indices into allcorners array. The
            ///                       indices must be given in lexicographical order
            ///                       by (kji), i.e. i running fastest.
            Geometry(const GlobalCoordinate& pos,
                     ctype vol,
                     const EntityVariable<cpgrid::Geometry<0, 3>, 3>& allcorners,
                     const int* corner_indices)
                : pos_(pos), vol_(vol), allcorners_(allcorners.data()), cor_idx_(corner_indices)
            {
                assert(allcorners_ && corner_indices);
            }

            /// @brief Construct from centroid and volume (1- and
            ///        0-moments).  Note that since corners are not
            ///        given, the geometry provides no mappings, and
            ///        some calls (corner(), global() etc.) will fail.
            ///        This possibly dangerous constructor is
            ///        available for the benefit of
            ///        CpGrid::readSintefLegacyFormat().
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            Geometry(const GlobalCoordinate& pos,
                     ctype vol)
                : pos_(pos), vol_(vol)
            {
            }

            /// Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0), vol_(0.0), allcorners_(0), cor_idx_(0)
            {
            }

            /// Provide a trilinear mapping.
            /// Note that this does not give a proper space-filling
            /// embedding of the cell complex in the general (faulted)
            /// case. We should therefore revisit this at some point.
            /// Map g from (local) reference domain to (global) cell
            GlobalCoordinate global(const LocalCoordinate& local_coord) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");
                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local_coord };
                uvw[0] -= local_coord;
                // Access pattern for uvw matching ordering of corners.
                const int pat[8][3] = { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } };
                GlobalCoordinate xyz(0.0);
                for (int i = 0; i < 8; ++i) {
                    GlobalCoordinate corner_contrib = corner(i);
                    double factor = 1.0;
                    for (int j = 0; j < 3; ++j) {
                        factor *= uvw[pat[i][j]][j];
                    }
                    corner_contrib *= factor;
                    xyz += corner_contrib;
                }
                return xyz;
            }

            /// Mapping from the cell to the reference domain.
            /// May be slow.
            LocalCoordinate local(const GlobalCoordinate& y) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");
                // This code is modified from dune/grid/genericgeometry/mapping.hh
                // \todo: Implement direct computation.
                const ctype epsilon = 1e-12;
                auto refElement = Dune::ReferenceElements<ctype, 3>::cube();
                LocalCoordinate x = refElement.position(0,0);
                LocalCoordinate dx;
                do {
                    // DF^n dx^n = F^n, x^{n+1} -= dx^n
                    JacobianTransposed JT = jacobianTransposed(x);
                    GlobalCoordinate z = global(x);
                    z -= y;
                    MatrixHelperType::template xTRightInvA<3, 3>(JT, z, dx );
                    x -= dx;
                } while (dx.two_norm2() > epsilon*epsilon);
                return x;
            }

            /// Equal to \sqrt{\det{J^T J}} where J is the Jacobian.
            /// J_{ij} = (dg_i/du_j)
            /// where g is the mapping from the reference domain,
            /// and {u_j} are the reference coordinates.
            double integrationElement(const LocalCoordinate& local_coord) const
            {
                JacobianTransposed Jt = jacobianTransposed(local_coord);
                return MatrixHelperType::template sqrtDetAAT<3, 3>(Jt);
            }

            /// Using the cube type for all entities now (cells and vertices),
            /// but we use the singular type for intersections.
            GeometryType type() const
            {
                return Dune::GeometryTypes::cube(mydimension);
            }

            /// The number of corners of this convex polytope.
            /// Returning 8, since we treat all cells as hexahedral.
            int corners() const
            {
                return 8;
            }

            /// @brief Get the i-th of 8 corners of the hexahedral base cell.
            GlobalCoordinate corner(int cor) const
            {
                assert(allcorners_ && cor_idx_);
                return allcorners_[cor_idx_[cor]].center();
            }

            /// Cell volume.
            ctype volume() const
            {
                return vol_;
            }

            void set_volume(ctype volume) {
                vol_ = volume;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// @brief Jacobian transposed.
            /// J^T_{ij} = (dg_j/du_i)
            /// where g is the mapping from the reference domain,
            /// and {u_i} are the reference coordinates.
            /// g = g(u) = (g_1(u), g_2(u), g_3(u)), u=(u_1,u_2,u_3)
            /// g = map from (local) reference domain to global cell
            const JacobianTransposed
            jacobianTransposed(const LocalCoordinate& local_coord) const
            {
                static_assert(mydimension == 3, "");
                static_assert(coorddimension == 3, "");

                // uvw = { (1-u, 1-v, 1-w), (u, v, w) }
                LocalCoordinate uvw[2] = { LocalCoordinate(1.0), local_coord };
                uvw[0] -= local_coord;
                // Access pattern for uvw matching ordering of corners.
                const int pat[8][3] = { { 0, 0, 0 },
                                        { 1, 0, 0 },
                                        { 0, 1, 0 },
                                        { 1, 1, 0 },
                                        { 0, 0, 1 },
                                        { 1, 0, 1 },
                                        { 0, 1, 1 },
                                        { 1, 1, 1 } };
                JacobianTransposed  Jt(0.0);
                for (int i = 0; i < 8; ++i) {
                    for (int deriv = 0; deriv < 3; ++deriv) {
                        // This part contributing to dg/du_{deriv}
                        double factor = 1.0;
                        for (int j = 0; j < 3; ++j) {
                            factor *= (j != deriv) ? uvw[pat[i][j]][j]
                                : (pat[i][j] == 0 ? -1.0 : 1.0);
                        }
                        GlobalCoordinate corner_contrib = corner(i);
                        corner_contrib *= factor;
                        Jt[deriv] += corner_contrib; // using FieldMatrix row access.
                    }
                }
                return Jt;
            }

            /// @brief Inverse of Jacobian transposed. \see jacobianTransposed().
            const JacobianInverseTransposed
            jacobianInverseTransposed(const LocalCoordinate& local_coord) const
            {
                JacobianInverseTransposed Jti = jacobianTransposed(local_coord);
                Jti.invert();
                return Jti;
            }

            /// The mapping implemented by this geometry is not generally affine.
            bool affine() const
            {
                return false;
            }

            /**
             * @brief Refine a single cell with regular intervals.
             * 
             * For each cell to be created, storage must be passed for its corners and the indices. That storage
             * must be externally managed, since the newly created geometry structures only store pointers and do
             * not free them on destruction.
             *
             * @param      cells_per_dim    The number of sub-cells in each direction,
             * @param[out] refined_geom     Geometry Policy for the refined geometries. Those will be added there.
             * @param[out] indices_storage  A vector of mutable references to storage for the indices of each new cell.
             * @return A vector with the created cells.
             * @todo We do not need to return anything here.
             */
            std::vector<Geometry<3, cdim>> refine(const std::array<int, 3>& cells_per_dim,
                                                  DefaultGeometryPolicy& all_geom,
                                                  std::vector<std::array<int, 8>>&  global_refined_cell8corners_indices_storage)
            // the 2nd and 3rd arguments are 'empty' and we populate them once we refine the global cell.
            // @todo QUESTION Why do we need to provide them as arguments? Shall we only return them instead?
            {  
                // Below are basically std::vector of Geometry. Use resize(), reserve(), push_back(), etc.
                EntityVariable<cpgrid::Geometry<0, 3>, 3>& global_refined_corners = all_geom.geomVector(std::integral_constant<int, 3>()); // was called corner_storage before
                EntityVariable<cpgrid::Geometry<2, 3>, 1>& global_refined_faces = all_geom.geomVector(std::integral_constant<int, 1>()); // Missed by Peter, we need to add the faces
                EntityVariable<cpgrid::Geometry<3, 3>, 0>& global_refined_cells = all_geom.geomVector(std::integral_constant<int, 0>()); // Put the refined cells here.
                
                // @todo Do we need "indices_storage" now called "global_refined_cell8corners_indices_storage"?
                // @todo DOUBLE CHECK COMMENTS. REDO THEM IF THEY ARE NOT CLEAR AFTER GROUPING FACE/CELL CODE.
                // @todo DOUBLE CHECK CONSISTANCY OF NOTATION (MANY CHANGES IN NAME VARIABLES)
                // @todo REMOVE CONTAINERS THAT WE DO NOT NEED

                /// --- CORNERS ---
                /// GLOBAL REFINED CORNERS 
                /// The strategy is to compute the local refined corners
                /// of the unit/reference cube, and then apply the map global().
                //
                // Refine corners of the (local) unit/reference cube.
                // Create a vector to store the refined corners of the unit/reference cube
                // such that each entry is one of the (local) refined corners.
                std::vector<std::array<double,3>> local_refined_corners;
                // Determine the size of the vector containing the (local) refined corners of the unit/reference cube. 
                local_refined_corners.resize((cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1));
                // Vector to store the indices of the (local) refined corners of the unit/reference cube.
                // Notice that we use the SAME INDEX for a local refined corners and its global refined version.
                // So we call the vector containing such indices "global_refined_corner_indices".
                std::vector<int> global_refined_corner_indices;
                // The size of "global_refined_corner_indices" is
                // (cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1).
                // The nummbering starts at the botton, so k=0 (z-axis), and j=0 (y-axis), i=0 (x-axis).
                // Then, increasing i, followed by increasing j, and finally, increasing k.
                // 'Right[increasing i]-Back[incresing j]-Up[increasing k]'.
                for (int k = 0; k < cells_per_dim[2] + 1; k++) {
                    for (int j = 0; j < cells_per_dim[1] + 1; j++) {
                        for (int i = 0; i < cells_per_dim[0] + 1; i++) {
                            // Compute the index of each (local/global) refined corner "kji".
                            // Recall that a local refined corner and its associated global version
                            // have the SAME INDEX.
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // Incorporate the index in "global_refined_corner_indices".
                            global_refined_corner_indices.push_back(kji_idx);
                            // Change 'int' type to 'double' for k,j,i. 
                            double kd = k;
                            double jd = j;
                            double id = i;
                            // Compute the 3 (local) coordinates of the "kji" refined corner of the unit/reference cube. 
                            local_refined_corners[kji_idx][0] = id/cells_per_dim[0];
                            local_refined_corners[kji_idx][1] = jd/cells_per_dim[1];
                            local_refined_corners[kji_idx][2] = kd/cells_per_dim[2];
                        } // end i-for-loop
                    } // end j-for-loop
                } // end k-for-loop  
                //
                // Get the global refined corners from the local refined corners,
                // applying global().
                for (const auto& corner : local_refined_corners) { 
                    global_refined_corners.push_back(Geometry<0, 3>(this->global(corner)));
                } // end GLOBAL REFINED CORNERS
                /// --- END CORNERS ---
                
                /// --- FACES ---
                // GLOBAL REFINED FACES
                // We want to populate "global_refined_faces".
                // The size of "global_refined_faces" must be
                // (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))   (botton and top faces)
                //  + (cells_per_dim[0]*cells_per_dim[2]*(cells_per_dim[1]+1)) (front and back faces)
                //  + (cells_per_dim[1]*cells_per_dim[2]*(cells_per_dim[0]+1))  (left and right faces).
                // Determine the size of "global_refined_faces".
                int global_refined_faces_size = (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) // 'botton-top faces'
                    + (cells_per_dim[0]*cells_per_dim[2]*(cells_per_dim[1]+1)) // 'front-back faces'
                    + (cells_per_dim[1]*cells_per_dim[2]*(cells_per_dim[0]+1));  // 'left-right faces'
                global_refined_faces.resize(global_refined_faces_size);
                // 
                // To understand the size, count the 'horizontal(botton-top)' refined faces on the botton (k=0),
                // that is cells_per_dim[0]*cells_per_dim[1]. Now, varying k between 0 and cell_per_dim[2] we
                // get the third factor of the first term. Similarly with the second and fisrt terms where we
                // count the 'vertical'-faces front-back(2nd coordinate constant in each face)
                // and left-right(1st coordinate constant constant in each face) respectively.
                //
                // To create a face as a Geometry<2,3> type object we need its centroid and its volume.
                // CENTROIDS of the faces of the refined global cells.
                // We store the centroids (and later on the faces) in the following order:
                // - botton-top faces
                // - front-back faces
                // - left-right faces
                // Vector to store the centroids of the faces of the global refined cells.
                std::vector<std::array<double,3>> global_refined_face_centroids;
                // Determine the size of "global_refined_face_centroids".
                global_refined_face_centroids.resize(global_refined_faces_size);
                // Vector to store refined faces indices (one index per face = one index per centroid).
                // We call it "global_refined_face_indices" since the INDEX is the SAME for a local
                // face/centroid and its global version(s).
                std::vector<int> global_refined_face_indices;
                // The size of "global_refined_face_indices" (it's the same as the size of "global_refined_faces").
                // One face has one centroid.(We do not 'resize' this vector, it will be 'pushed back' later on).
                //
                // Container to store, in each entry, the 4 indices of the 4 corners
                // of each global refined face. We call the container
                // "global_refined_face4corners_indices_storage" (same size as "global_refined_faces")
                // and, later on, we call each entry "global_refined_face4corners_indices".
                std::vector<std::array<int,4>> global_refined_face4corners_indices_storage;
                // Determine the size of "global_refined_face4corners_indices_storage".
                global_refined_face4corners_indices_storage.resize(global_refined_faces_size);
                //
                // Container to store, in each entry, the {edge_indix[0], edge_index[1]}
                // for each edge of each refined face. We call the container
                // "global_refined_face4edges_indices_storage" (same size as "global_refined_faces")
                // and, later on, we call each entry "global_refined_face4edges_indices".
                std::vector<std::array<std::array<int,2>,4>> global_refined_face4edges_indices_storage;
                // Determine the size of the vector "global_refined_face4edges_indices_storage".
                global_refined_face4edges_indices_storage.resize(global_refined_faces_size);
                //
                // REFINED FACE AREAS
                // To compute the area of each face, we divide it in 4 triangles,
                // compute the area of those with "simplex_volume()", where the arguments
                // are the 3 corners of each triangle. Then, sum them up to get the area
                // of the global refined face.
                //
                // Vector to store all the areas of all the global refined faces.
                std::vector<double> global_refined_face_areas;
                // Determine the size of "global_refined_face_areas"
                // (same as total amount of refined faces and their centroids).
                global_refined_face_areas.resize(global_refined_faces_size);
                //
                // For each face, we construct 4 triangles with
                // (1) centroid of the face, [available in "global_refined_face_centroids"]
                // (2) one of the edges of the face.
                //
                // A triangle has 3 edges. Once we choose a face to base a triangle on,
                // we choose an edge of that face as one of the edges of the triangle.
                // The other 2 edges are fixed, since the centroid of the face they are
                // based on is one of the corners of the triangle. That's why to identify
                // a triangle we only need two things:
                // (1) the face it's based on and
                // (2) one of the four edges of that face.
                //
                // For each face, we need
                // 1. index of the global refined face
                //    [available in "global_refined_face_indices"]
                //    [needed to access centroid of the global refined face in "global_refined_face_indices"]
                //    [needed to access indices of the 4 edges of the face in "global_refined_face4edges_indices_storage"]
                // 2. centroid of the face (common corner of the 4 triangles based on that face).
                //    [available in "global_refined_face_centroids"]
                // 3. container of 4 entries, each entry consists in the 2 indices defining
                //    one edge of the face.
                //    [available in "global_refined_face4edges_indices_storage"].
                //
                // Populate
                // "global_refined_face_centroids".
                // "global_refined_face_areas"
                //
                // Auxiliary 'popoluations'
                // "global_refined_face_indices"
                // "global_refined_face4edges_indices_storage"  @todo PROBABLY THIS ONE IS NOT NEEDED! 
                //
                // - BOTTON-TOP faces, horizontal, 3rd coordinate constant in each face.
                for (k = 0; k < cells_per_dim[2] + 1; k++) {
                    for (j = 0; j < cells_per_dim[1]; j++) {
                        for (i = 0; i < cells_per_dim[0]; i++) {
                            // Compute the index of the refined face (same index for its centroid).
                            int global_refined_face_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i;
                            // Add the index to "global_refined_face_indices".
                            global_refined_face_indices.push_back(refined_face_idx);
                            // Change type 'int' by 'double' (only needed for k).
                            double kd = k;
                            // Centroid of the face of the local refined cell.
                            std::array<double, 3> local_refined_face_centroid = {
                                (.5 + i)/cells_per_dim[0], (.5 + j)/cells_per_dim[1], kd/cells_per_dim[2]};
                            // Compute and add the (global) centroid of the face of the global refined cell.
                            global_refined_face_centroids[global_refined_face_idx] = Geometry<0, 3>(this->global(local_refined_face_centroid));
                            //
                            // We construct "global_refined_face4corners_indices"
                            // with the indices of the 4 corners of the refined face.
                            // The actual 3d corners are available (via these indices)
                            // in "global_refined_corners".
                            std::array<int, 4> global_refined_face4corners_indices = {
                                // fake '{0,1,2,3}' for 'botton' faces, fake '{4,5,6,7}' for 'top' faces
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '0'/'4'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake '1'/'5'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i,// fake '2'/'6'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1}; // fake '3'/'7'
                            // We add "global_refined_face4corners_indices" in the corresponding entry
                            // of "global_refined_face4corners_indices_storage".
                            global_refined_face4corners_indices_storage[global_refined_face_idx] = global_refined_face4corners_indices;
                            //
                            // We construct "global_refined_face4edges_indices"
                            // with the {edge_indix[0], edge_index[1]} for each edge of the refined face.
                            std::array<std::array<int,2>,4> global_refined_face4edges_indices = {
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[1]}, // fake '{0,1}'/'{4,5}'
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[2]}, // fake '{0,2}'/'{4,6}'
                                { global_refined_face4corners_indices[1], global_refined_face4corners_indices[3]}, // fake '{1,3}'/'{5,7}'
                                { global_refined_face4corners_indices[2], global_refined_face4corners_indices[3]}}; // fake '{2,3}'/'{6,7}'
                            // We add "global_refined_face4edges_indices" in the corresponding entry of
                            // "global_refined_face4edges_indices_storage".
                            global_refined_face4edges_indices_storage[global_refined_face_idx] = global_refined_face4edges_indices;
                            //
                            // Calculate the area of each face of a global refined cell,
                            // by adding the 4 areas of the triangles partitioning each face.
                            double global_refined_face_area = 0.0;
                            for (int edge = 0; edge < 4; edge++) {
                                // Construction of each triangle on the current face with one
                                // of its edges equal to "edge".
                                const Geometry<0, 3>::GlobalCoordinate trian_corners[3] = {
                                    global_refined_corners[global_refined_face4edges_indices[edge][0]],
                                    global_refined_corners[global_refined_face4edges_indices[edge][1]],
                                    global_refined_face_centroids[global_refined_face_idx]};
                                global_refined_face_area += std::fabs(simplex_volume(trian_corners));
                            } // end edge-for-loop
                            // We add the "global_refined_face_area" to the corresponding entry
                            // of "global_refined_face_areas".
                            // @todo REMOVE VARIABLES WE DO NOT NEED. DO WE NEED TO STORE ALL THE AREAS???
                            //       SAME QUESTION FOR ALL VOLUMES OF ALL REFINED CELLS
                            global_refined_face_areas[global_refined_face_idx] = global_refined_face_area;
                            //
                            // Construct the Geometry<2,3> of the global refined face. 
                            global_refined_faces[global_refined_face_idx] = Geometry<2,cdim>(global_refined_face_centorid,
                                                                                             global_refined_face_area);
                        } // end i-for-loop
                    } // end j-for-loop
                }; // end k-for-loop
                // - FRONT-BACK faces, vertical, 2nd coordinate constant in each face.
                for (j = 0; j < cells_per_dim[1] + 1; j++) {
                    for (k = 0; k < cells_per_dim[2]; k++) {
                        for (i = 0; i < cells_per_dim[0]; i++) {
                            // Compute the index of the refined face (same index for its centroid).
                            int global_refined_face_idx = (j*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i;
                            // Add the index to "global_refined_face_indices".
                            global_refined_face_indices.push_back(global_refined_face_idx);
                            // Change type 'int' by 'double' (only needed for j).
                            double jd = j;
                            // Centroid of the face of the local refined cell.
                            std::array<double, 3> local_refined_face_centroid = {
                                (.5 + i)/cells_per_dim[0], jd/cells_per_dim[1], (.5 + k)/cells_per_dim[2]};
                            // Compute and add the (global) centroid of the face of the global refined cell.
                            global_refined_face_centroids[global_refined_face_idx] = Geometry<0, 3>(this->global(local_refined_face_centroid));
                            //
                            // We construct "global_refined_face4corners_indices"
                            // with the indices of the 4 corners of the refined face.
                            // The actual 3d corners are available (via these indices)
                            // in "global_refined_corners".
                            std::array<int, 4> global_refined_face4corners_indices = {
                                // fake '{0,1,4,5}' for 'front' faces, fake '{2,3,6,7}' for 'back' faces
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '0'/'2'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake '1'/'3'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '4'/'6'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1}; // fake '5'/'7'
                            // We add "global_refined_face4corners_indices" in the corresponding entry
                            // of "global_refined_face4corners_indices_storage".
                            global_refined_face4corners_indices_storage[global_refined_face_idx] = global_refined_face4corners_indices;
                            //
                            // We construct "global_refined_face4edges_indices"
                            // with the {edge_indix[0], edge_index[1]} for each edge of the refined face.
                            std::array<std::array<int,2>,4> global_refined_face4edges_indices = {
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[1]}, // fake '{0,1}'/'{2,3}'
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[2]}, // fake '{0,4}'/'{2,6}'
                                { global_refined_face4corners_indices[1], global_refined_face4corners_indices[3]}, // fake '{1,5}'/'{3,7}'
                                { global_refined_face4corners_indices[2], global_refined_face4corners_indices[3]}}; // fake '{4,5}'/'{6,7}'
                            // We add "global_refined_face4edges_indices" in the corresponding entry of
                            // "global_refined_face4edges_indices_storage".
                            global_refined_face4edges_indices_storage[global_refined_face_idx] = global_refined_face4edges_indices;
                            //
                            // Calculate the area of each face of a global refined cell,
                            // by adding the 4 areas of the triangles partitioning each face.
                            double global_refined_face_area = 0.0;
                            for (int edge = 0; edge < 4; edge++) {
                                // Construction of each triangle on the current face with one
                                // of its edges equal to "edge".
                                const Geometry<0, 3>::GlobalCoordinate trian_corners[3] = {
                                    global_refined_corners[global_refined_face4edges_indices[edge][0]],
                                    global_refined_corners[global_refined_face4edges_indices[edge][1]],
                                    global_refined_face_centroids[global_refined_face_idx]};
                                global_refined_face_area += std::fabs(simplex_volume(trian_corners));
                            } // end edge-for-loop
                            // We add the "global_refined_face_area" to the corresponding entry
                            // of "global_refined_face_areas".
                            // @todo REMOVE VARIABLES WE DO NOT NEED. DO WE NEED TO STORE ALL THE AREAS???
                            //       SAME QUESTION FOR ALL VOLUMES OF ALL REFINED CELLS
                            global_refined_face_areas[global_refined_face_idx] = global_refined_face_area;
                            //
                            // Construct the Geometry<2,3> of the global refined face. 
                            global_refined_faces[global_refined_face_idx] = Geometry<2,cdim>(global_refined_face_centorid,
                                                                                             global_refined_face_area);
                        } // end i-for-loop
                    } // end k-for-loop
                }; // end j-for-loop
                // - LEFT-RIGHT faces, VERTICAL, 1nd coordinate constant in each face.
                for (i = 0; i < cells_per_dim[0] + 1; i++) {
                    for (k = 0; k < cells_per_dim[2]; k++) {
                        for (j = 0; j < cells_per_dim[1]; j++) {
                            // Compute the index of the refined face (same index for its centroid).
                            int global_refined_face_idx = (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j;
                            // Add the index to "global_refined_face_indices".
                            global_refined_face_indices.push_back(global_refined_face_idx);
                            // Change type 'int' by 'double' (only needed for i).
                            double id = i;
                            // Centroid of the face of the local refined cell.
                            std::array<double, 3> local_refined_face_centroid = {
                                id/cells_per_dim[0], (.5+j)/cells_per_dim[1], (.5 + k)/cells_per_dim[2]};
                            // Compute and add the (global) centroid of the face of the global refined cell.
                            global_refined_face_centroids[global_refined_face_idx] = Geometry<0, 3>(this->global(local_refined_face_centroid));
                            //
                            // We construct "global_refined_face4corners_indices"
                            // with the indices of the 4 corners of the refined face.
                            // The actual 3d corners are available (via these indices)
                            // in "global_refined_corners".
                            std::array<int, 4> global_refined_face4corners_indices = {
                                // fake {0,2,4,6} for 'left' faces, fake '{1,3,5,7}' for 'right' faces
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '0'/'1'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i,// fake '2'/'3'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '4'/'5'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i}; // fake '6'/'7'
                            // We add "global_refined_face4corners_indices" in the corresponding entry
                            // of "global_refined_face4corners_indices_storage".
                            global_refined_face4corners_indices_storage[global_refined_face_idx] = global_refined_face4corners_indices;
                            //
                            // We construct "global_refined_face4edges_indices"
                            // with the {edge_indix[0], edge_index[1]} for each edge of the refined face.
                            std::array<std::array<int,2>,4> global_refined_face4edges_indices = {
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[1]}, // fake '{0,2}'/'{1,3}'
                                { global_refined_face4corners_indices[0], global_refined_face4corners_indices[2]}, // fake '{0,4}'/'{1,5}'
                                { global_refined_face4corners_indices[1], global_refined_face4corners_indices[3]}, // fake '{2,6}'/'{3,7}'
                                { global_refined_face4corners_indices[2], global_refined_face4corners_indices[3]}}; // fake '{4,6}'/'{5,7}'
                            // We add "global_refined_face4edges_indices" in the corresponding entry of
                            // "global_refined_face4edges_indices_storage".
                            global_refined_face4edges_indices_storage[global_refined_face_idx] = global_refined_face4edges_indices;
                            //
                            // Calculate the area of each face of a global refined cell,
                            // by adding the 4 areas of the triangles partitioning each face.
                            double global_refined_face_area = 0.0;
                            for (int edge = 0; edge < 4; edge++) {
                                // Construction of each triangle on the current face with one
                                // of its edges equal to "edge".
                                const Geometry<0, 3>::GlobalCoordinate trian_corners[3] = {
                                    global_refined_corners[global_refined_face4edges_indices[edge][0]],
                                    global_refined_corners[global_refined_face4edges_indices[edge][1]],
                                    global_refined_face_centroids[global_refined_face_idx]};
                                global_refined_face_area += std::fabs(simplex_volume(trian_corners));
                            } // end edge-for-loop
                            // We add the "global_refined_face_area" to the corresponding entry
                            // of "global_refined_face_areas".
                            // @todo REMOVE VARIABLES WE DO NOT NEED. DO WE NEED TO STORE ALL THE AREAS???
                            //       SAME QUESTION FOR ALL VOLUMES OF ALL REFINED CELLS
                            global_refined_face_areas[global_refined_face_idx] = global_refined_face_area;
                            //
                            // Construct the Geometry<2,3> of the global refined face. 
                            global_refined_faces[global_refined_face_idx] = Geometry<2,cdim>(global_refined_face_centorid,
                                                                                             global_refined_face_area);
                        } // end j-for-loop
                    } // end k-for-loop
                }; // end i-for-loop
                // end GLOBAL REFINED FACES
                /// --- END FACES ---
                

                /// --- CELLS ---
                /// GLOBAL REFINED CELLS
                //
                // We need to populate the variable "global_refined_cells"
                // "global_refined_cells" should have size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                // To build each global refined cell, we need
                // 1. its global refined center [available in "global_refined_cell_centers"]
                // 2. its volume [available in "global_refined_cell_volumes"]
                // 3. all global refined corners [available in "global_refined_corners"]
                // 4. indices of its 8 corners [available in "global_refined_corner_indices"]
                //
                // Vector to store, in each entry, the 8 indices of the 8 corners
                // of each global refined cell.
                std::vector<std::array<int,8>> global_refined_cell8corners_indices_storage;
                // Determine the size of "global_refined_cell8corners_indices_storage".
                global_refined_cell8corners_indices_storage.resize(); // same size as "global_refined_cells"
                //
                //  CENTERS
                /// GLOBAL REFINED CELL CENTERS
                /// The strategy is to compute the  centers of the refined local
                /// unit/reference cube, and then apply the map global().
                //
                // We can associate each local/global refined center with an index.
                // As local/global refined corners, a local refined center and its global version
                // have the SAME INDEX. Further, since there is one center per cell, this index can be
                // the SAME INDEX of the new (local/global) refined cell 'kji'.
                // So we call "global_refined_cell_indices" the vector
                // where we store the INDICES for cells/centers of the (local/global) refined cells.
                // The numembering follows the rule of starting with index 0 for the refined cell
                // {0,0,0}, ...,{1/cells_per_dim[0], 1/cells_per_dim[1], 1/cells_per_dim[2]},
                // then the indices grow first picking the cells in the x-axis (Right, i), then y-axis (Back, j), and
                // finally, z-axis (Up, k).
                std::vector<int> global_refined_cell_indices; // same indices for (local/global) refined centers!
                // "global_refined_cell_indices"'s size is cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2].
                // We do not resize this vector, it will be 'pushed back' later on.
                //
                // Vector to store the centers of the global refined cells.
                // Get the global refined centers from the local refined centers.
                std::vector<std::array<double,3>> global_refined_cell_centers;
                // Determine the size of "global_refined_cell_centers" (same as total amount of refined cells).
                global_refined_cell_centers.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                //
                // VOLUMES
                 /// VOLUMES OF THE GLOBAL REFINED CELLS
                /// REMARK: Each global refined 'cell' is a hexahedron since it may not be cube-shaped
                /// since its a 'deformation' of unit/reference cube. We use 'hexahedron' to refer
                /// to the global refined cell in the computation of its volume.
                /// 
                /// The strategy is to construct 24 tetrahedra in each hexahedron.
                /// Each tetrahedron is built with
                /// (1) the center of the hexahedron,
                /// (2) the middle point of the face the tetrahedron is based on, and
                /// (3) one of the edges of the face mentioned in 2.
                /// Each face 'supports' 4 tetrahedra, and we have 6 faces per hexahedron, which
                /// gives us the 24 tetrahedra per 'cell' (hexahedron).
                ///
                /// To compute the volume of each tetrahedron, we use "simplex_volume()" with 
                /// the 6 corners of the tetrahedron as arguments. Summing up the 24 volumes,
                /// we get the volumne of the hexahedorn (global refined 'cell').
                ///
                // Vector to store the volumes of the global refined 'cells'.
                std::vector<double> global_refined_cell_volumes;
                // Determine the size of "global_refined_cell_volumes" (same as total amount of refined cells).
                global_refined_cell_volumes.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                //
                // For each (global refined 'cell') hexahedron, to create 24 tetrahedra and their volumes,
                // we introduce
                // Vol1. "hexa_face_0to5_indices" (needed to access face centroids).
                // Vol2. "hexa_face_centroids" (one of the 6 corners of all 4 tetrahedra based on that face).
                // Vol3. "global_refined_cell_center" the center of the global refined 'cell' (hexahedron) (common corner of the 24 tetrahedra).
                // Vol4. "tetra_edge_indices" indices of the 4x6 tetrahedra per 'cell',
                //     grouped by the face they are based on.
                // Then we construct and compute the volume of the 24 tetrahedra with mainly
                // "hexa_face_centroids" (Vol2.), "global_refined_cell_center" (Vol3.), and "tetra_edge_indices" (Vol4.).
                //
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // INDEX of the 'kji' global refined cell. 
                            int global_refined_cell_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // Incorporate the index in "global_refined_cell_indices".
                            global_refined_cell_indices.push_back(global_refined_cell_idx);
                            //
                            // 1. CENTER of the global refined 'kji' cell (Vol3.)
                            std::array<double, 3> global_refined_cell_center; 
                            // Compute the 3 (local) coordinates of the 'kji' refined center of the unit/reference cube. 
                            local_refined_cell_center = {
                                (.5 + i)/cells_per_dim[0], (.5 + j)/cells_per_dim[1],(.5 + k)/cells_per_dim[2]};
                            // Get the center of the global refined cell from the center of the local refined cell.
                            global_refined_cell_center = Geometry<0, 3>(this->global(local_refined_cell_center));
                            // Store it in the corresponding entry of "global_refined_cell_centers".
                            global_refined_cell_centers[global_refined_cell_idx] = global_refined_cell_center;
                            //
                            // 2. Volume of the global refined 'kji' cell
                            double global_refined_cell_volume = 0.0; // (computed below!)
                            // 3. All Global refined corners ("global_refined_corners")
                            // 4. Indices of the 8 corners of the global refined 'kji' cell
                            std::vector<int, 8> global_refined_cell_corners_indices = {
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '0'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake '1'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake '2'
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1, // fake '3'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake '4'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake '5'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake '6'
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1 // fake '7'
                            };
                            // Add this 8 corners to the corresponding entry of
                            // "global_refined_cell8corners_indices_storage"
                            global_refined_cell8corners_indices_storage[global_refined_cell_idx] = global_refined_cell_corners_indices;
                            //
                            /// VOLUME HEXAHEDRON (GLOBAL REFINED 'CELL')
                            ///
                             // Vol1. Indices ('from 0 to 5') of the faces of the hexahedron (needed to access face centroids).
                            std::vector<int, 6> hexa_face_0to5_indices = {
                                // index face '0' botton
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i,
                                // index face '1' front
                                (j*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i,
                                // index face '2' left
                                (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j,
                                // index face '3' right
                                ((i+1)*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j,
                                // index face '4' back
                                ((j+1)*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i,
                                // index face '5' top
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
                                };
                             //
                            // Vol2. Centroids of the faces of the hexahedron.
                            // (one of the 6 corners of all 4 tetrahedra based on that face).
                            std::vector<std::array<double,3>, 6> hexa_face_centroids;
                            for (auto& idx : hexa_face_0to5_indices) {
                                hexa_face_centroids.push_back(global_refined_face_centroids[idx]);
                            }
                            // Indices of the 4 edges of each face of the hexahedron.
                            // A tetrahedron has six edges. Once we choose a face to base a
                            // tetrahedron on, we choose an edge of that face as one of the
                            // edges of the tetrahedron. The other five edges are fixed, since
                            // the center of the hexahedron and the center of the face are fixed too.
                            // That's why to identify a tetrahedron we only need two things:
                            // (1) the face it's based on and
                            // (2) one of the four edges of that face.
                            //
                            // Vol4. Container with indices of the edges of the 4 tetrahedra per face
                            // [according to description above]
                            const int tetra_edge_indices[6][4][2] = {
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[0]]
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[1]],
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[2]],
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[3]],
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[4]],
                                global_refined_face4edges_indices_storage[hexa_face_0to5_indices[5]]};
                            //
                            // Calculate the volume of each hexahedron, by adding
                            // the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
                            for (int face = 0; face < 6; face++) {
                                for (int edge = 0; edge < 4; edge++) {
                                    // Construction of each tetrahedron based on "face" with one
                                    // of its edges equal to "edge".
                                    const Geometry<0, 3>::GlobalCoordinate tetra_corners[4] = {
                                        global_refined_corners[[tetra_edge_indices[face][edge][0]}],  // (see Vol4.)
                                        global_refined_corners[[tetra_edge_indices[face][edge][1]]],  // (see Vol4.)
                                        hexa_face_centroids[face],  // (see Vol2.)
                                        global_refined_cell_center};  // (see Vol3.)
                                     global_refined_cell_volume += std::fabs(simplex_volume(tetra_corners));
                                } // end edge-for-loop
                            } // end face-for-loop
                            // Sum of the 24 volumes to get the volume of the hexahedron, stored in "hexa__volume".
                            //
                            // Add the volume of the hexahedron (global refined 'cell')
                            // to the container with of all volumes of all the refined cells. 
                            global_refined_cell_volumes[global_refined_cell_idx] = global_refined_cell_volume;
                            //
                            // Construct the Geometry of the refined 'kji' cell 
                            global_refined_cells.push_back(Geometry<3,cdim>(global_refined_cell_center, global_refined_cell_volume,
                                                                     global_refined_corners, global_refined_cell_corners_indices));
                        } // end i-for-loop
                    }  // end j-for-loop
                } // end k-for-loop
                //
                // Rescale all volumes if the sum of volume of all the global refined 'cells' does not match the
                // volume of the 'parent cell'.
                // Sum of all the volumes of all the (children) global refined cells.
                double sum_all_global_refined_cell_volumes = 0.0;
                for (auto& volume : global_refined_cell_volumes) {
                    sum_all_global_refined_cell_volumes += volume;
                };
                // Compare the sum of all the volumes of all refined cells with 'parent cell' volume.
                if (std::fabs(sum_all_global_refined_cell_volumes - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / sum_all_global_refined_cell_volumes;
                    for (auto& volume : global_refined_cell_volumes) {
                        volume *= correction;
                    } 
                } // end if-statement
                // end GLOBAL REFINED CELLS
                /// --- END CELLS ---

        private:
            GlobalCoordinate pos_;
            double vol_;
            const cpgrid::Geometry<0, 3>* allcorners_; // For dimension 3 only
            const int* cor_idx_;               // For dimension 3 only
        };





        /// Specialization for 2 dimensional geometries, that is
        /// intersections (since codim 1 entities are not in CpGrid).
        template <int cdim> // GridImp arg never used
        class Geometry<2, cdim>
        {
            static_assert(cdim == 3, "");
        public:
            /// Dimension of underlying grid.
            enum { dimension = 3 };
            /// Dimension of domain space of \see global().
            enum { mydimension = 2 };
            /// Dimension of range space of \see global().
            enum { coorddimension = cdim };
            /// World dimension of underlying grid.
            enum { dimensionworld = 3 };

            /// Coordinate element type.
            typedef double ctype;

            /// Domain type of \see global().
            typedef FieldVector<ctype, mydimension> LocalCoordinate;
            /// Range type of \see global().
            typedef FieldVector<ctype, coorddimension> GlobalCoordinate;

            /// Type of Jacobian matrix.
            typedef FieldMatrix< ctype, coorddimension, mydimension >         Jacobian;
            /// Type of transposed Jacobian matrix.
            typedef FieldMatrix< ctype, mydimension, coorddimension >         JacobianTransposed;
            /// Type of the inverse of the transposed Jacobian matrix
            typedef FieldMatrix< ctype, coorddimension, mydimension >         JacobianInverseTransposed;

            /// @brief Construct from centroid and volume (1- and 0-moments).
            /// @param pos the centroid of the entity
            /// @param vol the volume(area) of the entity
            Geometry(const GlobalCoordinate& pos,
                     ctype vol)
                : pos_(pos), vol_(vol)
            {
            }

            /// Default constructor, giving a non-valid geometry.
            Geometry()
                : pos_(0.0), vol_(0.0)
            {
            }

            /// This method is meaningless for singular geometries.
            const GlobalCoordinate& global(const LocalCoordinate&) const
            {
                OPM_THROW(std::runtime_error, "Geometry::global() meaningless on singular geometry.");
            }

            /// This method is meaningless for singular geometries.
            LocalCoordinate local(const GlobalCoordinate&) const
            {
                OPM_THROW(std::runtime_error, "Geometry::local() meaningless on singular geometry.");
            }

            /// For the singular geometry, we return a constant
            /// integration element equal to the volume.
            double integrationElement(const LocalCoordinate&) const
            {
                return vol_;
            }

            /// We use the singular type (None) for intersections.
            GeometryType type() const
            {
                return Dune::GeometryTypes::none(mydimension);
            }

            /// The number of corners of this convex polytope.
            /// Since this geometry is singular, we have no corners as such.
            int corners() const
            {
                return 0;
            }

            /// This method is meaningless for singular geometries.
            GlobalCoordinate corner(int /* cor */) const
            {
                // Meaningless call to cpgrid::Geometry::corner(int):
                //"singular geometry has no corners.
                // But the DUNE tests assume at least one corner.
                return GlobalCoordinate( 0.0 );
            }

            /// Volume (area, actually) of intersection.
            ctype volume() const
            {
                return vol_;
            }

            /// Returns the centroid of the geometry.
            const GlobalCoordinate& center() const
            {
                return pos_;
            }

            /// This method is meaningless for singular geometries.
            const FieldMatrix<ctype, mydimension, coorddimension>&
            jacobianTransposed(const LocalCoordinate& /* local */) const
            {
                OPM_THROW(std::runtime_error, "Meaningless to call jacobianTransposed() on singular geometries.");
            }

            /// This method is meaningless for singular geometries.
            const FieldMatrix<ctype, coorddimension, mydimension>&
            jacobianInverseTransposed(const LocalCoordinate& /*local*/) const
            {
                OPM_THROW(std::runtime_error, "Meaningless to call jacobianInverseTransposed() on singular geometries.");
            }

            /// Since integrationElement() is constant, returns true.
            bool affine() const
            {
                return true;
            }

        private:
            GlobalCoordinate pos_;
            ctype vol_;
        };





    } // namespace cpgrid

    template< int mydim, int cdim >
    auto referenceElement(const cpgrid::Geometry<mydim,cdim>& geo) -> decltype(referenceElement<double,mydim>(geo.type()))
    {
        return referenceElement<double,mydim>(geo.type());
    }

} // namespace Dune

#endif // OPM_GEOMETRY_HEADER
