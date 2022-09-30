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
                                                  std::vector<std::array<int, 8>>& indices_storage)
            {

                //std::vector<EntityVariable<Geometry<0, 3>, 3>>& corner_storage,
                // Below are basically std::vector of Geometry. Use resize(), reserve(), push_back(), etc.
                EntityVariable<cpgrid::Geometry<0, 3>, 3>& global_refined_corners = all_geom.geomVector(std::integral_constant<int, 3>()); // was called corner_storage before
                EntityVariable<cpgrid::Geometry<2, 3>, 1>& refined_faces = all_geom.geomVector(std::integral_constant<int, 1>()); // Missed by Peter, we need to add the faces
                EntityVariable<cpgrid::Geometry<3, 3>, 0>& refined_cells = all_geom.geomVector(std::integral_constant<int, 0>()); // Put the refined cells here.

                // @todo Maybe use vector::reserve to prevent allocations (or not), and later use push_back to populate.
                //
                // @todo What do we storage in indices_storage?? Is it what now we call refined_cells_indices?
                // @todo When simplifying code grouping stuff according to the kji-for-loops, take into account that
                //       A. "global_refined_faces_indices" (@todo USE THIS NAME FOR CONSISTENCY) can be used for
                //           1. "global_refined_faces" (@todo USE THIS NAME FOR CONSISTANCY)
                //           2. "global_refined_faces_centroids" (@todo USE THIS NAME FOR CONSISTENCY)
                //           3. "global_refined_faces_areas"
                //       B. "global_refined_cells_indices" (@todo USE THIS NAME FOR CONSISTENCY) can be used for
                //           1. "global_refined_cells" (@todo USE THIS NAME FOR CONSISTENCY)
                //           2. "local_refined_cells_centers" (@todo USE THIS NAME FOR CONSISTENCY)
                //           3. "global_refined_cells_centers" (@todo USE THIS NAME FOR CONSISTENCY)
                //           4. "global_refined_cells_volumes" (@todo USE THIS NAME FOR CONSISTENCY)

                
                /// GLOBAL REFINED CORNERS 
                /// The strategy is to compute the local refined corners
                /// of the unit/reference cube, and then apply the map global()
                //
                // Refine corners the (local) unit/reference cube
                // Create a vector to storage the refined corners of the unit/reference cube
                // such that each entry is one of the (local) refined corners.
                std::vector<std::array<double,3>> refined_reference_corners;
                // Determine the right size of the vector containing the refiend corners of the unit/reference cube 
                refined_reference_corners.resize((cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1));
                // Vector to storage the indices of the (local) refined corners of the unit/reference cube
                std::vector<int> refined_ref_corn_idx;
                // The size of "refined_ref_cor_idx" is (cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1)
                // The nummbering starts at the botton, so k=0 (z-axis), and j=0 (y-axis), i=0 (x-axis).
                // Then, increasing i, followed by increasing j, and finally, increasing k.
                // "Right[increasing i]-Back[incresing j]-Up[increasing k]"
                for (int k = 0; k < cells_per_dim[2] + 1; k++) {
                    for (int j = 0; j < cells_per_dim[1] + 1; j++) {
                        for (int i = 0; i < cells_per_dim[0] + 1; i++) {
                            // Compute the index of each (local) refined corner "kji"
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // Incorporate the index in "refined_ref_corn_idx"
                            refined_ref_corn_idx.push_back(kji_idx);
                            // Change int type to double for k,j,i 
                            double kd = k;
                            double jd = j;
                            double id = i;
                            // Compute the 3 (local) coordinates of the "kji" refined corner of the unit/reference cube 
                            refined_reference_corners[kji_idx][0] = id/cells_per_dim[0];
                            refined_reference_corners[kji_idx][1] = jd/cells_per_dim[1];
                            refined_reference_corners[kji_idx][2] = kd/cells_per_dim[2];
                        } // end i-for-loop
                    } // end j-for-loop
                } // end k-for-loop  
                //
                // Get the global refined corners from the (local) refined reference corners
                // applying global().
                for (const auto& corner : refined_reference_corners) { 
                    global_refined_corners.push_back(Geometry<0, 3>(this->global(corner)));
                } // end GLOBAL REFINED CORNERS

                
                /// GLOBAL REFINED CENTERS
                /// The strategy is to compute the  centers of the refined (local)
                /// unit/reference cube, and then apply the map global()
                //
                // Refined (local) reference centers
                std::vector<std::array<double,3>> refined_reference_centers;
                // Must have size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                refined_reference_centers.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                // We can associate each center with an index.
                // This index can be the same index of the new refined cell kji.
                // The numembering follows the rule of starting with
                // index 0 for the refined cell {0,0,0}, ...,{1/cells_per_dim[0], 1/cells_per_dim[1], 1/cells_per_dim[2]}
                // then the indices grow first picking the cells in the x-axis (Right, i), then y-axis (Back, j), and
                // finally, z-axis (Up, k)
                std::vector<int> refined_cells_indices; // same indices for refined centers!
                // "refined_cells_indices"-size is cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // Compute the index of each refined cell(or its center) "kji"
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // Incorporate the index in "refined_cells_indices"
                            refined_cells_indices.push_back(kji_idx);
                            // Compute the 3 (local) coordinates of the "kji" refined center of the unit/reference cube 
                            refined_reference_centers[kji_idx][0] = (.5 + i)/cells_per_dim[0];
                            refined_reference_centers[kji_idx][1] = (.5 + j)/cells_per_dim[1];
                            refined_reference_centers[kji_idx][2] = (.5 + k)/cells_per_dim[2];
                        } // end i-for-loop
                    } // end j-for-loop
                } // end k-for-loop            
                //
                // Get the global refined centers from the refined reference centers
                std::vector<std::array<double,3>> global_refined_centers;
                for (const auto& ref_center : refined_reference_centers) {
                    global_refined_centers.push_back(Geometry<0, 3>(this->global(ref_center)));
                } // end GLOBAL REFINED CENTERS

               
                // GLOBAL REFINED FACES
                // We want to populate "refined_faces"
                // The size of "refined_face" must be
                // (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1))   "botton and top faces"
                //  + (cells_per_dim[0]*cells_per_dim[2]*(cells_per_dim[1]+1)) "front and back faces"
                //  + (cells_per_dim[1]*cells_per_dim[2]*(cells_per_dim[0]+1))  "left and right faces"
                // 
                // To understand the size, count the "horizontal(botton-top)" refined faces on the botton (k=0),
                // that is cells_per_dim[0]*cells_per_dim[1]. Now, varying k between 0 and cell_per_dim[2] we
                // get the third factor of the first term. Similar with the second and fisrt terms where we
                // count the "vertical"-faces front-back(2nd coordinate constant) and left-right(1st coordinate
                // constant) respectively.
                //
                // To create a face as a Geometry<2,3> type object we need its centroid and its volume.
                // CENTROIDS of the faces of the refined global cells.
                // We storage the centroids (and later on the faces) in the following order:
                // - botton-top faces
                // - front-back faces
                // - left-right faces
                // Vector to storage the centroids of the faces of the global refined cells
                std::vector<std::array<double,3>> centroids_refined_faces;
                // Determine the size of "centroid_refined_faces".
                int centroids_refined_faces_size = (cells_per_dim[0]*cells_per_dim[1]*(cells_per_dim[2]+1)) //   "botton and top faces"
                    + (cells_per_dim[0]*cells_per_dim[2]*(cells_per_dim[1]+1)) // "front and back faces"
                    + (cells_per_dim[1]*cells_per_dim[2]*(cells_per_dim[0]+1));  // "left and right faces"
                // Resize the vector "centroids_refined_faces"'s size.
                centroids_refined_faces.resize(centroids_refined_faces_size);
                // Vector to storage refined faces indices (one index per face = one index per centroid)
                std::vector<int> refined_faces_indices;
                // The size of refined_faces_indices (it's the same as the size of centroids_refined_faces),
                // that is, "refined_faces_indices"-size is equal to "centroids_refined_faces_size".
                // Populate "centroids_refined_faces".
                // - BOTTON-TOP faces, horizontal, 3rd coordinate constant (in each face).
                for (k = 0; k < cells_per_dim[2] + 1; k++) {
                    for (j = 0; j < cells_per_dim[1]; j++) {
                        for (i = 0; i < cells_per_dim[0]; i++) {
                            // Compute the index of the refined face (same index for its centroid).
                            int refined_face_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i;
                            // Add the index to refined_faces_indices.
                            refined_faces_indices.push_back(refined_face_idx);
                            // Change type int by double (only needed for k).
                            double kd = k;
                            // Centroid of the face of the (local) refined reference cell.
                            std::array<double, 3> local_centroid = {
                                (.5 + i)/cells_per_dim[0], (.5 + j)/cells_per_dim[1], kd/cells_per_dim[2]};
                            // Add the centroid of the face of the global refined cell.
                            centroids_refined_faces[refined_face_idx] = Geometry<0, 3>(this->global(local_centroid));
                        } // end i-for-loop
                    } // end j-for-loop
                }; // end k-for-loop
                // - FRONT-BACK faces, vertical, 2nd coordinate constant (in each face).
                for (j = 0; j < cells_per_dim[1] + 1; j++) {
                    for (k = 0; k < cells_per_dim[2]; k++) {
                        for (i = 0; i < cells_per_dim[0]; i++) {
                            // Compute the index of the refined face (same index for its centroid).
                            int refined_face_idx = (j*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i;
                            // Add the index to refined_faces_indices.
                            refined_faces_indices.push_back(refined_face_idx);
                            // Change type int by double (only needed for j).
                            double jd = j;
                            // Centroid of the face of the (local) refined reference cell.
                            std::array<double, 3> local_centroid = {
                                (.5 + i)/cells_per_dim[0], jd/cells_per_dim[1], (.5 + k)/cells_per_dim[2]};
                            // Add the centroid of the face of the global refined cell.
                            centroids_refined_faces[refined_face_idx] = Geometry<0, 3>(this->global(local_centroid));
                        } // end i-for-loop
                    } // end k-for-loop
                }; // end j-for-loop
                // - LEFT-RIGHT faces, VERTICAL, 1nd coordinate constant (in each face).
                for (i = 0; i < cells_per_dim[0] + 1; i++) {
                    for (k = 0; k < cells_per_dim[2]; k++) {
                        for (j = 0; j < cells_per_dim[1]; j++) {
                             // Compute the index of the refined face (same index for its centroid).
                            int refined_face_idx = (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j;
                            // Add the index to refined_faces_indices.
                            refined_faces_indices.push_back(refined_face_idx);
                            // Change type int by double (only needed for i).
                            double id = i;
                            // Centroid of the face of the (local) refined reference cell.
                            std::array<double, 3> local_centroid = {
                                id/cells_per_dim[0], (.5+j)/cells_per_dim[1], (.5 + k)/cells_per_dim[2]};
                            // Add the centroid of the face of the global refined cell.
                            centroids_refined_faces[refined_face_idx] = Geometry<0, 3>(this->global(local_centroid));
                        } // end j-for-loop
                    } // end k-for-loop
                }; // end i-for-loop
                //
                // -------------- In progress
                //
                // REFINED FACE VOLUME (AREA)
                // To compute the volume (area) of each face, we divide it in 4 triangles,
                // compute the volume (area) of those with "simplex_volume()", where the arguments
                // are the 4 corners of each triangle. Then, sum them up to get the volume (area)
                // of the global refined face.
                // Check and correct with total volume (area) of 'parent faces'.
                //
                // Vector to storage all the volumes (areas) of each global refined cell.
                std::vector<double> global_refined_faces_areas;
                // Determine the size of "global_refined_faces_areas"
                // (same as total amount of refined faces and their centroids).
                global_refined_faces_areass.resize(centroids_refined_faces_size);
                //
                // For each face, we construct 4 triangles with
                // 1. centroid of the face,
                // 2. one of the edges of the face.
                // We can re-use "tetra_edge_indices" (defined in VOLUME GLOBAL REFINED CELLS)
                // which consists in 6 arrays, each of them containing the indices of the
                // 4 edges of each face of the refined cell.
                // @todo RECICLE CODE - REDUCE AMOUNT OF FOR-LOOPS
                // Calculate the area of each face of a global refined cell,
                // by adding the 4 areas of the triangles partitioning each face
                // (4x6 = 24 triangles per global refined cell).
                //
                // @TODO INCLUDE WHAT WE NEED TO BE ABLE TO DO THE FOR-LOOP
                // THEN CHECK HOW TO SIMPLIFY THE CODE
                // 
                Geometry<3, cdim>::ctype refined_face_area = 0.0;
                for (int face = 0; face < 6; face++) {
                    for (int edge = 0; edge < 4; edge++) {
                        // Construction of each triangle on "face" with one
                        // of its edges equal to "edge".
                        const Geometry<0, 3>::GlobalCoordinate trian_corners[3] = {
                            global_refined_corners[[tetra_edge_indices[face][edge][0]}],
                            global_refined_corners[[tetra_edge_indices[face][edge][1]]],
                            refined_faces_centroids[face]};
                        refined_face_area += std::fabs(simplex_volume(trian_corners));
                    } // end edge-for-loop
                } // end face-for-loop
                //
                // Rescale all areas if the sum of areas of all the
                // (0) botton faces areas
                
                // (1) front faces areas
                
                // (2) left faces areas
                
                // (3) right faces areas
                
                // (4) back faces areas
                
                // (5) top faces areas
                
                // do not coincide with the respectively areas of the 6 faces of the 'parent cell'
                //-----------------------
                // Sum of all the areas of all the (children) global refined faces.
                // 
                /*  double sum_all_volumes_refined_cells = 0.0;
                for (auto& volume : refined_cells_volume) {
                    sum_all_volumes_refined_cells += volume;
                };
                // Compare the sum of all the volumes of all refined cells with 'parent cell' volume.
                if (std::fabs(sum_all_volumes_refined_cells - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / sum_all_volumes_refined_cells;
                    for (auto& volume : global_refined_volumes) {
                        volume *= correction;
                    } 
                    } // end if-statement */
                // end GLOBAL REFINED FACES




                /// VOLUME OF THE GLOBAL REFINED CELLS
                /// REMARK: Each global refined cell is a hexahedron. We use "hexahedron" to refer
                /// to the global refined cell in the computation of its volume.
                /// 
                /// The strategy is to construct 24 tetrahedra in each hexahedron.
                /// Each tetrahedron is built with
                /// 1. the center of the hexahedron,
                /// 2. the middle point of the face the tetrahedron is based on, and
                /// 3. one of the edges of the face mentioned in 2.
                /// Each face 'supports' 4 tetrahedra, and we have 6 faces per hexahedron, which
                /// gives us the 24 tetrahedra per "cell" (hexahedron).
                ///
                /// To compute the volume of each tetrahedron, we use "simplex_volume()" with 
                /// the 6 corners of the tetrahedron as arguments. Summing up the 24 volumes,
                /// we get the volumne of the hexahedorn (global refined cell).
                ///
                // Vector to storage the volume of the refined cells.
                std::vector<double> global_refined_volumes;
                // Determine the size of "refined_cells_volumes" (same as total amount of refined cells)
                global_refined_volumes.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                //
                // For each refined cell (hexahedron), to create 24 tetrahedra and their volumes,
                // we introduce
                // 1. "refined_cell_faces_indices" (needed to access face centroids).
                // 2. "refined_faces_centroids" (one of the 6 corners of all 4 tetrahedra based on that face)
                // 3. "kji_idx" index of the refined cell (needed to access its center).
                // 4. "center_refined_cell" the center of the global refined cell (common corner of the 24 tetrahedra).
                // 5. "refined_cell_corns_indices" indices of the 8 corners of the global refined cell.
                // 6. "refined_cell_faces_indices" indices of the faces of the global refined cell (5. is used here).
                //     (6. is implicitly used to see the edges of each face).
                // 7. "tetra_edge_indices" indices of the 4x6 tetrhedra per cell,
                //     grouped by the face they are based on (5. is used here).
                // Then we construct and compute the volume of the 24 tetrahedra with mainly
                // "refined_faces_centroids" (2.), "center_refined_cell" (4.), and "tetra_edge_indices" (7.).
                for (k = 0; k < cells_per_dim[2]; k++) {
                    for (j = 0; j < cells_per_dim[1]; j++) {
                        for (i = 0; i < cells_per_dim[0]; i++) {
                            // 1. Indices of the faces of the refined cell (hexahedron) (needed to access face centroids).
                            std::vector<int, 6> refined_cell_faces_indices = {
                                // index face "0" botton
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i,
                                // index face "1" front
                                (j*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i,
                                // index face "2" left
                                (i*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j,
                                // index face "3" right
                                ((i+1)*cells_per_dim[1]*cells_per_dim[2]) + (k*cells_per_dim[0]) + j,
                                // index face "4" back
                                ((j+1)*cells_per_dim[0]*cells_per_dim[2]) + (k*cells_per_dim[0]) + i,
                                // index face "5" top
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) + i
                                };
                            // 2. Centroids of the faces of the refined cell (hexahedron).
                            // (one of the 6 corners of all 4 tetrahedra based on that face).
                            std::vector<std::array<double,3>, 6> refined_faces_centroids;
                            for (auto& idx : refined_cell_faces_indices) {
                                refined_faces_centroids.push_back(centroids_refined_faces[idx]);
                            }
                            // 3. Index of the refined cell (hexahedron) (needed to access its center).
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // 4. Center of the refined cell (hexahedron) (common corner of all 24 tetrahedra).
                            std::array<double, 3> center_refined_cell =  global_refined_centers[kji_idx];
                            // 5. Container with the 8 refined corners indices.
                            std::array<int, 8> refined_cell_corns_indices = {
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake "0"
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake "1"
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake "2"
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1, // fake "3"
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake "4"
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake "5"
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake "6"
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1 // fake "7"
                                };
                            // 6. Container with the faces-indices of the refined cell.
                            std::array<std::array<int, 4>, 6> refined_cell_faces_indices = {
                                // fake face0 indices (botton), fake "{0, 1, 2, 3}"
                                { refined_cell_corns_indices[0], refined_cell_corns_indices[1],
                                  refined_cell_corns_indices[2], refined_cell_corns_indices[3]},
                                // fake face1 indices (front), fake "{0, 1, 4, 5}"
                                { refined_cell_corns_indices[0], refined_cell_corns_indices[1],
                                  refined_cell_corns_indices[4], refined_cell_corns_indices[5]},
                                // fake face2 indices (left), fake "{0, 2, 4, 6}"
                                { refined_cell_corns_indices[0], refined_cell_corns_indices[2],
                                  refined_cell_corns_indices[4], refined_cell_corns_indices[6]},
                                // fake face3 indices (right), fake "{1, 3, 5 ,7}"
                                { refined_cell_corns_indices[1], refined_cell_corns_indices[3],
                                  refined_cell_corns_indices[5], refined_cell_corns_indices[7]},
                                // fake face4 indices (back), fake "{2, 3, 6, 7}"
                                { refined_cell_corns_indices[2], refined_cell_corns_indices[3],
                                  refined_cell_corns_indices[6], refined_cell_corns_indices[7]},
                                // fake face5 indices (top), fake "{4, 5, 6, 7}"
                                { refined_cell_corns_indices[4], refined_cell_corns_indices[5],
                                  refined_cell_corns_indices[6], refined_cell_corns_indices[7]}
                            };
                            // A tetrahedron has six edges. Once we choose a face to base a
                            // tetrahedron on, we choose an edge of that face as one of the
                            // edges of the tetrahedron. The other five edges are fixed, since
                            // the center of the hexahedron and the center of the face are fixed too.
                            // That's why to identify a tetrahedron we only need two things:
                            // (1) the face it's based on and
                            // (2) one of the four edges of that face.
                            //
                            // 7. Container with indices of the edges of the 4 tetrahedra per face
                            // [according to description above]
                            const int tetra_edge_indices[6][4][2] = {
                                // tetra-edge-indices FACE "0" BOTTON
                                {{refined_cell_corns_indices[0], refined_cell_corns_indices[1]},
                                 {refined_cell_corns_indices[0], refined_cell_corns_indices[2]},
                                 {refined_cell_corns_indices[1], refined_cell_corns_indices[3]},
                                 {refined_cell_corns_indices[2], refined_cell_corns_indices[3]}},
                                // tetra-edge-indices FACE "1" FRONT
                                {{refined_cell_corns_indices[0], refined_cell_corns_indices[1]},
                                 {refined_cell_corns_indices[0], refined_cell_corns_indices[4]},
                                 {refined_cell_corns_indices[1], refined_cell_corns_indices[5]},
                                 {refined_cell_corns_indices[4], refined_cell_corns_indices[5]}},
                                // tetra-edge-indices FACE "2" LEFT
                                {{refined_cell_corns_indices[0], refined_cell_corns_indices[2]},
                                 {refined_cell_corns_indices[0], refined_cell_corns_indices[4]},
                                 {refined_cell_corns_indices[2], refined_cell_corns_indices[6]},
                                 {refined_cell_corns_indices[4], refined_cell_corns_indices[6]}},
                                // tetra-edge-indices FACE "3" RIGHT
                                {{refined_cell_corns_indices[1], refined_cell_corns_indices[3]},
                                 {refined_cell_corns_indices[1], refined_cell_corns_indices[5]},
                                 {refined_cell_corns_indices[3], refined_cell_corns_indices[7]},
                                 {refined_cell_corns_indices[5], refined_cell_corns_indices[7]}},
                                // tetra-edge-indices FACE "4" BACK
                                {{refined_cell_corns_indices[2], refined_cell_corns_indices[3]},
                                 {refined_cell_corns_indices[2], refined_cell_corns_indices[6]},
                                 {refined_cell_corns_indices[3], refined_cell_corns_indices[7]},
                                 {refined_cell_corns_indices[6], refined_cell_corns_indices[7]}},
                                // tetra-edge-indices FACE "5" TOP
                                {{refined_cell_corns_indices[4], refined_cell_corns_indices[5]},
                                 {refined_cell_corns_indices[4], refined_cell_corns_indices[6]},
                                 {refined_cell_corns_indices[5], refined_cell_corns_indices[7]},
                                 {refined_cell_corns_indices[6], refined_cell_corns_indices[7]}},
                            };
                            ///
                            /// VOLUME GLOBAL REFINED CELL
                            ///
                            // Calculate the volume of each global refined cell (hexahedron),
                            // by adding the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
                            Geometry<3, cdim>::ctype refined_cell_volume = 0.0;
                            for (int face = 0; face < 6; face++) {
                                for (int edge = 0; edge < 4; edge++) {
                                    // Construction of each tetrahedron based on "face" with one
                                    // of its edges equal to "edge".
                                    const Geometry<0, 3>::GlobalCoordinate tetra_corners[4] = {
                                        global_refined_corners[[tetra_edge_indices[face][edge][0]}],
                                        global_refined_corners[[tetra_edge_indices[face][edge][1]]],
                                        refined_faces_centroids[face],
                                        center_refined_cell};
                                    refined_cell_volume += std::fabs(simplex_volume(tetra_corners));
                                } // end edge-for-loop
                            } // end face-for-loop
                            // Sum of the 24 volumes to get the volume of the hexahedron is storaged in "refined_cel_volume"
                            //
                            // Add the volume of the refined cell to the container
                            // of all volumes of all the refined cells. 
                            global_refined_volumes[kji_idx] = refined_cell_volume;
                        } // end i-for-loop   
                    } // end j-for-lopp
                } // end k-for-loop
                //
                // Rescale all volumes if the sum of volume of all the global refined cells does not match the
                // volume of the 'parent cell'.
                // Sum of all the volumes of all the (children) global refined cells.
                double sum_all_volumes_refined_cells = 0.0;
                for (auto& volume : refined_cells_volume) {
                    sum_all_volumes_refined_cells += volume;
                };
                // Compare the sum of all the volumes of all refined cells with 'parent cell' volume.
                if (std::fabs(sum_all_volumes_refined_cells - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / sum_all_volumes_refined_cells;
                    for (auto& volume : global_refined_volumes) {
                        volume *= correction;
                    } 
                } // end if-statement
            // end VOLUME OF THE GLOBAL REFINED CELLS
                
                /// GLOBAL REFINED CELLS
                //
                // We need to populate the variable "refined_cells"
                // "refined_cells" should have size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                // To build each global refined cell, we need
                // 1. its global refined center [available in "global_refined_centers"]
                // 2. its volume [available in "global_refined_volumes"]
                // 3. all global refined corners [available in "global_refined_corners"]
                // 4. indices of its 8 corners [available in "refined_ref_corn_idx"]
                //
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // Index of the "kji" global refined cell (needed to acces its center and volume).
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // 1. Center of the global refined "kji" cell
                            std::array<double, 3> global_refined_center = global_refined_centers[kji_idx];
                            // 2. Volume of the global refined "kji" cell
                            double volume = global_refined_volumnes[kji_idx];
                            // 3. All Global refined corners ("global_refined_corners")
                            // 4. Indices of the 8 corners of the global refined "kji" cell
                            std::vector<int, 8> refined_corns_indices = {
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake 0
                                (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake 1
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake 2
                                (k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1, // fake 3
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i, // fake 4
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1, // fake 5
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i, // fake 6
                                ((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1 // fake 7
                            };
                            // Construct the Geometry of the refined "kji" cell 
                            refined_cells.push_back(Geometry<3,cdim>(global_refined_center, volume,
                                                                     global_refined_corners, refined_corns_indices));
                        } // end i-for-loop
                    }  // end j-for-loop
                } // end k-for-loop
                // end GLOBAL REFINED CELLS

                // ----------------------------------------------------------------------------

                    
                // The center of the parent in local coordinates.
                const Geometry<3, cdim>::LocalCoordinate parent_center(this->local(this->center()));

                // Corners of the parent hexahedron in order, in local coordinates.
                const Geometry<3, cdim>::LocalCoordinate parent_corners[8] = {
                    {0.0, 0.0, 0.0},
                    {1.0, 0.0, 0.0},
                    {0.0, 1.0, 0.0},
                    {1.0, 1.0, 0.0},
                    {0.0, 0.0, 1.0},
                    {1.0, 0.0, 1.0},
                    {0.0, 1.0, 1.0},
                    {1.0, 1.0, 1.0},
                };

                // Indices of the corners of the 6 faces of the hexahedrons.
                const int face_corner_indices[6][4] = {
                    {0, 1, 2, 3},
                    {0, 1, 4, 5},
                    {0, 2, 4, 6},
                    {1, 3, 5, 7},
                    {2, 3, 6, 7},
                    {4, 5, 6, 7},
                };

                // To calculate a refined cell's volume, the hexahedron is
                // divided in 24 tetrahedrons, each of which is defined by the
                // center of the cell, the center of one face, and by one edge
                // of that face. This struct defines that edge for each face,
                // for each of the four possible tetrahedrons that are based on
                // that face.
                // In other words, the center of the cell and the center of each
                // face are fixed. A tetrahedron has six edges. Once we choose a
                // face to base a tetrahedron on, we choose an edge of that face 
                // as one of the edges of the tetrahedron. The other five edges
                // are fixed, since the center of the cell and the center
                // of the face are fixed too. That's why to identify a tetrahedron
                // we only need two things: the face it's based on and one of the
                // four edges of that face.
                const int tetra_edge_indices[6][4][2] = {
                    {{0, 1}, {0, 2}, {1, 3}, {2, 3}},
                    {{0, 1}, {0, 4}, {1, 5}, {4, 5}},
                    {{0, 2}, {0, 4}, {2, 6}, {4, 6}},
                    {{1, 3}, {1, 5}, {3, 7}, {5, 7}},
                    {{2, 3}, {2, 6}, {3, 7}, {6, 7}},
                    {{4, 5}, {4, 6}, {5, 7}, {6, 7}},
                };


                ///    Apart from re-ordering code, there is no substantial change in the rest of the code  [Anto.20.09.2022]///

                std::vector<Geometry<3, cdim>> result;
                result.reserve(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);

                auto pis = indices_storage.begin();

                Geometry<3, cdim>::ctype total_volume = 0.0;
                // Each refined cell has kji values associated with it, where
                // k=0, ..., cells_per_dim[2] (z-axis)
                // j=0, ..., cells_per_dim[1] (y-axis)
                // i=0, ..., cells_per_dim[0] (x-axis)
                // We can think as the parent cell as a collection of horizontal
                // cells_per_dim[2] slices. Each slice has height 1/cells_per_dim[2]
                // and cell_per_dim[1]*cell_per_dim[0] refined cells.
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    Geometry<3, cdim>::LocalCoordinate refined_corners[8];
                    Geometry<3, cdim>::LocalCoordinate refined_center(0.0);

                    refined_center[2] = (parent_center[2] + k) / cells_per_dim[2];
                    // 3rd (local) coordinate of the 8 corners of the refined cell kji
                    for (int h = 0; h < 8; h++) {
                        refined_corners[h][2] = (parent_corners[h][2] + k) / cells_per_dim[2];
                    }
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        // 2nd coordinate of the center of the refined cell kji
                        refined_center[1] = (parent_center[1] + j) / cells_per_dim[1];
                        // 2nd (local) coordinate of the 8 corners of the refined cell kji
                        for (int h = 0; h < 8; h++) {
                            refined_corners[h][1] = (parent_corners[h][1] + j) / cells_per_dim[1];
                        }
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // 1st coordinate of the center of the refined cell kji
                            refined_center[0] = (parent_center[0] + i) / cells_per_dim[0];
                            // 1st coordinate of the 8 corners of the refined cell kji
                            for (int h = 0; h < 8; h++) {
                                refined_corners[h][0] = (parent_corners[h][0] + i) / cells_per_dim[0];
                            }  // end h-for-lopp

                            for (const auto& corner : refined_corners) {
                                // @todo Only push new corners.
                                global_refined_corners.push_back(Geometry<0, 3>(this->global(corner)));
                                // @todo add the correct index to indices!
                            }

                            // The indices must match the order of the constant
                            // arrays containing unit corners, face indices, and
                            // tetrahedron edge indices. Do not reorder.
                            // Recall that pis = indices_storage.begin()
                            auto& indices = *pis++;
                            // @todo use the correct indices for the corner lookup!
                            // global_refined_corner[indices[0]] has to be the Geometry of the first corner of the cell!
                            indices = {0, 1, 2, 3, 4, 5, 6, 7};

                            // Get the center of the cell.
                            const Geometry<3, cdim>::GlobalCoordinate global_refined_center(
                                this->global(refined_center));

                            // Get the 8 corners of the global refined cell.
                            const auto& hex_corners = global_refined_corners.data();
                            Geometry<0, 3>::GlobalCoordinate face_centers[6];
                            // Calculate the centers of the 6 faces.
                            for (int f = 0; f < 6; f++) {
                                face_centers[f] = hex_corners[face_corner_indices[f][0]].center();
                                face_centers[f] += hex_corners[face_corner_indices[f][1]].center();
                                face_centers[f] += hex_corners[face_corner_indices[f][2]].center();
                                face_centers[f] += hex_corners[face_corner_indices[f][3]].center();
                                face_centers[f] /= 4;
                            }

                            // @todo Calculate face volume and add the face geometries to refined_faces!

                            // Calculate the volume of the global refined cell by adding the 4 tetrahedrons at each face.
                            Geometry<3, cdim>::ctype volume = 0.0;
                            for (int f = 0; f < 6; f++) {
                                for (int e = 0; e < 4; e++) {
                                    const Geometry<0, 3>::GlobalCoordinate tetra_corners[4]
                                        = {hex_corners[tetra_edge_indices[f][e][0]].center(),
                                           hex_corners[tetra_edge_indices[f][e][1]].center(),
                                           face_centers[f],
                                           global_refined_center};
                                    volume += std::fabs(simplex_volume(tetra_corners));
                                }
                            }
                            total_volume += volume;

                            // @todo the geometries should go to refined_cells instead
                            result.push_back(Geometry<3, cdim>(
                                global_refined_center, volume, global_refined_corners, indices.data()));
                        } // end i-for-loop
                     } // end j-for-loop
                } // end k-for-loop

                // Rescale all volumes if the sum of volumes does not match the parent.
                if (std::fabs(total_volume - this->volume())
                    > std::numeric_limits<Geometry<3, cdim>::ctype>::epsilon()) {
                    Geometry<3, cdim>::ctype correction = this->volume() / total_volume;
                    for (auto& r : result) {
                        r.set_volume(r.volume() * correction);
                    }
                }

                return result;
            }

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
