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

                // "global_refined_corners" has size (cells_per_dim[0] + 1)*(cells_per_dim[1] + 1)*(cells_per_dim[2] + 1)
                //                                   [without repeating corners]
                // "refined_faces" has size (3 * cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2])
                //                            + (cells_per_dim[0] * cells_per_dim[1])
                //                            + (cells_per_dim[0] * cells_per_dim[2])
                //                            + (cells_per_dim[1] * cells_per_dim[2])   [without repeating faces]
                // "refined_cells" has size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                

                /// GLOBAL REFINED CORNERS
                /// The strategy is to compute the local refined corners
                /// of the reference cube, and then apply the map global()
                
                // Refine corners the (local) unit/reference cube
                // Create a vector to storade the refined corners of the unit/reference cube
                // such that each entry of this vector is one of them.
                std::vector<std::array<double,3>> refined_reference_corners;
                // Determine the right size of the vector containing the refiend corners of the unit/reference cube 
                refined_reference_corners.resize((cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1));
                // Vector to storage the indices of the refined corners of the unit/reference cube
                std::vector<int> refined_ref_corn_idx;
                // The size of refined_ref_cor_idx is (cells_per_dim[0] + 1) *(cells_per_dim[1] + 1) * (cells_per_dim[2] + 1)
                // The nummbering starts at the botton, so k=0 (z-axis), and j=0 (y-axis), i=0 (x-axis).
                // Then, increasing i, followed by increasing j, and finally, increasing k.
                // "Right[increasing i]-Back[incresing j]-Up[increasing k]"
                for (int k = 0; k < cells_per_dim[2] + 1; k++) {
                    for (int j = 0; j < cells_per_dim[1] + 1; j++) {
                        for (int i = 0; i < cells_per_dim[0] + 1; i++) {
                            // change int type to double for k,j,i
                            double kd = k;
                            double jd = j;
                            double id = i;
                            // Compute the index of each refined corner "kji"
                            int kji_idx = (kd*cells_per_dim[0]*cells_per_dim[1]) + (jd*cells_per_dim[0]) +id;
                            refined_ref_corn_idx.push_back(kji_idx);
                            // Compute the 3 (local) coordinates of the "kji" refined corner of the unit/reference cube 
                            refined_reference_corners[kji_idx][0] = i/cells_per_dim[0];
                            refined_reference_corners[kji_idx][1] = j/cells_per_dim[1];
                            refined_reference_corners[kji_idx][2] = k/cells_per_dim[2];
                        } // end i-for-loop
                    } // end j-for-loop
                } // end k-for-loop  

                // Get the global refined corners from the refined reference corners
                for (const auto& corner : refined_reference_corners) { 
                    global_refined_corners.push_back(Geometry<0, 3>(this->global(corner)));
                    // @todo add the correct index to indices!
                }

                
                /// GLOBAL REFINED CENTERS
                /// The strategy is to compute the local refined centers
                /// of the reference cube, and then apply the map global()
                
                // Refine reference centers
                std::vector<std::array<double,3>> refined_reference_centers;
                // must have size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                refined_reference_centers.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                // We can associate each center with an index
                // This index can be the same index we create for the new refined cell kji
                // The numembering follows the rule of starting with
                // index 0 for the refined cell {0,0,0}, ...,{1/cells_per_dim[0], 1/cells_per_dim[1], 1/cells_per_dim[2]}
                // then the indices grow first picking the cells in the x-axis (right, i), then y-axis (back, j), and
                // finally, z-axis (up, k)
                std::vector<int> refined_cells_indices; // same indices for refined centers!
                //refined_cells_indices.resize(cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]);
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // change int type to double for k,j,i
                            double kd = k;
                            double jd = j;
                            double id = i;
                            // Compute the index of each refined corner "kji"
                            int kji_idx = (kd*cells_per_dim[0]*cells_per_dim[1]) + (jd*cells_per_dim[0]) +id;
                            refined_cells_indices.push_back(kji_idx);
                            // Compute the 3 (local) coordinates of the "kji" refined corner of the unit/reference cube 
                            refined_reference_centers[kji_idx][0] = (.5 + id)/cells_per_dim[0];
                            refined_reference_centers[kji_idx][1] = (.5 + jd)/cells_per_dim[1];
                            refined_reference_centers[kji_idx][2] = (.5 + kd)/cells_per_dim[2];
                        } // end i-for-loop
                    } // end j-for-loop
                } // end k-for-loop            
               
                // Get the global refined centers from the refined reference centers
                std::vector<std::array<double,3>> global_refined_centers;
                for (const auto& ref_center : refined_reference_centers) {
                    global_refined_centers.push_back(Geometry<0, 3>(this->global(ref_center)));
                }

                /// GLOBAL REFINED CELLS
                /// The strategy is to ...

                // We need to populate the variable "refined_cells"
                // refined_cells has size cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2]
                // We need FOR EACH REFINED CELL "c":
                // 1. global_refined_center_c  [can be deduced from "global_refined_centers"]
                // 2. volume_c [done? check Peter's code]  @TODO DOUBLE-CHECK VOLUME COMPUTATION OF REFINED CELL
                // 3. global_refined_corners_c [can be deduced from "global_refined_corners"]
                // 4. indices_c indices of its 8 corners [can be deduced from "refined_ref_corn_idx"]
                //
                for (int k = 0; k < cells_per_dim[2]; k++) {
                    for (int j = 0; j < cells_per_dim[1]; j++) {
                        for (int i = 0; i < cells_per_dim[0]; i++) {
                            // Index of the "kji" global refined cell
                            int kji_idx = (k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i;
                            // 1. Center of the global refined "kji" cell
                            std::array<double, 3> global_refined_center = global_refined_centers[kji_idx];
                            // 2. Volume of the global refined "kji" cell
                            double volume; // TO BE DOUBLE-CHECKED!
                            // 3. 8 corners of the global refined "kji" cell
                            std::vector<std::array<double,3>, 8> global_refined_corns = {
                                global_refined_corners[(k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i], //fake {0,0,0}
                                global_refined_corners[(k*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1], //fake {1,0,0}
                                global_refined_corners[(k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i], //fake {0,1,0}
                                global_refined_corners[(k*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1], //fake {1,1,0}
                                global_refined_corners[((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i], //fake {0,0,1}
                                global_refined_corners[((k+1)*cells_per_dim[0]*cells_per_dim[1]) + (j*cells_per_dim[0]) +i+1], //fake {1,0,1}
                                global_refined_corners[((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i], //fake {0,1,1}
                                global_refined_corners[((k+1)*cells_per_dim[0]*cells_per_dim[1]) + ((j+1)*cells_per_dim[0]) +i+1] //fake{1,1,1}
                            };
                            // 4. 8 indices of the 8 corners of the global refined "kji" cell
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
                                                                     global_refined_corns, refined_corns_indices));
                        } // end i-for-loop
                    }  // end j-for-loop
                } // end k-for-loop
                
                    
                        
                //  // [...].push_back(Geometry<3, cdim>(
                //                global_refined_center, volume, global_refined_corners, indices?);
                

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
