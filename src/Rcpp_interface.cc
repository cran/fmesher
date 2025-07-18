/*
 *  Copyright Finn Lindgren (2010-2025)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifdef FMESHER_WITH_R

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

// Order must be RcppFmesher, Rcpp to allow the Rcpp classes
// to find the fmesher types
#include "RcppFmesher.h"
#include "Rcpp.h"

#include "fmesher.h"
#include "fmesher_helpers.h"

using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::string;

using fmesh::constrListT;
using fmesh::constrMetaT;
using fmesh::constrT;
using fmesh::Dart;
using fmesh::DartList;
using fmesh::DartPair;
using fmesh::Int3;
using fmesh::Int3Raw;
using fmesh::Int4;
using fmesh::Int4Raw;
using fmesh::Matrix;
using fmesh::Matrix3double;
using fmesh::Matrix3int;
using fmesh::Matrix4double;
using fmesh::Matrix4int;
using fmesh::Matrix1int;
using fmesh::MatrixC;
using fmesh::Mesh;
using fmesh::Mesh3;
using fmesh::Dart3;
using fmesh::MeshC;
using fmesh::Point;
using fmesh::PointRaw;
using fmesh::SparseMatrix;
using fmesh::TriangleLocator;
using fmesh::Vector3;
using fmesh::vertexListT;

#ifdef FMESHER_WITH_EIGEN
template <class T> using EigenM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <class T> using EigenM1 = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <class T> using EigenSM = Eigen::SparseMatrix<T>;
template <class T> using EigenMM = Eigen::Map<EigenM<T>>;
template <class T> using EigenMM1 = Eigen::Map<EigenM1<T>>;
template <class T> using EigenMSM = Eigen::Map<EigenSM<T>>;
#endif

// const bool useVT = true;
// const bool useTTi = true;

template <typename T>
bool Rcpp_is_element(const Rcpp::List& list, std::string name) {
  if (!list.containsElementNamed(name.c_str()))
    return false;

  if (Rf_isNull(list[name.c_str()]))
    return false;

  return Rcpp::is<T>(list[name.c_str()]);
}


class Options {
public:
  bool delaunay;
  double cutoff;
  double sphere_tolerance;
  int cet_sides;
  double cet_margin;
  double rcdt_min_angle;
  double rcdt_max_edge;
  Matrix<double> quality;
  int rcdt_max_n0;
  int rcdt_max_n1;
  bool rcdt;

public:
  Options(Rcpp::List& options, size_t rows) :
  delaunay(true),
  cutoff(1.0e-12),
  sphere_tolerance(1.0e-7),
  cet_sides(8),
  cet_margin(-0.1),
  rcdt_min_angle(21),
  rcdt_max_edge(-1.0),
  quality(1),
  rcdt_max_n0(-1),
  rcdt_max_n1(-1),
  rcdt(true) {
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "cutoff"))
      cutoff = options["cutoff"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "sphere_tolerance"))
      sphere_tolerance = options["sphere_tolerance"];
    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "cet_sides"))
      cet_sides = options["cet_sides"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "cet_margin"))
      cet_margin = options["cet_margin"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "rcdt_min_angle"))
      rcdt_min_angle = options["rcdt_min_angle"];
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "rcdt_max_edge"))
      rcdt_max_edge = options["rcdt_max_edge"];
    if (Rcpp_is_element<Rcpp::LogicalVector>(options, "rcdt"))
      rcdt = options["rcdt"];
    if (Rcpp_is_element<Rcpp::LogicalVector>(options, "delaunay"))
      delaunay = options["delaunay"];

    /* Construct quality info */
    if (Rcpp_is_element<Rcpp::NumericVector>(options, "quality")) {
      quality = Rcpp::as<Rcpp::NumericVector>(options["quality"]);
      for (size_t r = quality.rows(); r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
      quality.rows(rows); /* Make sure we have the right number of rows */
    } else {
      quality.rows(rows);
      for (size_t r = 0; r < rows; r++)
        quality(r, 0) = rcdt_max_edge;
    }

    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "rcdt_max_n0"))
      rcdt_max_n0 = options["rcdt_max_n0"];
    if (Rcpp_is_element<Rcpp::IntegerVector>(options, "rcdt_max_n1"))
      rcdt_max_n1 = options["rcdt_max_n1"];
  };

};



Mesh Rcpp_import_mesh(Rcpp::NumericMatrix mesh_loc,
                      Rcpp::IntegerMatrix mesh_tv,
                      MatrixC & matrices,
                      Rcpp::List options) {
  const bool useVT = true;
  const bool useTTi = true;

  matrices.attach("mesh_loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(mesh_loc))));
  FMLOG("'mesh_loc' points imported." << std::endl);
  matrices.attach("mesh_tv", std::make_unique<Matrix<int>>(mesh_tv));
  FMLOG("'mesh_tv' points imported." << std::endl);

  Matrix<double>& iS0 = matrices.DD("mesh_loc");
  Matrix<int>& TV0 = matrices.DI("mesh_tv");

  /* Initialise mesh structure */
  Mesh M(Mesh::Mtype::Plane, 0, useVT, useTTi);
  //  if ((iS0.rows() > 0) && (iS0.cols() < 2)) {
  //    /* 1D data. Not implemented */
  //    FMLOG("1D data not implemented." << std::endl);
  //    return Rcpp::List();
  //  }

  if (iS0.rows() > 0) {
    //    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(iS0);
  }

  Options the_options(options, iS0.rows());

  //  double sphere_tolerance = 1e-10;
  (void)M.auto_type(the_options.sphere_tolerance);

  M.TV_set(TV0);

  return M;
}


Mesh3 Rcpp_import_mesh3d(Rcpp::NumericMatrix mesh_loc,
                         Rcpp::IntegerMatrix mesh_tv,
                         MatrixC & matrices,
                         Rcpp::List options) {
  const bool useVT = true;
  const bool useTTi = true;

  matrices.attach("mesh_loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(mesh_loc))));
  FMLOG("'mesh_loc' points imported." << std::endl);
  matrices.attach("mesh_tv", std::make_unique<Matrix<int>>(mesh_tv));
  FMLOG("'mesh_tv' points imported." << std::endl);

  Matrix<double>& iS0 = matrices.DD("mesh_loc");
  Matrix<int>& TV0 = matrices.DI("mesh_tv");

  /* Initialise mesh structure */
  Mesh3 M(Mesh3::Mtype::Plane, 0, useVT, useTTi);

  if (iS0.rows() > 0) {
    //    Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
    M.S_append(iS0);
  }

  Options the_options(options, iS0.rows());

  /* Make sure all tetra have positive volume, i.e. are oriented appropriately,
   * with the first three vertices forming a right-handed coordinate system as
   *  seen from the fourth vertex.
   */
  for (size_t t = 0; t < TV0.rows(); t++) {
    if (M.tetraVolume(M.S(TV0[t][0]), M.S(TV0[t][1]), M.S(TV0[t][2]), M.S(TV0[t][3])) < 0.0) {
      std::swap(TV0(t)[0], TV0(t)[1]);
    }
  }
  M.TV_set(TV0);

  return M;
}







#include "qtool.h"

//' @title Compute sparse matrix inverse
//'
//' @description
//' Requires RcppEigen which is not compiled in by default. Enable with
//' `PKG_CPPFLAGS=-DFMESHER_WITH_EIGEN` in `src/Makevars` and add `RcppEigen`
//' to the `DESCRIPTION` `LinkingTo` field.
//'
//' @param AA A sparse matrix
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List fmesher_qinv(SEXP AA) {
#ifdef FMESHER_WITH_EIGEN
  const EigenMSM<double> A(Rcpp::as<EigenMSM<double>>(AA));

  QTool<double> Q;
  Q.Q(A);

  Rcpp::List ret;
  ret["Qinv"] = Q.S();
  return ret;
#else
  Rcpp::stop("Unsupported method fmesher_qinv; fmesher was built without FMESHER_WITH_EIGEN");
  Rcpp::List ret;
  return ret;
#endif
}



//' @title Globe points
//'
//' @description
//' Create points on a globe
//'
//' @param globe integer; the number of edge subdivision segments, 1 or higher.
//' @returns A matrix of points on a unit radius globe
//' @examples
//' fmesher_globe_points(1)
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix fmesher_globe_points(Rcpp::IntegerVector globe) {
  MatrixC matrices;

  int num = globe[0];
  if (num < 1) {
    num = 1;
  }

  matrices.attach(
    ".globe",
    fmesh::make_globe_points(num, 1.0));
  FMLOG("globe points constructed." << std::endl);

  return Rcpp::wrap(matrices.DD(".globe"));
}




//' @title Refined Constrained Delaunay Triangulation
//'
//' @description
//' (...)
//'
//' @param options list of triangulation options
//' @param loc numeric matrix; initial points to include
//' @param tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param boundary 2-column integer matrix with 0-based vertex indices for each
//' boundary edge constraint
//' @param interior 2-column integer matrix with 0-based vertex indices for each
//' interior edge constraint
//' @param boundary_grp integer vector with group labels
//' @param interior_grp integer vector with group labels
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' @returns A list of information objects for a generated triangulation
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_rcdt(Rcpp::List options,
                           Rcpp::NumericMatrix loc,
                           Rcpp::Nullable<Rcpp::IntegerMatrix> tv = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerMatrix> boundary = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerMatrix> interior = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerVector> boundary_grp = R_NilValue,
                           Rcpp::Nullable<Rcpp::IntegerVector> interior_grp = R_NilValue) {
     const bool useVT = true;
     const bool useTTi = true;

     MatrixC matrices;

     matrices.attach("loc", std::make_unique<Matrix<double>>(loc));
     FMLOG("'loc' points imported." << std::endl);

     Matrix<double>& iS0 = matrices.DD("loc");
     Matrix<int>* TV0 = NULL;
     if (!tv.isNull()) {
       matrices.attach("tv0",
                       std::make_unique<Matrix<int>>(Rcpp::as<Rcpp::IntegerMatrix>(tv)));
       FMLOG("'tv0' points imported." << std::endl);
       TV0 = &matrices.DI("tv0");
     }

     Options rcdt_options(options, iS0.rows());
     FMLOG("rcdt_options parsed" << std::endl);

     /* Prepare boundary/interior edges */
     matrices.attach("boundary", std::make_unique<Matrix<int>>(2));
     matrices.attach("interior", std::make_unique<Matrix<int>>(2));
     matrices.attach("boundary_grp", std::make_unique<Matrix<int>>(1));
     matrices.attach("interior_grp", std::make_unique<Matrix<int>>(1));
     if (!boundary.isNull()) {
       matrices.DI("boundary") = Rcpp::as<Rcpp::IntegerMatrix>(boundary);
     }
     if (!interior.isNull()) {
       matrices.DI("interior") = Rcpp::as<Rcpp::IntegerMatrix>(interior);
     }
     if (!boundary_grp.isNull()) {
       matrices.DI("boundary_grp") = Rcpp::as<Rcpp::IntegerVector>(boundary_grp);
     } else {
       matrices.DI("boundary_grp")(0, 0) = 1;
     }
     if (!interior_grp.isNull()) {
       matrices.DI("interior_grp") = Rcpp::as<Rcpp::IntegerVector>(interior_grp);
     } else {
       matrices.DI("interior_grp")(0, 0) = 1;
     }

     constrListT cdt_boundary;
     constrListT cdt_interior;
     if (!boundary.isNull()) {
       prepare_cdt_input(matrices.DI("boundary"),
                         matrices.DI("boundary_grp"),
                         cdt_boundary);
     }
     if (!interior.isNull()) {
       prepare_cdt_input(matrices.DI("interior"),
                         matrices.DI("interior_grp"),
                         cdt_interior);
     }

     /* Prepare to filter out points at distance not greater than 'cutoff' */
     matrices.attach("idx", std::make_unique<Matrix<int>>(iS0.rows(), 1));
     matrices.output("idx");
     Matrix<int> &idx = matrices.DI("idx").clear();

     filter_locations(iS0, idx, rcdt_options.cutoff);

     /* Remap vertex input references */
     if (TV0) {
       remap_vertex_indices(idx, *TV0);
     }
     remap_vertex_indices(idx, cdt_boundary);
     remap_vertex_indices(idx, cdt_interior);

     /* Initialise mesh structure */
     Mesh M(Mesh::Mtype::Plane, 0, useVT, useTTi);
     if ((iS0.rows() > 0) && (iS0.cols() < 2)) {
       /* 1D data. Not implemented */
       FMLOG("1D data not implemented." << std::endl);
       return Rcpp::List();
     }

     if (iS0.rows() > 0) {
       FMLOG("Append S0." << std::endl);
       Matrix3double S0(iS0); /* Make sure we have a Nx3 matrix. */
       M.S_append(S0);
     }

     FMLOG("Auto-detect manifold type." << std::endl);
     //  double sphere_tolerance = 1e-10;
     (void)M.auto_type(rcdt_options.sphere_tolerance);

     if (TV0) {
       FMLOG("Set TV0." << std::endl);
       M.TV_set(*TV0);
     }

     FMLOG("Attach 's'." << std::endl);
     matrices.attach(string("s"), &M.S());
     FMLOG("Attach 'tv'." << std::endl);
     matrices.attach("tv", &M.TV());
     FMLOG("Set output of 's' and 'tv'." << std::endl);
     matrices.output("s").output("tv");

     FMLOG("Create MeshC helper." << std::endl);
     MeshC MC(&M);
     MC.setOptions(MC.getOptions() | MeshC::Option_offcenter_steiner);

     FMLOG("rcdt.options.delaunay = " << rcdt_options.delaunay << std::endl);
     FMLOG("MC state = " << MC.getState() << std::endl);
     if (!rcdt_options.delaunay ||
         ((M.type() != Mesh::Mtype::Plane) &&
         (M.type() != Mesh::Mtype::Sphere))) {
       FMLOG("MC state = " << MC.getState() << std::endl);
       if (rcdt_options.delaunay && (M.nT() == 0)) {
         FMLOG("MC state = " << MC.getState() << std::endl);
         FMLOG(
           "Points not in the plane or on a sphere, and triangulation empty."
           << std::endl);
       }
       FMLOG("MC state = " << MC.getState() << std::endl);
       if (rcdt_options.delaunay) {
         /* Remove everything outside the boundary segments, if any. */
         FMLOG("Prune exterior." << std::endl);
         MC.PruneExterior();
         FMLOG("MC state = " << MC.getState() << std::endl);
       }
       FMLOG("Invalidate unused vertex indices." << std::endl);
       invalidate_unused_vertex_indices(M, idx);
       FMLOG("MC state = " << MC.getState() << std::endl);

       /* Add the boundary segments, if any. */
       MC.make_boundary_segments();

       /* Nothing more to do here.  Cannot refine non R2/S2 meshes. */
     } else {
       /* If we don't already have a triangulation, we must create one. */
       if (M.nT() == 0) {
         FMLOG("Create covering triangulation." << std::endl);
         FMLOG("cet_sides = " << rcdt_options.cet_sides << std::endl);
         FMLOG("cet_margin = " << rcdt_options.cet_margin << std::endl);
         if (!MC.CET(rcdt_options.cet_sides, rcdt_options.cet_margin)) {
           FMLOG_("CET creation failed, exiting." << std::endl);
           return Rcpp::wrap(matrices);
         }
       }

       /* It is more robust to add the constraints before the rest of the
        nodes are added.  This allows points to fall onto constraint
        segments, subdividing them as needed. */
       if (cdt_boundary.size() > 0)
         MC.CDTBoundary(cdt_boundary);
       if (cdt_interior.size() > 0)
         MC.CDTInterior(cdt_interior);

       /* Add the rest of the nodes. */
       vertexListT vertices;
       for (size_t v = 0; v < iS0.rows(); v++)
         vertices.push_back(v);

       FMLOG("DT(vertices)" << std::endl);
       MC.DT(vertices);

       /* Remove everything outside the boundary segments, if any. */
       FMLOG("Prune exterior." << std::endl);
       MC.PruneExterior();
       FMLOG("Invalidate unused vertex indices." << std::endl);
       invalidate_unused_vertex_indices(M, idx);

       if ((rcdt_options.rcdt) &&
           (rcdt_options.rcdt_max_edge > 0)) {
         /* Calculate the RCDT: */
         FMLOG("Construct refinement." << std::endl);
         MC.RCDT(rcdt_options.rcdt_min_angle,
                 rcdt_options.rcdt_max_edge,
                 rcdt_options.quality.raw(),
                 rcdt_options.quality.rows(),
                 rcdt_options.rcdt_max_n0,
                 rcdt_options.rcdt_max_n1);
         FMLOG(MC << endl);
       }
       /* Done constructing the triangulation. */
     }

     FMLOG("MC state = " << MC.getState() << std::endl);

     /* Calculate and collect output. */

     matrices.attach("segm.bnd.idx", std::make_unique<Matrix<int>>(2),
                     fmesh::IOMatrixtype::General);
     matrices.attach("segm.bnd.grp", std::make_unique<Matrix<int>>(1),
                     fmesh::IOMatrixtype::General);
     MC.segments(
       true,
       &matrices.DI("segm.bnd.idx"),
       &matrices.DI("segm.bnd.grp")
     );

     matrices.output("segm.bnd.idx").output("segm.bnd.grp");

     matrices.attach("segm.int.idx", std::make_unique<Matrix<int>>(2),
                     fmesh::IOMatrixtype::General);
     matrices.attach("segm.int.grp", std::make_unique<Matrix<int>>(1),
                     fmesh::IOMatrixtype::General);
     MC.segments(
       false,
       &matrices.DI("segm.int.idx"),
       &matrices.DI("segm.int.grp")
     );

     matrices.output("segm.int.idx").output("segm.int.grp");

     matrices.attach("tt", &M.TT());
     M.useVT(true);
     //  matrices.attach("vt", &M.VT());
     M.useTTi(true);
     matrices.attach("tti", &M.TTi());
     matrices.attach("vv", std::make_unique<SparseMatrix<int>>(M.VV()),
                     fmesh::IOMatrixtype::Symmetric);

     matrices.output("tt").output("tti").output("vt").output("vv");

     //  FMLOG("Manifold output." << std::endl);
     //  /* Output the manifold type. */
     //  matrices.attach("manifold", std::make_unique<Matrix<int>>(1),
     //                  fmesh::IOMatrixtype::General);
     //  Matrix<int> &manifold = matrices.DI("manifold");
     //  manifold(0, 0) = static_cast<int>(M.type());
     //  matrices.output("manifold");

     Rcpp::List out = Rcpp::wrap(matrices);

     // Add VT information
     out["vt"] = Rcpp::wrap(M.VT());
     // Rcpp::List vt(M.VT().size());
     // for (size_t i = 0; i < M.VT().size(); i++) {
     //   vt[i] = Rcpp::transpose(Rcpp::as<Rcpp::IntegerMatrix>(Rcpp::wrap(M.VT()[i])));
     // }
     // out["vt"] = vt;

     switch (M.type()) {
     case Mesh::Mtype::Manifold:
       out["manifold"] = "M2";
       break;
     case Mesh::Mtype::Plane:
       out["manifold"] = "R2";
       break;
     case Mesh::Mtype::Sphere:
       out["manifold"] = "S2";
       break;
     }

     return out;
   }




//' @title Barycentric coordinate computation
//'
//' @description
//' Locate points and compute triangular barycentric coordinates
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param loc numeric matrix; coordinates of points to locate in the mesh
//' @param options list of triangulation options
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' b <- fmesher_bary(m$s,
//'                   m$tv,
//'                   matrix(c(0.5, 0.5), 1, 2),
//'                   list())
//' @returns A list with vector `index` (triangle index) and matrix `where`
//' (3-column barycentric matrix)
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_bary(Rcpp::NumericMatrix mesh_loc,
                        Rcpp::IntegerMatrix mesh_tv,
                        Rcpp::NumericMatrix loc,
                        Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);
  Options rcdt_options(options, M.nV());

  FMLOG("barycentric coordinate output." << std::endl);
  if ((M.type() != Mesh::Mtype::Plane) &&
      (M.type() != Mesh::Mtype::Sphere)) {
    FMLOG_("Cannot currently calculate points2mesh mapping for non R2/S2 manifolds"
             << std::endl);
    return Rcpp::List();
  }

  matrices.attach("loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(loc))));
  Matrix<double>& points2mesh = matrices.DD("loc");

  size_t points_n = points2mesh.rows();
  Matrix<int> &points2mesh_t =
    matrices.attach(string("index"), std::make_unique<Matrix<int>>(points_n, 1));
  Matrix<double> &points2mesh_b = matrices.attach(
    string("where"), std::make_unique<Matrix<double>>(points_n, 3));
  matrices.matrixtype("index", fmesh::IOMatrixtype::General);
  matrices.matrixtype("where", fmesh::IOMatrixtype::General);
  matrices.output("index").output("where");

  FMLOG("map_points_to_mesh start" << std::endl);
  map_points_to_mesh(M, points2mesh, points2mesh_t, points2mesh_b);
  FMLOG("map_points_to_mesh done" << std::endl);

  return Rcpp::wrap(matrices);
}

//' @title Barycentric coordinate computation
//'
//' @description
//' Locate points and compute triangular barycentric coordinates
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param loc numeric matrix; coordinates of points to locate in the mesh
//' @param options list of triangulation options
//' @examples
//' m <- fmesher_mesh3d(list(cet_margin = 1),
//'                     matrix(rnorm(15), 5, 3),
//'                     matrix(c(0,1,2,3), 1, 4))
//' b <- fmesher_bary3d(m$loc,
//'                     m$tv,
//'                     matrix(c(0.5, 0.5, 0.5), 1, 3),
//'                     list())
//' @returns A list with vector `index` (tetra index) and matrix `where`
//' (4-column barycentric matrix)
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_bary3d(Rcpp::NumericMatrix mesh_loc,
                          Rcpp::IntegerMatrix mesh_tv,
                          Rcpp::NumericMatrix loc,
                          Rcpp::List options) {
  MatrixC matrices;
  Mesh3 M = Rcpp_import_mesh3d(mesh_loc, mesh_tv, matrices, options);

  FMLOG("barycentric coordinate output." << std::endl);

  matrices.attach("loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(loc))));
  Matrix<double>& points2mesh = matrices.DD("loc");

  size_t points_n = points2mesh.rows();
  Matrix<int> &points2mesh_t =
    matrices.attach(string("index"), std::make_unique<Matrix<int>>(points_n, 1));
  Matrix<double> &points2mesh_b = matrices.attach(
    string("where"), std::make_unique<Matrix<double>>(points_n, 4));
  matrices.matrixtype("index", fmesh::IOMatrixtype::General);
  matrices.matrixtype("where", fmesh::IOMatrixtype::General);
  matrices.output("index").output("where");

  FMLOG("map_points_to_mesh start" << std::endl);
  map_points_to_mesh3d(M, points2mesh, points2mesh_t, points2mesh_b);
  FMLOG("map_points_to_mesh done" << std::endl);

  return Rcpp::wrap(matrices);
}




//' @title Rotationally invariant spherical B-splines
//'
//' @description
//' Compute rotationally invariant spherical B-splines on the unit sphere
//'
//' @param loc numeric vector/matrix; coordinates of points to locate in the mesh,
//' only the z-coordinates are used (`sin(latitude)`)
//' @param n The number of basis functions
//' @param degree The polynomial basis degree
//' @param uniform logical; If `TRUE`, the knots are spaced uniformly by latitude,
//' if `FALSE`, the knots are spaced uniformly by `sin(latitude)`
//' @rdname fmesher_spherical_bsplines
//' @examples
//' m <- fm_rcdt_2d(globe = 1)
//' fmesher_spherical_bsplines(m$loc, n = 3, degree = 2, uniform = FALSE)
//' fmesher_spherical_bsplines1(m$loc[, 3], n = 3, degree = 2, uniform = FALSE)
//' @export
//' @keywords internal
//' @returns A matrix of evaluated b-spline basis functions
// [[Rcpp::export]]
SEXP fmesher_spherical_bsplines1(Rcpp::NumericVector loc,
                                 int n,
                                 int degree,
                                 Rcpp::LogicalVector uniform) {
  if (n < 0) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (degree < 1) {
    Rcpp::stop("'degree' must be at least 0.");
  }
  if (n <= degree) {
    Rcpp::stop("'n' must be larger than 'degree'");
  }

  MatrixC matrices;
  matrices.attach("loc", std::make_unique<Matrix<double>>(loc));

  FMLOG("bspline output." << std::endl);

  bool bool_uniform = Rcpp::is_true(Rcpp::all(uniform));
  if (bool_uniform) {
    FMLOG("uniform = TRUE" << std::endl);
  } else {
    FMLOG("uniform = FALSE" << std::endl);
  }
  matrices.attach(
    string("bspline"),
    spherical_bsplines1(matrices.DD("loc"), n, degree, bool_uniform));
  matrices.matrixtype("bspline", fmesh::IOMatrixtype::General);
  matrices.output("bspline");

  return Rcpp::wrap(matrices.DD("bspline"));
}

//' @rdname fmesher_spherical_bsplines
//' @export
// [[Rcpp::export]]
SEXP fmesher_spherical_bsplines(Rcpp::NumericMatrix loc,
                                int n,
                                int degree,
                                Rcpp::LogicalVector uniform) {
  if (n < 0) {
    Rcpp::stop("'n' must be at least 1.");
  }
  if (degree < 1) {
    Rcpp::stop("'degree' must be at least 0.");
  }
  if (n <= degree) {
    Rcpp::stop("'n' must be larger than 'degree'");
  }
  if (loc.cols() < 3) {
    Rcpp::stop("'ncol(loc)' must be at least 3.");
  }

  MatrixC matrices;
  matrices.attach("loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(loc))));

  FMLOG("bspline output." << std::endl);

  bool bool_uniform = Rcpp::is_true(Rcpp::all(uniform));
  if (bool_uniform) {
    FMLOG("uniform = TRUE" << std::endl);
  } else {
    FMLOG("uniform = FALSE" << std::endl);
  }
  matrices.attach(
    string("bspline"),
    spherical_bsplines(matrices.DD("loc"), n, degree, bool_uniform));
  matrices.matrixtype("bspline", fmesh::IOMatrixtype::General);
  matrices.output("bspline");

  return Rcpp::wrap(matrices.DD("bspline"));
  //  return Rcpp::wrap(matrices);
}




//' @title Finite element matrix computation
//'
//' @description
//' Construct finite element structure matrices
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param fem_order_max integer; the highest operator order to compute
//' @param aniso If non-NULL, a `list(gamma, v)`. Calculates anisotropic structure
//' matrices (in addition to the regular) for \eqn{\gamma}{gamma} and \eqn{v}{v} for
//' an anisotropic operator \eqn{\nabla\cdot H \nabla}{div H grad}, where
//' \eqn{H=\gamma I + v v^\top}{H = gamma I + v v'}.
//' Currently (2023-08-05) the fields need to be given per vertex.
//' @param options list of triangulation options (`sphere_tolerance`)
//' @examples
//' m <- fmesher_rcdt(list(cet_margin = 1), matrix(0, 1, 2))
//' b <- fmesher_fem(m$s, m$tv, fem_order_max = 2, aniso = NULL, options = list())
//' @returns A list of matrices
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_fem(Rcpp::NumericMatrix mesh_loc,
                       Rcpp::IntegerMatrix mesh_tv,
                       int fem_order_max,
                       Rcpp::Nullable<Rcpp::List> aniso,
                       Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);

  FMLOG("Compute finite element matrices." << std::endl);

  if (fem_order_max >= 0) {
    FMLOG("fem output." << std::endl)
    SparseMatrix<double> &C0 = matrices.SD("c0").clear();
    SparseMatrix<double> &C1 = matrices.SD("c1").clear();
    SparseMatrix<double> &B1 = matrices.SD("b1").clear();
    SparseMatrix<double> &G = matrices.SD("g1").clear();
    SparseMatrix<double> &K = matrices.SD("k1").clear();
    /* K1=G1-B1, K2=K1*inv(C0)*K1, ... */
    Matrix<double> &Tareas = matrices.DD("ta").clear();

    M.calcQblocks(C0, C1, G, B1, Tareas);

    matrices.attach(string("va"), std::make_unique<Matrix<double>>(diag(C0)));

    K = G - B1;

    matrices.matrixtype("c0", fmesh::IOMatrixtype::Diagonal);
    matrices.matrixtype("c1", fmesh::IOMatrixtype::Symmetric);
    matrices.matrixtype("b1", fmesh::IOMatrixtype::General);
    matrices.matrixtype("g1", fmesh::IOMatrixtype::Symmetric);
    matrices.matrixtype("k1", fmesh::IOMatrixtype::General);
    matrices.output("c0");
    matrices.output("c1");
    matrices.output("b1");
    matrices.output("g1");
    matrices.output("k1");
    matrices.output("va");
    matrices.output("ta");

    SparseMatrix<double> C0inv = inverse(C0, true);
    // Protect temporary local variables
    {
      SparseMatrix<double> tmp = G * C0inv;
      SparseMatrix<double> *a;
      SparseMatrix<double> *b = &G;
      for (size_t i = 1; int(i) < fem_order_max; i++) {
        std::stringstream ss;
        ss << i + 1;
        std::string Gname = "g" + ss.str();
        a = b;
        b = &(matrices.SD(Gname).clear());
        *b = tmp * (*a);
        matrices.matrixtype(Gname, fmesh::IOMatrixtype::Symmetric);
        matrices.output(Gname);
      }
      tmp = C0inv * K;
      b = &K;
      for (size_t i = 1; int(i) < fem_order_max; i++) {
        std::stringstream ss;
        ss << i + 1;
        std::string Kname = "k" + ss.str();
        a = b;
        b = &(matrices.SD(Kname).clear());
        *b = (*a) * tmp;
        matrices.matrixtype(Kname, fmesh::IOMatrixtype::General);
        matrices.output(Kname);
      }
    }

    if (!aniso.isNull()) {
      FMLOG("Compute anisotropic finite element matrices." << std::endl);
      if (Rcpp::as<Rcpp::List>(aniso).size() < 2) {
        Rcpp::stop("'aniso' list must have at least two elements.");
      }
      matrices.attach("gamma_field",
                      std::make_unique<Matrix<double>>(
                          Rcpp::as<Rcpp::NumericVector>(
                            Rcpp::as<Rcpp::List>(aniso)[0]
                          )
                      ));
      FMLOG("'gamma_field' imported." << std::endl);
      matrices.attach("vector_field",
                      std::make_unique<Matrix<double>>(
                          Rcpp::as<Rcpp::NumericMatrix>(
                            Rcpp::as<Rcpp::List>(aniso)[1]
                          )
                      ));
      FMLOG("'vector_field' imported." << std::endl);
      if (matrices.DD("gamma_field").rows() < M.nV()) {
        Rcpp::stop("'aniso[[1]]' length should match the number of vertices.");
      }
      if (matrices.DD("vector_field").rows() < M.nV()) {
        Rcpp::stop("'aniso[[2]]' rows should match the number of vertices.");
      }

      SparseMatrix<double> &Gani = matrices.SD("g1aniso").clear();
      M.calcQblocksAni(Gani,
                       matrices.DD("gamma_field"),
                       matrices.DD("vector_field"));
      matrices.output("g1aniso");

      // Protect temporary local variables
      {
        SparseMatrix<double> tmp = Gani * C0inv;
        SparseMatrix<double> *a;
        SparseMatrix<double> *b = &Gani;
        for (size_t i = 1; int(i) < fem_order_max; i++) {
          std::stringstream ss;
          ss << i + 1;
          std::string Gname = "g" + ss.str() + "aniso";
          a = b;
          b = &(matrices.SD(Gname).clear());
          *b = tmp * (*a);
          matrices.matrixtype(Gname, fmesh::IOMatrixtype::Symmetric);
          matrices.output(Gname);
        }
      }
    }

  }

  return Rcpp::wrap(matrices);
}



//' @title Split lines at triangle edges
//'
//' @description
//' Split a sequence of line segments at triangle edges
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param loc numeric coordinate matrix
//' @param idx 2-column integer matrix
//' @param options list of triangulation options (`sphere_tolerance`)
//' @export
//' @returns A list of line splitting information objects
//' @seealso [fm_split_lines()]
//' @examples
//' mesh <- fm_mesh_2d(
//'   boundary = fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0, 1)), is.bnd = TRUE)
//' )
//' splitter <- fm_segm(rbind(c(0.8, 0.2), c(0.2, 0.8)))
//' segm_split <- fm_split_lines(mesh, splitter)
// [[Rcpp::export]]
Rcpp::List fmesher_split_lines(
    Rcpp::NumericMatrix mesh_loc,
    Rcpp::IntegerMatrix mesh_tv,
    Rcpp::NumericMatrix loc,
    Rcpp::IntegerMatrix idx,
    Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, options);

  FMLOG("Compute line splitting." << std::endl);

  matrices.attach("loc",
                  std::make_unique<Matrix<double>>(Matrix3double(Matrix<double>(loc))));
  matrices.attach("idx",
                  std::make_unique<Matrix<int>>(idx));

  /* Make sure we have a Nx3 matrix: */
  auto splitloc1 = std::make_unique<Matrix<double>>(3);
  auto splitidx1 = std::make_unique<Matrix<int>>(2);
  auto splittriangle1 = std::make_unique<Matrix<int>>(1);
  auto splitbary1 = std::make_unique<Matrix<double>>(3);
  auto splitbary2 = std::make_unique<Matrix<double>>(3);
  auto splitorigin1 = std::make_unique<Matrix<int>>(1);

  split_line_segments_on_triangles(
    M, matrices.DD("loc"), matrices.DI("idx"), *splitloc1,
    *splitidx1, *splittriangle1, *splitbary1, *splitbary2, *splitorigin1);

  /* Now it's ok to overwrite potential input split* matrices. */
  matrices.attach("split.loc", std::move(splitloc1));
  matrices.attach("split.idx", std::move(splitidx1));
  matrices.attach("split.t", std::move(splittriangle1));
  matrices.attach("split.b1", std::move(splitbary1));
  matrices.attach("split.b2", std::move(splitbary2));
  matrices.attach("split.origin", std::move(splitorigin1));
  matrices.output("split.loc").output("split.idx");
  matrices.output("split.b1").output("split.b2");
  matrices.output("split.t").output("split.origin");

  return Rcpp::wrap(matrices);
}



typedef struct {
  int start;
  int finish;
  int sequence_index; // 0, ..., subdivisions - 1
} edge_point_t;
bool operator==(const edge_point_t& lhs, const edge_point_t& rhs) {
  return lhs.start == rhs.start && lhs.finish == rhs.finish &&
    lhs.sequence_index == rhs.sequence_index;
}
template<>
struct std::hash<edge_point_t>
{
  std::size_t operator()(const edge_point_t& ep) const noexcept
  {
    std::size_t h1 = std::hash<int>{}(ep.start);
    std::size_t h2 = std::hash<int>{}(ep.finish);
    std::size_t h3 = std::hash<int>{}(ep.sequence_index);
    return ((h1 ^ (h2 << 1)) >>  1) ^ (h3 << 1);
  }
};




//' @title Subdivide triangles
//'
//' @description
//' Subdivide a mesh with congruent and anti-congruent subtriangles
//'
//' @param mesh_loc numeric matrix; mesh vertex coordinates
//' @param mesh_tv 3-column integer matrix with 0-based vertex indices for each triangle
//' @param mesh_boundary 2-column integer matrix with 0-based vertex indices for
//' boundary constraints, currently ignored
//' @param mesh_interior 2-column integer matrix with 0-based vertex indices for
//' interior constraints, currently ignored
//' @param subdivisions integer; number of new points along each edge.
//' @param options list of triangulation options (`sphere_tolerance`)
//' @returns A list of new `loc` and `tv` information
//' @keywords internal
//' @seealso [fm_subdivide()]
//' @examples
//' mesh <- fm_mesh_2d(
//'   boundary = fm_segm(rbind(c(0,0), c(1,0), c(1,1), c(0, 1)), is.bnd = TRUE)
//' )
//' new_mesh <- fm_subdivide(mesh, n = 3)
//' plot(new_mesh, edge.color = 2)
//' plot(mesh, add = TRUE, edge.color = 1)
// [[Rcpp::export]]
Rcpp::List fmesher_subdivide(
    Rcpp::NumericMatrix mesh_loc,
    Rcpp::IntegerMatrix mesh_tv,
    Rcpp::IntegerMatrix mesh_boundary,
    Rcpp::IntegerMatrix mesh_interior,
    int subdivisions,
    Rcpp::List options) {
  MatrixC matrices;
  Mesh M = Rcpp_import_mesh(mesh_loc, mesh_tv, matrices, Rcpp::List());

  FMLOG("Subdivide triangles." << std::endl);

  matrices.attach("loc",
                  std::make_unique<Matrix<double>>(Matrix3double(mesh_loc)));
  matrices.attach("tv",
                  std::make_unique<Matrix<int>>(Matrix3int()));
  Matrix<double> &loc = matrices.DD("loc");
  Matrix<int> &tv = matrices.DI("tv");

  loc.capacity(M.nT() * (subdivisions + 2) * (subdivisions + 3) / 2);
  tv.rows(M.nT() * (subdivisions + 1) * (subdivisions + 1));

  FMLOG("Construct new points." << std::endl);
  // Points along edge
  std::unordered_map<edge_point_t, int> edge_to_point;

  FMLOG("Construct triangles." << std::endl);
  std::vector<std::vector<int>> triangle_to_point(subdivisions + 2);
  for (int i = 0; i < subdivisions + 2; i++) {
    triangle_to_point[i].resize(subdivisions + 2 - i, 0);
  }

  size_t t_offset = 0;
  for (size_t t = 0; t < M.nT(); t++) {
    FMLOG("Corner points" << std::endl);

    triangle_to_point[0][0] = M.TV()(t, 0);
    triangle_to_point[subdivisions + 1][0] = M.TV()(t, 1);
    triangle_to_point[0][subdivisions + 1] = M.TV()(t, 2);

    FMLOG("Edge points" << std::endl);

    for (int i = 0; i < 3; i++) {
      for (int k = 0; k < subdivisions; ++k) {
        int new_vtx = loc.rows();
        const std::unordered_map<edge_point_t, int>::const_iterator new_vtx_iter =
          edge_to_point.find(
          edge_point_t{M.TV()(t, (i+1) % 3), M.TV()(t, i), subdivisions - 1 - k});
        if (new_vtx_iter != edge_to_point.end()) {
          new_vtx = new_vtx_iter->second;
        } else {
          FMLOG("vtx 0: " << loc(M.TV()(t, i),0) << ", "
                           << loc(M.TV()(t, i),1) << ", "
                           << loc(M.TV()(t, i),2) << std::endl);
          FMLOG("vtx 1: " << loc(M.TV()(t, (i+1) % 3),0) << ", "
                           << loc(M.TV()(t, (i+1) % 3),1) << ", "
                           << loc(M.TV()(t, (i+1) % 3),2) << std::endl);
          FMLOG("subs - k = " << subdivisions - k << std::endl);
          FMLOG("k + 1    = " << k + 1 << std::endl);
          loc.rows(new_vtx + 1);
          for (int dim = 0; dim < 3; ++dim) {
            loc(new_vtx, dim) = (loc(M.TV()(t, i), dim) * (subdivisions - k) +
              loc(M.TV()(t, (i+1) % 3), dim) * (k + 1)) / (subdivisions + 1);
          }
          FMLOG("interp: " << loc(new_vtx,0) << ", "
                            << loc(new_vtx,1) << ", "
                            << loc(new_vtx,2) << std::endl);
          // TODO: project onto sphere for S2 meshes
          // Now done in the R interface fm_subdivide
          // FMLOG("TODO: fmesher_subdivide: handle spherical meshes." << std::endl);
        }
        edge_to_point.emplace(
              edge_point_t{M.TV()(t, i), M.TV()(t, (i+1) % 3), k},
              new_vtx);

        if (i == 0) {
          triangle_to_point[k + 1][0] =  new_vtx;
        } else if (i == 1) {
          triangle_to_point[subdivisions - k][k + 1] =  new_vtx;
        } else { // i == 2
          triangle_to_point[0][subdivisions - k] =  new_vtx;
        }
      }
    }

    FMLOG("Vertices = " << loc.rows() << std::endl);

    FMLOG("Interior points" << std::endl);

    // interior points:
    // j = 1, ..., subs - 1
    // i = 1, ..., subs - j
    // total : subs * (subs - 1) / 2
    int new_vtx = loc.rows();
    loc.rows(loc.rows() + subdivisions * (subdivisions - 1) / 2);
    for (int j = 1; j < subdivisions; ++j) {
      for (int i = 1; i < subdivisions - j + 1; ++i) {
        int k = subdivisions + 1 - i - j;
        for (int dim = 0; dim < 3; ++dim) {
          loc(new_vtx, dim) = (
            loc(M.TV(t)[0], dim) * k +
              loc(M.TV(t)[1], dim) * i +
              loc(M.TV(t)[2], dim) * j) / (subdivisions + 1);
        }
        triangle_to_point[i][j] = new_vtx;
        // TODO: project onto sphere for S2 meshes
        // Now done in the R interface fm_subdivide
        // FMLOG("TODO: fmesher_subdivide: handle spherical meshes." << std::endl);
        ++new_vtx;
      }
    }

    FMLOG("Triangle " << t << " = (" << M.TV()(t,0) <<
      "," << M.TV()(t,1) << "," << M.TV()(t,2) << ")" << std::endl);

    FMLOG("Triangles" << std::endl);

    if (subdivisions <= 0) {
      tv(t_offset, 0) = M.TV(t)[0];
      tv(t_offset, 1) = M.TV(t)[1];
      tv(t_offset, 2) = M.TV(t)[2];
      ++t_offset;
    } else {
      // congruent triangles; (subs + 1) * (subs + 2) / 2
      // converse triangles; subs * (subs + 1) / 2
      // total : (subs + 1)^2
      for (int j = 0; j < subdivisions + 1; ++j) {
        // congruent triangles; subs + 1 - j
        // converse triangles; subs - j
        for (int i = 0; i < subdivisions - j + 1; ++i) {
          // congruent triangle:
          tv(t_offset, 0) = triangle_to_point[i][j];
          tv(t_offset, 1) = triangle_to_point[i + 1][j];
          tv(t_offset, 2) = triangle_to_point[i][j + 1];
          ++t_offset;
        }
        for (int i = 0; i < subdivisions - j; ++i) {
          // converse triangle:
          tv(t_offset, 0) = triangle_to_point[i + 1][j + 1];
          tv(t_offset, 1) = triangle_to_point[i][j + 1];
          tv(t_offset, 2) = triangle_to_point[i + 1][j];
          ++t_offset;
        }
      }
    }
  }

  FMLOG("Done" << std::endl);

  FMLOG("Vertices = " << loc.rows() << std::endl);
  FMLOG("Triangles = " << tv.rows() << std::endl);

  matrices.output("loc").output("tv");

  return Rcpp::wrap(matrices);
}


// //' @title Test the matrix I/O system
// //'
// //' @param args_input Input argument list
// //' @examples
// //' \dontrun{
// //' A <- Matrix::sparseMatrix(i=1:4,j=4:1,x=2:5,dims=c(4,4))
// //' inp <- list(
// //'   A = fm_as_dgTMatrix(A),
// //'   Bd = matrix((11:22)+0.5,4,3),
// //'   Bi = matrix(121L:132L,4,3),
// //'   B1d=as.matrix((31:34)+0.5),
// //'   B1i=as.matrix(41L:44L),
// //'   Ad = fm_as_fmesher_sparse(A)
// //' )
// //' inp[["BdM"]] <- as(inp[["Bd"]], "unpackedMatrix")
// //' out <- C_matrixio_test2(args_input = inp)
// //' str(out)
// //' }
// //' @keywords internal
// // [[Rcpp::export]]
// Rcpp::List C_matrixio_test2(Rcpp::List args_input) {
//   MatrixC matrices(args_input);
//   Rcpp::List ret = Rcpp::wrap(matrices.output("-"));
//   return (ret);
// }



// //' @title Test the matrix I/O system
// //'
// //' @param args_input Input argument list
// //' @examples
// //' \dontrun{
// //' A <- Matrix::sparseMatrix(i=1:4,j=4:1,x=2:5,dims=c(4,4))
// //' out <- C_matrixio_test(args_input=list(
// //'   A = fm_as_dgTMatrix(A),
// //'   Bd = matrix((11:22)+0.5,4,3),
// //'   Bi = matrix(121L:132L,4,3),
// //'   B1d=as.matrix((31:34)+0.5),
// //'   B1i=as.matrix(41L:44L),
// //'   Ad = fm_as_fmesher_sparse(A)
// //' ))
// //' Aout <- fm_as_dgTMatrix(out[["Ad"]])
// //' A
// //' Aout
// //' }
// //' @keywords internal
// // [[Rcpp::export]]
//      Rcpp::List C_matrixio_test(Rcpp::List args_input) {
//        MatrixC matrices;
//
//        //  matrices.attach("loc", std::make_unique<Matrix<double>>(Rcpp::as<EigenMM<double>>(args_input["loc"])));
//        //  matrices.attach("tv", std::make_unique<Matrix<int>>(Rcpp::as<EigenMM<int>>(args_input["tv"])));
//
//        bool is_list = Rcpp::is<Rcpp::List>(args_input);
//        bool is_numeric_matrix = Rcpp::is<Rcpp::NumericMatrix>(args_input["A"]);
//        bool is_numeric_vector = Rcpp::is<Rcpp::NumericVector>(args_input["A"]);
//        bool is_integer_matrix = Rcpp::is<Rcpp::IntegerMatrix>(args_input["A"]);
//        bool is_integer_vector = Rcpp::is<Rcpp::IntegerVector>(args_input["A"]);
//
//        Rcpp::NumericMatrix Bd = Rcpp::as<Rcpp::NumericMatrix>(args_input["Bd"]);
//        Rcpp::IntegerMatrix Bi = Rcpp::as<Rcpp::IntegerMatrix>(args_input["Bi"]);
//        Rcpp::NumericVector B1d = Rcpp::as<Rcpp::NumericVector>(args_input["B1d"]);
//        Rcpp::IntegerVector B1i = Rcpp::as<Rcpp::IntegerVector>(args_input["B1i"]);
//
//
//        fmesh::Matrix<double> Bdd = Bd;
//        fmesh::Matrix<int> Bdi(Rcpp::as<Rcpp::IntegerMatrix>(Bd));
//        fmesh::Matrix<double> Bid(Rcpp::as<Rcpp::NumericMatrix>(Bi));
//        fmesh::Matrix<int> Bii(Bi);
//
//        FMLOG_("Bdd: " << Bdd << std::endl);
//        FMLOG_("Bdi: " << Bdi << std::endl);
//        FMLOG_("Bid: " << Bid << std::endl);
//        FMLOG_("Bii: " << Bii << std::endl);
//
//        fmesh::Matrix1<double> Bdd1 = B1d;
//        fmesh::Matrix1<double> Bdd_1 = Rcpp::NumericVector(Bd(Rcpp::_, 1));
//        fmesh::Matrix1<int> Bdi1(Rcpp::as<Rcpp::IntegerVector>(B1d));
//        fmesh::Matrix1<double> Bid1(Rcpp::as<Rcpp::NumericVector>(B1i));
//        fmesh::Matrix1<int> Bii1(B1i);
//
//        FMLOG_("Bdd1: " << Bdd1 << std::endl);
//        FMLOG_("Bdd_1: " << Bdd_1 << std::endl);
//        FMLOG_("Bdi1: " << Bdi1 << std::endl);
//        FMLOG_("Bid1: " << Bid1 << std::endl);
//        FMLOG_("Bii1: " << Bii1 << std::endl);
//
//        fmesh::Matrix3<double> Bdd3 = Bd;
//        fmesh::Matrix3<int> Bdi3(Rcpp::as<Rcpp::IntegerMatrix>(Bd));
//        fmesh::Matrix3<double> Bid3(Rcpp::as<Rcpp::NumericMatrix>(Bi));
//        fmesh::Matrix3<int> Bii3(Bi);
//
//        FMLOG_("Bdd3: " << Bdd3 << std::endl);
//        FMLOG_("Bdi3: " << Bdi3 << std::endl);
//        FMLOG_("Bid3: " << Bid3 << std::endl);
//        FMLOG_("Bii3: " << Bii3 << std::endl);
//
//        const Rcpp::List Ad(Rcpp::as<Rcpp::List>(args_input["Ad"]));
//
//        fmesh::SparseMatrix<double> Ad_fm(Ad);
//
//        //  const EigenMSM<int> Ai(Rcpp::as<EigenMSM<int>>(args_input["Ai"]));
//
//        //  bool is_msm = Rcpp::is<Eigen::SparseMatrix<double>>(args_input["a"]);
//
//        matrices.attach("Ad_fm", &Ad_fm);
//        matrices.output("Ad_fm");
//
//        MatrixC mat2(args_input);
//
//        FMLOG_(mat2.DI("tv"))
//
//          Rcpp::List ret;
//        ret["is_list"] = is_list;
//        ret["is_numeric_matrix"] = is_numeric_matrix;
//        ret["is_numeric_vector"] = is_numeric_vector;
//        ret["is_integer_matrix"] = is_integer_matrix;
//        ret["is_integer_vector"] = is_integer_vector;
//        //  ret["A"] = A;
//        ret["Ad"] = Ad;
// #ifdef FMESHER_WITH_EIGEN
//        ret["Ad_fm"] = Ad_fm.EigenSparseMatrix();
// #endif
//        ret["Ad_fm_auto"] = Ad_fm;
//        ret["Ad_fm_ijx"] = Ad_fm.fmesher_sparse();
//        ret["Bid3"] = Bid3;
//        ret["Bdi3"] = Bdi3;
//        ret["matrices"] = matrices;
//        ret["mat2"] = mat2.output("-");
//        return (ret);
//      }



//' @title 3D tetrahedralisation storage
//'
//' @description
//' (...)
//'
//' @param options list of triangulation options
//' @param loc numeric matrix; initial points to include
//' @param tv 4-column integer matrix with 0-based vertex indices for each triangle
//' @examples
//' m <- fmesher_mesh3d(list(),
//'                     matrix(c(1,0,0,0,1,0,0,0,1,0,0,0), 4, 3, byrow=TRUE),
//'                     matrix(c(0,1,2,3), 1, 4, byrow=TRUE))
//' @returns A list of information objects for a generated tetrahedralisation
//' @export
// [[Rcpp::export]]
Rcpp::List fmesher_mesh3d(Rcpp::List options,
                           Rcpp::NumericMatrix loc,
                           Rcpp::IntegerMatrix tv) {
   MatrixC matrices;
   Mesh3 M = Rcpp_import_mesh3d(loc, tv, matrices, options);

   FMLOG("Attach 'loc'." << std::endl);
   matrices.attach(string("loc"), &M.S());
   FMLOG("Attach 'tv'." << std::endl);
   matrices.attach("tv", &M.TV());
   FMLOG("Set output of 'loc' and 'tv'." << std::endl);
   matrices.output("loc").output("tv");

   FMLOG("M = " << M << std::endl);
   FMLOG("M.TT = " << M.TT() << std::endl);
   FMLOG("M.TTi = " << M.TTi() << std::endl);
   FMLOG("M.VV = " << M.VV() << std::endl);

   matrices.attach("tt", &M.TT());
   M.useVT(true);
   //  matrices.attach("vt", &M.VT());
   M.useTTi(true);
   matrices.attach("tti", &M.TTi());
   matrices.attach("vv", std::make_unique<SparseMatrix<int>>(M.VV()),
                   fmesh::IOMatrixtype::Symmetric);

   matrices.output("tt").output("tti").output("vt").output("vv");

   Rcpp::List out = Rcpp::wrap(matrices);

   switch (M.type()) {
   case Mesh3::Mtype::Manifold:
     out["manifold"] = "M3";
     break;
   case Mesh3::Mtype::Plane:
     out["manifold"] = "R3";
     break;
   }

   return out;
 }





#endif
