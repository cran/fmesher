/*
 *  Copyright Finn Lindgren (2010-2024)
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public License,
 *  v. 2.0. If a copy of the MPL was not distributed with this file, You can
 *  obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef _FMESH_MESH_
#define _FMESH_MESH_ 1

#include <cstddef>
//#include <cstring>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#ifdef FMESHER_WITH_X
#include "x11utils.h"
#endif
#include "fmesher_debuglog.h"
#include "vector.h"

#define MESH_EPSILON 1e-15

namespace fmesh {

class Xtmpl;
class Mesh;
class Dart;
class MOAint;
class MOAVTMap;
class MOAVTMapV;
class MOAint3;
class MOAdouble3;
class MeshC;

typedef std::pair<int, int> IntPair;
typedef std::list<int> vertexListT;
typedef std::set<int> triangleSetT;
typedef std::list<Dart> DartList;
typedef std::pair<Dart, Dart> DartPair;
typedef std::map<int,int> IntMap;
typedef std::vector<IntMap> VTMapT;

class Mesh {
  friend class Dart;
  friend std::ostream &operator<<(std::ostream &output, const Mesh &M);

public:
  enum class Mtype : int { Manifold = 0, Plane, Sphere };
  friend std::ostream &operator<<(std::ostream &output, const Mesh::Mtype &type);

private:
  Mtype type_;
  double sphere_radius_;
  bool use_VT_;
  bool use_TTi_;
  Matrix3int TV_;  /* TV[t]  : {v1,v2,v3} */
  Matrix3int TT_;  /* TT[t]  : {t1,t2,t3} */
  VTMapT VT_mapping_;  /* VT[v] : map from t to the vi (0,1,2) value for v */
  /* TV[ VT[v][i].t, VT[v][i].vi ] == v */
  Matrix3int TTi_; /* TTi[t] : {vi1,vi2,vi3},
                     t == TT[ TT[t][i] ][ TTi[t][i] ] */
  Matrix3double S_;
#ifdef FMESHER_WITH_X
  Xtmpl *X11_;
  int X11_v_big_limit_;
  int verbose_;
#endif

private:
  Mesh &rebuildTT();

  /*! Add triangle to VT[v] */
  Mesh &add_VT(const int v, const int t);
  /*! Add triangle to VT[v] with precomputed vi information */
  Mesh &add_VT(const int v, const int t, const int vi);
  /*! Remove triangle to VT[v] */
  Mesh &remove_VT(const int v, const int t);
  /*! Clear VT[v] info */
  Mesh &clear_VT(const int v);
  /*! Set VT[v]=empty for v>=v_start */
  Mesh &reset_VT(const int v_start = 0);

  Mesh &add_VT_triangle(const int t);
  Mesh &remove_VT_triangle(const int t);
  Mesh &add_VT_triangles(const int t_start = 0);
  Mesh &remove_VT_triangles(const int t_start = 0);
  void check_VT_mapping_consistency() const;

  Mesh &rebuild_VT();
  Mesh &rebuildTTi();

#ifdef FMESHER_WITH_X
public:
  void drawX11point(int v, bool fg);
  void drawX11triangle(int t, bool fg);
  void redrawX11(std::string str);
#endif

public:
  Mesh(void)
      : type_(Mtype::Manifold), use_VT_(false), use_TTi_(true), TV_(), TT_(),
        VT_mapping_(), TTi_(), S_()
#ifdef FMESHER_WITH_X
        ,
        X11_(NULL), X11_v_big_limit_(0), verbose_(0)
#endif
                                             {};
  Mesh(Mtype manifold_type, size_t Vcapacity, bool use_VT = true,
       bool use_TTi = false);
  Mesh(const Mesh &M)
      : type_(Mtype::Manifold), sphere_radius_(1.0), use_VT_(true),
        use_TTi_(false), TV_(), TT_(), VT_mapping_(), TTi_(), S_()
#ifdef FMESHER_WITH_X
        ,
        X11_(NULL), X11_v_big_limit_(0), verbose_(0)
#endif
  {
    *this = M;
  };
  Mesh &operator=(const Mesh &M);
  ~Mesh();
  Mesh &clear();
  Mesh &empty();

  /*!
    \brief Check the storage capacity, and increase if necessary
  */
  Mesh &check_capacity(size_t nVc, size_t nTc);
  size_t Vcap() const { return S_.capacity(); }

  bool useVT() const { return use_VT_; };
  Mesh &useVT(bool use_VT);
  bool useTTi() const { return use_TTi_; };
  Mesh &useTTi(bool use_TTi);

#ifdef FMESHER_WITH_X
  bool useX11() const {
    return (X11_ != NULL);
    return false;
  };
  void setX11VBigLimit(int lim) { X11_v_big_limit_ = lim; };
  void setX11delay(double set_delay);
  Mesh &useX11(bool use_X11, bool draw_text, int sx = 500, int sy = 500,
               double minx = -0.05, double maxx = 1.05, double miny = -0.05,
               double maxy = 1.05, std::string name = "fmesher::Mesh");
#endif

  Mtype type() const { return type_; };
  void type(Mtype set_type) { type_ = set_type; };
  double sphere_radius() const { return sphere_radius_; };
  void sphere_radius(double default_sphere_radius) {
    if (S_.rows() == 0) {
      sphere_radius_ = default_sphere_radius;
    } else {
      sphere_radius_ = S_(0).length();
    }
  };

  Mtype determine_type(double sphere_tolerance) {
    bool issphere = false;
    bool isflat = false;
    if ((*this).S().rows() == 0) {
      FMLOG("Empty manifold detected, treating as flat." << std::endl);
      return Mtype::Plane;
    }
    isflat = (std::fabs((*this).S(0)[2]) < 1.0e-10);
    double radius = (*this).S(0).length();
    issphere = (radius > sphere_tolerance);
    for (size_t i = 1; i < (*this).nV(); i++) {
      isflat = (isflat && (std::fabs((*this).S(i)[2]) < 1.0e-10));
      if (issphere) {
        issphere =
          (std::fabs((*this).S(i).length() / radius - 1.0) < sphere_tolerance);
      }
    }
    if (isflat) {
      FMLOG("Plane detected" << std::endl);
      return Mtype::Plane;
    } else if (issphere) {
      FMLOG("Sphere detected" << std::endl);
      return Mtype::Sphere;
    } else {
      FMLOG("General manifold detected" << std::endl);
      return Mtype::Manifold;
    }
  }


  Mtype auto_type(double sphere_tolerance) {
    Mtype mtype = determine_type(sphere_tolerance);
    (*this).type(mtype);
    if (mtype == fmesh::Mesh::Mtype::Sphere) {
      // The radius argument gives a default radius,
      // used only if there are no points.
      sphere_radius(1.0);
    }
    return mtype;
  }




  size_t nV() const { return S_.rows(); };
  size_t nT() const { return TV_.rows(); };
  const Matrix3int &TV() const { return TV_; };
  const Matrix3int &TT() const { return TT_; };
  const VTMapT &VT() const { return VT_mapping_; };
  const Matrix3int &TTi() const { return TTi_; };
  const Matrix3double &S() const { return S_; };

  Matrix3int &TV() { return TV_; };
  Matrix3int &TT() { return TT_; };
  VTMapT &VT() { return VT_mapping_; };
  Matrix3int &TTi() { return TTi_; };
  Matrix3double &S() { return S_; };

  SparseMatrix<int> VV() const;
  const Int3 &TV(int t) const { return TV_[t]; };
  const Int3 &TT(int t) const { return TT_[t]; };
  const IntMap &VT(int v) const { return VT_mapping_[v]; };
  const Int3 &TTi(int t) const { return TTi_[t]; };
  const Point &S(int v) const { return S_[v]; };

#ifdef FMESHER_WITH_X
  Xtmpl *X11() { return X11_; };
#endif

  MOAint3 TVO() const;
  MOAint3 TTO() const;
  MOAVTMap VTO() const;
  MOAVTMapV VTO(const int v) const;
  MOAint3 TTiO() const;
  MOAdouble3 SO() const;

  Mesh &S_set(const Matrix3double &S);
  Mesh &TV_set(const Matrix3int &TV);
  Mesh &S_append(const Point &s);
  Mesh &S_append(const Matrix3double &S);
  Mesh &TV_append(const Matrix3int &TV);

  Dart find_path_direction(const Dart &d0, const Point &s1,
                           const int v1 = -1) const;
  Dart find_path_direction(const Point &s0, const Point &s1,
                           const Dart &d0) const;
  DartPair trace_path(const Dart &d0, const Point &s1, const int v1 = -1,
                      DartList *trace = NULL) const;
  DartPair trace_path(const Point &s0, const Point &s1, const Dart &d0,
                      DartList *trace = NULL) const;

  Dart find_path_direction_new(const Dart &d0, const Point &s1,
                           const int v1 = -1) const;
  Dart find_path_direction_new(const Point &s0, const Point &s1,
                           const Dart &d0) const;
  DartPair trace_path_new(const Dart &d0, const Point &s1, const int v1 = -1,
                      DartList *trace = NULL) const;
  DartPair trace_path_new(const Point &s0, const Point &s1, const Dart &d0,
                      DartList *trace = NULL) const;

  Dart locate_point(const Dart &d0, const Point &s, const int v = -1) const;
  Dart locate_vertex(const Dart &d0, const int v) const;

  Dart swapEdge(const Dart &d);
  Dart splitEdge(const Dart &d, int v);
  Dart splitTriangle(const Dart &d, int v);

  Mesh &removeLastVertex() { /* Does not check that the vertex is unused! */
    if (nV() > 0) {
      S_.rows(nV() - 1);
      if (use_VT_)
        VT_mapping_.pop_back();
    }
    return *this;
  };
  Mesh &unlinkEdge(const Dart &d);
  Mesh &unlinkTriangle(const int t);
  Mesh &relocateTriangle(const int t_source, const int t_target);
  int removeTriangle(const int t);

  Mesh &quad_tesselate(const Mesh &M);
  Mesh &make_globe(int subsegments, double radius);

  /* Traits: */
  double edgeLength(const Point &s0, const Point &s1) const;
  double edgeLength(const Dart &d) const;
  void barycentric(const Dart &d, const Point &s, Point &bary) const;
  void triangleBoundingBox(const Point &s0, const Point &s1, const Point &s2,
                           Point &mini, Point &maxi) const;
  double triangleArea(const Point &s0, const Point &s1, const Point &s2) const;
  double triangleCircumcircleRadius(const Point &s0,
                                    const Point &s1,
                                    const Point &s2) const;
  double edgeIntersection(const Point &s00, const Point &s01, const Point &s10,
                          const Point &s11, Point &c) const;

  void triangleBoundingBox(int t, Point &mini, Point &maxi) const;
  double triangleArea(int t) const;
  void triangleCircumcenter(int t, Point &c) const;
  double triangleCircumcircleRadius(int t) const;
  bool triangleEdgeLengths(int t, Point &len) const;
  int triangleEdgeLengthsArgMin(int t, Point &len) const;
  int triangleEdgeLengthsArgMax(int t, Point &len) const;
  double triangleLongestEdge(int t) const;
  double triangleShortestEdge(int t) const;
  double edgeEncroached(const Dart &d, const Point &s) const;

  /*!
    \brief Compute dart half-space test for a point.

    Positive if s is to the left of the edge defined by d.
   */
  double inLeftHalfspace(const Point &s0, const Point &s1,
                         const Point &s) const;

  /*!
    \brief Calculate FEM matrices.
   */
  void calcQblocks(SparseMatrix<double> &C0, SparseMatrix<double> &C1,
                   SparseMatrix<double> &G1, SparseMatrix<double> &B1,
                   Matrix<double> &Tareas) const;
  void calcQblocksAni(SparseMatrix<double> &G1, const Matrix<double> &gamma,
                      const Matrix<double> &vec) const;
  std::vector<SparseMatrix<double>> calcGradientMatrices() const;

  // No need for IOHeader and IOHelper classes when using Rcpp
#ifndef FMESHER_WITH_R
  /*! \brief Store the mesh in files. */
  bool save(std::string filename_s, std::string filename_tv,
            bool binary = true) const;
  /*! \brief Read a mesh from files. */
  bool load(std::string filename_s, std::string filename_tv,
            bool binary = true);
  /*! \brief Store the mesh in files in old headerless ascii format. */
  bool save_ascii_2009(std::string filename_s, std::string filename_tv) const;
  /*! \brief Read a mesh from files in old headerless ascii format. */
  bool load_ascii_2009(std::string filename_s, std::string filename_tv);
#endif // not FMESHER_WITH_R
};

std::unique_ptr<Matrix<double>> make_globe_points(int subsegments, double radius);

class MOAint {
  friend std::ostream &operator<<(std::ostream &output, const MOAint &MO);

private:
  size_t n_;
  const Matrix1int &M_;

public:
  MOAint(const Matrix1int &M, size_t n) : n_(n), M_(M){};
};

class MOAVTMap {
  friend std::ostream &operator<<(std::ostream &output, const MOAVTMap &MO);

private:
  size_t n_;
  const VTMapT &M_;

public:
  MOAVTMap(const VTMapT &M, size_t n) : n_(n), M_(M){};
};

class MOAVTMapV {
  friend std::ostream &operator<<(std::ostream &output, const MOAVTMapV &MO);

private:
  int v_;
  const VTMapT &M_;

public:
  MOAVTMapV(const VTMapT &M, int v) : v_(v), M_(M) {};
};

class MOAint3 {
  friend std::ostream &operator<<(std::ostream &output, const MOAint3 &MO);

private:
  size_t n_;
  const Matrix3int &M_;

public:
  MOAint3(const Matrix3int &M, size_t n) : n_(n), M_(M){};
};

class MOAdouble3 {
  friend std::ostream &operator<<(std::ostream &output, const MOAdouble3 &MO);

private:
  size_t n_;
  const Matrix3double &M_;

public:
  MOAdouble3(const Matrix3double &M, size_t n) : n_(n), M_(M){};
};

std::ostream &operator<<(std::ostream &output, const Point &MO);

/*! \brief Darts */
class Dart {
  friend std::ostream &operator<<(std::ostream &output, const Dart &d);

private:
  const Mesh *M_;
  size_t vi_;
  int edir_;
  int t_;

public:
  Dart(void) : M_(NULL), vi_(0), edir_(1), t_(0){};
  Dart(const Mesh &M, int t = 0, int edir = 1, size_t vi = 0)
    : M_(&M), vi_(vi), edir_(edir), t_(t){};
  Dart(const Dart &d) : M_(d.M_), vi_(d.vi_), edir_(d.edir_), t_(d.t_){};
  Dart &operator=(const Dart &d) {
    M_ = d.M_;
    vi_ = d.vi_;
    edir_ = d.edir_;
    t_ = d.t_;
    return *this;
  };

  const Mesh *M() const { return M_; };
  int vi() const { return vi_; };
  int edir() const { return edir_; };
  int t() const { return t_; };
  int v() const {
    if (!M_)
      return -1;
    else
      return M_->TV_[t_][vi_];
  };
  /* Opposite vertex; alpha0().v() */
  int vo() const {
    if (!M_)
      return -1;
    else
      return M_->TV_[t_][(vi_ + (3 + edir_)) % 3];
  };
  /* Adjacent triangle; alpha2().t() */
  int tadj() const {
    if (!M_)
      return -1;
    else
      return M_->TT_[t_][(vi_ + (3 - edir_)) % 3];
  };

  bool isnull() const { return (!M_); };
  bool operator==(const Dart &d) const {
    return ((d.t_ == t_) && (d.vi_ == vi_) && (d.edir_ == edir_));
  };
  bool operator<(const Dart &d) const {
    /* TODO: Add debug check for M_==d.M_ */
    return ((d.t_ < t_) ||
            ((d.t_ == t_) &&
             ((d.edir_ < edir_) || ((d.edir_ == edir_) && (d.vi_ < vi_)))));
  };
  bool operator!=(const Dart &d) const { return !(d == *this); };

  bool onBoundary() const { return (M_->TT_[t_][(vi_ + (3 - edir_)) % 3] < 0); }

  double inLeftHalfspace(const Point &s) const;
  double inCircumcircle(const Point &s) const;
  bool circumcircleOK(void) const;

  bool isSwapable() const {
    if (onBoundary())
      return false; /* Not swapable. */
    Dart dh(*this);
    const Point &s00 = M_->S_[dh.v()];
    dh.orbit2();
    const Point &s01 = M_->S_[dh.v()];
    dh.orbit2();
    const Point &s10 = M_->S_[dh.v()];
    dh.orbit2().orbit0rev().orbit2();
    const Point &s11 = M_->S_[dh.v()];
    /* Do both diagonals cross? Swapable. */
    return (((M_->inLeftHalfspace(s00, s01, s10) *
              M_->inLeftHalfspace(s00, s01, s11)) < 0.0) &&
            ((M_->inLeftHalfspace(s10, s11, s00) *
              M_->inLeftHalfspace(s10, s11, s01)) < 0.0));
  };

  bool isSwapableD() const { return (!circumcircleOK()); };

  /* Graph traversal algebra. */
  Dart &alpha0(void);
  Dart &alpha1(void);
  Dart &alpha2(void);
  Dart &orbit0(void);
  Dart &orbit1(void);
  Dart &orbit2(void);
  Dart &orbit0rev(void);
  Dart &orbit1rev(void);
  Dart &orbit2rev(void);
};

} /* namespace fmesh */

#endif
