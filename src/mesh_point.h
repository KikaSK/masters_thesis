#ifndef MESH_POINT_H
#define MESH_POINT_H

#include "assertm.h"
#include "constants.h"
#include "point.h"
#include "vector.h"

using HalfEdgeIndex = int;

class Point;
class Vector;

class MeshPoint : public Point {
 private:
  Point _point;
  // Outgoing halfedge of boundary vertice must be boundary edge
  HalfEdgeIndex _outgoing;

 public:
  MeshPoint(const numeric x, const numeric y, const numeric z,
            HalfEdgeIndex outgoing = -1);
  MeshPoint(const Point &A, HalfEdgeIndex outgoing = -1);
  MeshPoint(const Point &A, const Vector &u, HalfEdgeIndex outgoing = -1);
  MeshPoint(const MeshPoint &A, const Vector &u, HalfEdgeIndex outgoing = -1);
  MeshPoint() = delete;

  numeric x() const;
  numeric y() const;
  numeric z() const;

  friend std::ostream &operator<<(std::ostream &os, const MeshPoint &A) {
    os << '[' << A._point.x() << ',' << A._point.y() << ',' << A._point.z()
       << ']';
    return os;
  }
  friend Point operator+(const MeshPoint &A, const MeshPoint &B) {
    return Point(A._point.x() + B._point.x(), A._point.y() + B._point.y(),
                 A._point.z() + B._point.z());
  }
  friend bool operator==(const MeshPoint &A, const MeshPoint &B) {
    auto diff = GiNaC::sqrt(GiNaC::pow(B._point.x() - A._point.x(), 2) +
                            GiNaC::pow(B._point.y() - A._point.y(), 2) +
                            GiNaC::pow(B._point.z() - A._point.z(), 2));
    return diff < k_eps;
  }
  friend bool operator!=(const MeshPoint &A, const MeshPoint &B) {
    return !(A == B);
  }
};

#endif
