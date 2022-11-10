#ifndef MESH_POINT_H
#define MESH_POINT_H

#include "assertm.h"
#include "point.h"
#include "vector.h"

class Point;
class Vector;

class MeshPoint : public Point {
 private:
  numeric _x, _y, _z;

 public:
  MeshPoint(const numeric x, const numeric y, const numeric z);
  MeshPoint(const Point &A);
  MeshPoint(const Point &A, const Vector &u);
  MeshPoint(const MeshPoint &A);
  MeshPoint(const MeshPoint *const A);
  MeshPoint(const MeshPoint &A, const Vector &u);
  MeshPoint(const MeshPoint *const A, const Vector &u);
  MeshPoint() = delete;

  numeric x() const;
  numeric y() const;
  numeric z() const;

  friend std::ostream &operator<<(std::ostream &os, const MeshPoint &A) {
    os << '[' << A._x << ',' << A._y << ',' << A._z << ']';
    return os;
  }
  friend std::ostream &operator<<(std::ostream &os, const MeshPoint *const A) {
    os << '[' << A->_x << ',' << A->_y << ',' << A->_z << ']';
    return os;
  }
  friend Point operator+(const MeshPoint &A, const MeshPoint &B) {
    return Point(A._x + B._x, A._y + B._y, A._z + B._z);
  }
  friend bool operator==(const MeshPoint &A, const MeshPoint &B) {
    auto diff =
        GiNaC::sqrt(GiNaC::pow(B._x - A._x, 2) + GiNaC::pow(B._y - A._y, 2) +
                    GiNaC::pow(B._z - A._z, 2));
    return diff < eps;
  }
  friend bool operator!=(const MeshPoint &A, const MeshPoint &B) {
    return !(A == B);
  }
};

#endif
