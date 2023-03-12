#pragma once

#include "assertm.h"
#include "constants.h"
#include "point.h"
#include "vector.h"

using HalfEdgeIndex = int;

class MeshPoint : public Point {
 private:
  Point _point;
  // Outgoing halfedge of boundary vertice must be boundary edge
  std::vector<HalfEdgeIndex> _outgoing;
  MeshPointIndex _self;
  MeshPoint(const Point &A, HalfEdgeIndex outgoing = kInvalidEdgeIndex);

 public:
  MeshPoint(const numeric x, const numeric y, const numeric z,
            HalfEdgeIndex outgoing = kInvalidEdgeIndex);
  MeshPoint(const Point &A, const Vector &u,
            HalfEdgeIndex outgoing = kInvalidEdgeIndex);
  MeshPoint(const MeshPoint &A, const Vector &u,
            HalfEdgeIndex outgoing = kInvalidEdgeIndex);
  MeshPoint() = default;
  MeshPoint(const MeshPoint &A) = default;

  template <std::size_t Index>
  std::tuple_element_t<Index, MeshPoint> get() const {
    if constexpr (Index == 0) return _point.x();
    if constexpr (Index == 1) return _point.y();
    if constexpr (Index == 2) return _point.z();
  }

  numeric x() const;
  numeric y() const;
  numeric z() const;

  Point get_point() const;
  std::vector<HalfEdgeIndex> get_outgoing() const;
  MeshPointIndex get_index() const;
  void add_outgoing(HalfEdgeIndex outgoing);
  bool has_outgoing(HalfEdgeIndex outgoing) const;
  void set_index(MeshPointIndex index);
  static MeshPoint create_from_point(
      const Point A, const HalfEdgeIndex outgoing = kInvalidEdgeIndex) {
    MeshPoint P(A.x(), A.y(), A.z(), outgoing);
    return P;
  }

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
    auto diff = GiNaC::pow(B._point.x() - A._point.x(), 2) +
                GiNaC::pow(B._point.y() - A._point.y(), 2) +
                GiNaC::pow(B._point.z() - A._point.z(), 2);
    return diff < kEps * kEps;
  }
  friend bool operator!=(const MeshPoint &A, const MeshPoint &B) {
    return !(A == B);
  }
};

namespace std {
template <>
struct tuple_size<::MeshPoint> {
  static constexpr size_t value = 3;
};
template <>
struct tuple_element<0, ::MeshPoint> {
  using type = numeric;
};
template <>
struct tuple_element<1, ::MeshPoint> {
  using type = numeric;
};
template <>
struct tuple_element<2, ::MeshPoint> {
  using type = numeric;
};
}  // namespace std
