#include "mesh_point.h"

MeshPoint::MeshPoint(const numeric x, const numeric y, const numeric z,
                     HalfEdgeIndex outgoing)
    : _point(Point(x, y, z)) {
  if (outgoing != kInvalidEdgeIndex) _outgoing.push_back(outgoing);
};
MeshPoint::MeshPoint(const Point &A, HalfEdgeIndex outgoing) : _point(A) {
  if (outgoing != kInvalidEdgeIndex) _outgoing.push_back(outgoing);
};
MeshPoint::MeshPoint(const Point &A, const Vector &u, HalfEdgeIndex outgoing)
    : _point(Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z())) {
  if (outgoing != kInvalidEdgeIndex) _outgoing.push_back(outgoing);
};
MeshPoint::MeshPoint(const MeshPoint &A, const Vector &u,
                     HalfEdgeIndex outgoing)
    : _point(Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z())) {
  if (outgoing != kInvalidEdgeIndex) _outgoing.push_back(outgoing);
};

numeric MeshPoint::x() const { return _point.x(); }
numeric MeshPoint::y() const { return _point.y(); }
numeric MeshPoint::z() const { return _point.z(); }

Point MeshPoint::get_point() const { return _point; }

bool MeshPoint::has_outgoing(HalfEdgeIndex outgoing) const {
  for (auto edge : _outgoing)
    if (edge == outgoing) return true;

  return false;
}

void MeshPoint::add_outgoing(HalfEdgeIndex outgoing) {
  if (!has_outgoing(outgoing)) _outgoing.push_back(outgoing);
}

void MeshPoint::set_index(MeshPointIndex index) {
  assertm(index != kInvalidPointIndex, "Setting invalid index!");
  _self = index;
}

std::vector<HalfEdgeIndex> MeshPoint::get_outgoing() const { return _outgoing; }

MeshPointIndex MeshPoint::get_index() const { return _self; }