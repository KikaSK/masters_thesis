#include "mesh_point.h"

MeshPoint::MeshPoint(const numeric x, const numeric y, const numeric z,
                     HalfEdgeIndex outgoing)
    : _point(Point(x, y, z)), _outgoing(outgoing){};
MeshPoint::MeshPoint(const Point &A, HalfEdgeIndex outgoing)
    : _point(A), _outgoing(outgoing){};
MeshPoint::MeshPoint(const Point &A, const Vector &u, HalfEdgeIndex outgoing)
    : _point(Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z())),
      _outgoing(outgoing){};
MeshPoint::MeshPoint(const MeshPoint &A, const Vector &u,
                     HalfEdgeIndex outgoing)
    : _point(Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z())),
      _outgoing(outgoing){};

numeric MeshPoint::x() const { return _point.x(); }
numeric MeshPoint::y() const { return _point.y(); }
numeric MeshPoint::z() const { return _point.z(); }
