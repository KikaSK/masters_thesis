#include "mesh_point.h"

MeshPoint::MeshPoint(const numeric x, const numeric y, const numeric z)
    : _x(x), _y(y), _z(z){};
MeshPoint::MeshPoint(const Point &A) : _x(A.x()), _y(A.y()), _z(A.z()){};
MeshPoint::MeshPoint(const Point &A, const Vector &u)
    : _x(A.x() + u.x()), _y(A.y() + u.y()), _z(A.z() + u.z()){};

MeshPoint::MeshPoint(const MeshPoint &A, const Vector &u)
    : _x(A.x() + u.x()), _y(A.y() + u.y()), _z(A.z() + u.z()){};
MeshPoint::MeshPoint(const MeshPoint *const A, const Vector &u)
    : _x(A->x() + u.x()), _y(A->y() + u.y()), _z(A->z() + u.z()){};
// copy constructor
MeshPoint::MeshPoint(const MeshPoint &A)
    : Point(A), _x(A.x()), _y(A.y()), _z(A.z()){};
MeshPoint::MeshPoint(const MeshPoint *const A)
    : _x(A->x()), _y(A->y()), _z(A->z()){};

numeric MeshPoint::x() const { return _x; }
numeric MeshPoint::y() const { return _y; }
numeric MeshPoint::z() const { return _z; }
