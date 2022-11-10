#include "point.h"

Point::Point(const numeric x, const numeric y, const numeric z)
    : _x(x), _y(y), _z(z){};
// copy constructor
Point::Point(const Point &A) : _x(A.x()), _y(A.y()), _z(A.z()){};
Point::Point(const Point &A, const Vector &u)
    : _x(A.x() + u.x()), _y(A.y() + u.y()), _z(A.z() + u.z()){};

numeric Point::x() const { return _x; }
numeric Point::y() const { return _y; }
numeric Point::z() const { return _z; }
