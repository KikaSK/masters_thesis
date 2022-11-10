#include "edge.h"

#include <exception>

Edge::Edge(Point A, Vector u) : _A(A) {
  assertm(!u.is_zero(), "Trying to make edge from the same points!");
  _B = Point(A.x() + u.x(), A.y() + u.y(), A.z() + u.z());
}

Edge::Edge(Point A, Point B) : _A(A), _B(B) {
  assertm(_A != _B, "Trying to make edge from the same points!");
}

Point Edge::A() const { return _A; }
Point Edge::B() const { return _B; }

numeric Edge::get_length() const {
  numeric length =
      sqrt(pow(_A.x() - _B.x(), numeric(2)) + pow(_A.y() - _B.y(), numeric(2)) +
           pow(_A.z() - _B.z(), numeric(2)));
  return length;
}

Point Edge::get_midpoint() const {
  Point mid((_A.x() + _B.x()) / 2, (_A.y() + _B.y()) / 2,
            (_A.z() + _B.z()) / 2);
  return mid;
}