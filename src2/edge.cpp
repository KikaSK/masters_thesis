#include "edge.h"

#include <exception>

Edge::Edge(MeshPoint *A, Vector u, Edge *previous, Edge *next, Edge *opposite)
    : _A(A), _previous(previous), _next(next), _opposite(opposite) {
  assertm(!u.is_zero(), "Trying to make edge from the same points!");
  MeshPoint B(A->x() + u.x(), A->y() + u.y(), A->z() + u.z());
  _B = &B;
}

Edge::Edge(MeshPoint *A, MeshPoint *B, Edge *previous, Edge *next,
           Edge *opposite)
    : _A(A), _B(B), _previous(previous), _next(next), _opposite(opposite) {
  assertm(*_A != *_B, "Trying to make edge from the same points!");
}

MeshPoint *Edge::A() const { return _A; }
MeshPoint *Edge::B() const { return _B; }

Edge *Edge::get_opposite() const { return _opposite; }
void Edge::set_opposite(Edge *const opposite) {
  _opposite = opposite;
  return;
}
Edge *Edge::get_previous() const { return _previous; }
void Edge::set_previous(Edge *const previous) {
  _previous = previous;
  return;
}
Edge *Edge::get_next() const { return _next; }
void Edge::set_next(Edge *const next) {
  _next = next;
  return;
}

numeric Edge::get_length() const {
  numeric length = sqrt(pow(_A->x() - _B->x(), numeric(2)) +
                        pow(_A->y() - _B->y(), numeric(2)) +
                        pow(_A->z() - _B->z(), numeric(2)));
  return length;
}

Point Edge::get_midpoint() const {
  Point mid((_A->x() + _B->x()) / 2, (_A->y() + _B->y()) / 2,
            (_A->z() + _B->z()) / 2);
  return mid;
}