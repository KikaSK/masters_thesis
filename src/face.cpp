#include "face.h"

Face::Face(Triangle triangle, HalfEdgeIndex halfedge)
    : _triangle(triangle), _halfedge(halfedge) {}

Triangle Face::get_triangle() const { return _triangle; }

HalfEdgeIndex Face::get_halfedge() const { return _halfedge; }
void Face::set_halfedge(HalfEdgeIndex index) { _halfedge = index; }

Point Face::A() const { return _triangle.A(); }
Point Face::B() const { return _triangle.B(); }
Point Face::C() const { return _triangle.C(); }

Edge Face::AB() const { return _triangle.AB(); }
Edge Face::BC() const { return _triangle.BC(); }
Edge Face::CA() const { return _triangle.CA(); }

Point Face::get_gravity_center() const {
  return _triangle.get_gravity_center();
}
bool Face::is_triangle() const { return _triangle.is_triangle(); }
Point Face::get_circumcenter() const { return _triangle.get_circumcenter(); }
Vector Face::get_normal() const { return _triangle.get_normal(); }
bool Face::is_in_triangle(Point P) const { return _triangle.is_in_triangle(P); }