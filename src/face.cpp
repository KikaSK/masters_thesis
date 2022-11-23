#include "face.h"

Face::Face(Triangle triangle, HalfEdgeIndex halfedge)
    : _triangle(triangle), _halfedge(halfedge) {}

Triangle Face::get_triangle() const { return _triangle; }

HalfEdgeIndex Face::get_halfedge() const { return _halfedge; }
void Face::set_halfedge(HalfEdgeIndex index) { _halfedge = index; }