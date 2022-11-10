#include "face.h"

Face::Face(HalfEdgeIndex halfedge, Triangle triangle)
    : _halfedge(halfedge), _triangle(triangle) {}

Triangle Face::triangle() const { return _triangle; }

HalfEdgeIndex Face::halfedge_index() const { return _halfedge; }