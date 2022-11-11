#include "half_edge.h"

HalfEdge::HalfEdge(Edge edge, MeshPointIndex A, MeshPointIndex B,
                   HalfEdgeIndex previous, HalfEdgeIndex next,
                   HalfEdgeIndex opposite, FaceIndex incident)
    : _edge(edge),
      _A(A),
      _B(B),
      _previous(previous),
      _next(next),
      _opposite(opposite),
      _incident(incident) {}

MeshPointIndex HalfEdge::A() const { return _A; }
MeshPointIndex HalfEdge::B() const { return _B; }

HalfEdgeIndex HalfEdge::get_opposite() const { return _opposite; }
void HalfEdge::set_opposite(HalfEdgeIndex opposite) { _opposite = opposite; }
HalfEdgeIndex HalfEdge::get_previous() const { return _previous; }
void HalfEdge::set_previous(HalfEdgeIndex previous) { _previous = previous; }
HalfEdgeIndex HalfEdge::get_next() const { return _next; }
void HalfEdge::set_next(HalfEdgeIndex next) { _next = next; }
FaceIndex HalfEdge::get_incident() const { return _incident; }
void HalfEdge::set_incident(FaceIndex incident) { _incident = incident; }

numeric HalfEdge::get_length() const { return _edge.get_length(); }
Point HalfEdge::get_midpoint() const { return _edge.get_midpoint(); }
