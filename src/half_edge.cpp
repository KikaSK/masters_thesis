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

MeshPointIndex HalfEdge::get_A() const { return _A; }
MeshPointIndex HalfEdge::get_B() const { return _B; }

HalfEdgeIndex HalfEdge::get_opposite() const { return _opposite; }
void HalfEdge::set_opposite(HalfEdgeIndex opposite) { _opposite = opposite; }
HalfEdgeIndex HalfEdge::get_previous() const { return _previous; }
void HalfEdge::set_previous(HalfEdgeIndex previous) { _previous = previous; }
HalfEdgeIndex HalfEdge::get_next() const { return _next; }
void HalfEdge::set_next(HalfEdgeIndex next) { _next = next; }
FaceIndex HalfEdge::get_incident() const { return _incident; }
void HalfEdge::set_incident(FaceIndex incident) { _incident = incident; }
HalfEdgeIndex HalfEdge::get_index() const { return _self; }
void HalfEdge::set_index(HalfEdgeIndex index) { _self = index; }

numeric HalfEdge::get_length() const { return _edge.get_length(); }
Point HalfEdge::get_midpoint() const { return _edge.get_midpoint(); }
Edge HalfEdge::get_edge() const { return _edge; }
Point HalfEdge::get_point_A() const { return _edge.A(); }
Point HalfEdge::get_point_B() const { return _edge.B(); }

void HalfEdge::set_active() {
  _is_active = true;
  _is_checked = false;
}
void HalfEdge::set_checked() {
  _is_active = false;
  _is_checked = true;
}
void HalfEdge::set_inside() {
  _is_active = false;
  _is_checked = false;
}

bool HalfEdge::is_active() const { return _is_active; }
bool HalfEdge::is_checked() const { return _is_checked; }
bool HalfEdge::is_boundary() const { return _is_active || _is_checked; }
bool HalfEdge::is_inside() const { return !is_boundary(); }