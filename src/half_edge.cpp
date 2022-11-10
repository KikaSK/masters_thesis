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