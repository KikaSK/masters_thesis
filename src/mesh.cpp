#include "mesh.h"

Mesh::Mesh(Face F) {
  _mesh_triangles.push_back(F);
  MeshPoint A(F.triangle().A(), 0);  // the outgoing halfedge will be AB
  MeshPoint B(F.triangle().B(), 1);  // the outgoing halfedge will be BC
  MeshPoint C(F.triangle().C(), 2);  // the outgoing halfedge will be CA

  _mesh_points.push_back(A);
  _mesh_points.push_back(B);
  _mesh_points.push_back(C);

  HalfEdge AB(F.triangle().AB(),
              0,   // A index is 0
              1,   // B index is 1
              2,   // previous is CA
              1,   // next is BC
              -1,  // no opposite
              0    // incident face is F
  );
  HalfEdge BC(F.triangle().BC(),
              1,   // B index is 1
              2,   // C index is 2
              0,   // previous is AB
              2,   // next is CA
              -1,  // no opposite
              0    // incident face is F
  );
  HalfEdge CA(F.triangle().CA(),
              2,   // C index is 2
              0,   // A index is 0
              1,   // previous is BC
              0,   // next is AB
              -1,  // no opposite
              0    // incident face is F
  );
}

void Mesh::cout_triangles() const {
  for (auto T : _mesh_triangles) {
    std::cout << T.triangle() << endl;
  }
}

bool Mesh::is_active(HalfEdgeIndex index) const {
  return _active_edges.find(_mesh_edges[index]) != _active_edges.end();
}
bool Mesh::is_checked(HalfEdgeIndex index) const {
  return _checked_edges.find(_mesh_edges[index]) != _checked_edges.end();
}
bool Mesh::is_boundary(HalfEdgeIndex index) const {
  return is_active(index) || is_checked(index);
}
bool Mesh::is_active(HalfEdge halfedge) const {
  return _active_edges.find(halfedge) != _active_edges.end();
}
bool Mesh::is_checked(HalfEdge halfedge) const {
  return _checked_edges.find(halfedge) != _checked_edges.end();
}
bool Mesh::is_boundary(HalfEdge halfedge) const {
  return is_active(halfedge) || is_checked(halfedge);
}
bool Mesh::is_in_mesh(HalfEdge halfedge) const {
  for (auto edge : _mesh_edges) {
    if (edge == halfedge) {
      return true;
    }
  }
  return false;
}

void Mesh::_bound_consecutive(HalfEdge *previous, HalfEdgeIndex i_previous,
                              HalfEdge *next, HalfEdgeIndex i_next) const {
  previous->set_next(i_next);
  next->set_previous(i_previous);
}

void Mesh::_bound_opposite(HalfEdge *edge1, HalfEdgeIndex i_edge1,
                           HalfEdge *edge2, HalfEdgeIndex i_edge2) const {
  edge1->set_opposite(i_edge2);
  edge2->set_opposite(i_edge1);
  edge1->set_inside();
  edge2->set_inside();
}

void Mesh::_bound_face(FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                       HalfEdge *edge3) const {
  edge1->set_incident(i_face);
  edge2->set_incident(i_face);
  edge3->set_incident(i_face);
}

// having edge BA, we are trying to bound it with edge AB if it exists, setting
// AB inside edge
void Mesh::_bound_opposite_outgoing(const MeshPoint &A, MeshPointIndex i_B,
                                    HalfEdge *BA, HalfEdgeIndex i_BA) {
  const std::vector<HalfEdgeIndex> A_outgoing =
      A.get_outgoing();  // all edges outgoing from A
  for (auto outgoing : A_outgoing) {
    if (_mesh_edges[outgoing].get_B() /*second point of outgoing edge*/ ==
        i_B)  // we found edge AB
    {
      _bound_opposite(&(_mesh_edges[outgoing]), outgoing, BA, i_BA);
    }
  }
}

// Adding triangle ABP where P is new_point
void Mesh::add_triangle(HalfEdgeIndex i_AB, Point new_point,
                        std::string type /*new, next, previous, overlap, fill*/,
                        MeshPointIndex i_P) {
  /*
  Get AB from the list of halfedges (edges[i_AB])
  i = edges.size
  j = points.size
  k = triangles.size
  Create halfedge BA (i)
  bound_opposite(BA(i), AB(i_AB))
  if P is new point
    Create meshpoint P (j)
  */

  HalfEdgeIndex i_edge = _mesh_edges.size();
  MeshPointIndex i_point = _mesh_points.size();
  FaceIndex i_face = _mesh_triangles.size();

  HalfEdge AB = _mesh_edges[i_AB];
  MeshPointIndex i_A = AB.get_A();
  MeshPointIndex i_B = AB.get_B();
  MeshPoint A = _mesh_points[i_A];
  MeshPoint B = _mesh_points[i_B];

  Edge e_BA(AB.get_point_B(), AB.get_point_A());
  HalfEdge BA(e_BA, i_B, i_A);  // i_edge
  HalfEdgeIndex i_BA = i_edge;
  _bound_opposite(&AB, i_AB, &BA, i_BA);

  MeshPoint P(new_point);
  if (type == "new")
    i_P = i_point;
  else
    P = _mesh_points[i_P];

  /*
  Create halfedge AP (i+1) (boundary edge)
  bound_consecutive(BA(i), AP(i+1)
  Create halfedge PB (i+2) (boundary edge)
  bound_consecutive(PB(i+2), BA(i))
  */

  Edge e_AP(AB.get_point_A(), new_point);
  HalfEdge AP(e_AP, AB.get_A(), i_P);
  HalfEdgeIndex i_AP = i_edge + 1;

  Edge e_PB(new_point, AB.get_point_B());
  HalfEdge PB(e_PB, i_P, i_B);
  HalfEdgeIndex i_PB = i_edge + 2;

  _bound_consecutive(&BA, i_BA, &AP, i_AP);
  _bound_consecutive(&AP, i_AP, &PB, i_PB);
  _bound_consecutive(&PB, i_PB, &BA, i_BA);

  A.add_outgoing(i_AP);
  P.add_outgoing(i_PB);
  B.add_outgoing(i_BA);

  /*
    create face F(triangle ABP, halfedge AP) (k)
    bound_face(F(k), AP(i+1), PB(i+2), BA(i));
  */

  Triangle ABP(A.get_point(), B.get_point(), P.get_point());
  Face F(i_AP, ABP);
  _bound_face(i_face, &BA, &AP, &PB);

  /*
    Find out which edges are the new active edges
  */

  if (type == "fill") {  // the edges are inside and we need to bound them with
                         // opposite edges
    _bound_opposite_outgoing(P, i_A, &AP, i_AP);  // try to bound AP with PA
    _bound_opposite_outgoing(B, i_P, &PB, i_PB);  // try to bound PB with BP
  } else if (type == "previous") {                // only edge PB is active edge
    _bound_opposite_outgoing(P, i_A, &AP, i_AP);  // try to bound AP with PA
    PB.set_active();
    _active_edges.insert(PB);
  } else if (type == "next") {                    // only edge AP is active edge
    _bound_opposite_outgoing(B, i_P, &PB, i_PB);  // try to bound PB with BP
    AP.set_active();
    _active_edges.insert(AP);
  } else if (type == "new" ||
             type == "overlap") {  // both AP and PB are new active edges
    AP.set_active();
    PB.set_active();
    _active_edges.insert(PB);
    _active_edges.insert(AP);
  }

  /*
    Insert the new edges, faces and vertices to lists
  */

  _mesh_edges.push_back(BA);
  _mesh_edges.push_back(AP);
  _mesh_edges.push_back(PB);

  _mesh_points.push_back(P);
  _mesh_triangles.push_back(F);
}