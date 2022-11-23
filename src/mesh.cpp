#include "mesh.h"

Mesh::Mesh(Face F) {
  F.set_halfedge(0);
  _mesh_triangles.push_back(F);
  MeshPoint A(F.get_triangle().A(), 0);  // the outgoing halfedge will be AB
  MeshPoint B(F.get_triangle().B(), 1);  // the outgoing halfedge will be BC
  MeshPoint C(F.get_triangle().C(), 2);  // the outgoing halfedge will be CA

  _mesh_points.push_back(A);
  _mesh_points.push_back(B);
  _mesh_points.push_back(C);

  HalfEdge AB(F.get_triangle().AB(),
              0,   // A index is 0
              1,   // B index is 1
              2,   // previous is CA
              1,   // next is BC
              -1,  // no opposite
              0    // incident face is F
  );
  HalfEdge BC(F.get_triangle().BC(),
              1,   // B index is 1
              2,   // C index is 2
              0,   // previous is AB
              2,   // next is CA
              -1,  // no opposite
              0    // incident face is F
  );
  HalfEdge CA(F.get_triangle().CA(),
              2,   // C index is 2
              0,   // A index is 0
              1,   // previous is BC
              0,   // next is AB
              -1,  // no opposite
              0    // incident face is F
  );

  _mesh_edges.push_back(AB);
  _mesh_edges.push_back(BC);
  _mesh_edges.push_back(CA);

  _active_edges.insert(AB);
  _active_edges.insert(BC);
  _active_edges.insert(CA);
}

void Mesh::cout_triangles() const {
  for (auto T : _mesh_triangles) {
    std::cout << T.get_triangle() << endl;
  }
}

bool Mesh::is_active(HalfEdgeIndex index) const {
  return _mesh_edges[index].is_active();
}
bool Mesh::is_checked(HalfEdgeIndex index) const {
  return _mesh_edges[index].is_checked();
  //_checked_edges.find(_mesh_edges[index]) != _checked_edges.end();
}
bool Mesh::is_boundary(HalfEdgeIndex index) const {
  return is_active(index) || is_checked(index);
}
bool Mesh::is_in_mesh(const HalfEdge &halfedge) const {
  for (auto edge : _mesh_edges) {
    if (edge == halfedge) {
      return true;
    }
  }
  return false;
}

HalfEdge Mesh::get_halfedge(HalfEdgeIndex index) const {
  return _mesh_edges[index];
}
MeshPoint Mesh::get_meshpoint(MeshPointIndex index) const {
  return _mesh_points[index];
}
Face Mesh::get_face(FaceIndex index) const { return _mesh_triangles[index]; }

HalfEdge Mesh::get_previous_halfedge(const HalfEdge &halfedge) const {
  auto i_previous = halfedge.get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  return _mesh_edges[i_previous];
}
HalfEdge Mesh::get_previous_halfedge(HalfEdgeIndex index) const {
  auto i_previous = _mesh_edges[index].get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  return _mesh_edges[i_previous];
}
HalfEdge Mesh::get_next_halfedge(const HalfEdge &halfedge) const {
  auto i_next = halfedge.get_next();
  assertm(i_next != -1, "No next halfedge");
  return _mesh_edges[i_next];
}
HalfEdge Mesh::get_next_halfedge(HalfEdgeIndex index) const {
  auto i_next = _mesh_edges[index].get_next();
  assertm(i_next != -1, "No next halfedge");
  return _mesh_edges[i_next];
}
std::optional<HalfEdge> Mesh::get_opposite_halfedge(
    const HalfEdge &halfedge) const {
  auto i_opposite = halfedge.get_opposite();
  if (i_opposite == -1) return std::nullopt;
  return _mesh_edges[i_opposite];
}
std::optional<HalfEdge> Mesh::get_opposite_halfedge(HalfEdgeIndex index) const {
  auto i_opposite = _mesh_edges[index].get_opposite();
  if (i_opposite == -1) return std::nullopt;
  return _mesh_edges[i_opposite];
}

HalfEdgeIndex Mesh::get_previous_index(HalfEdgeIndex index) const {
  auto i_previous = _mesh_edges[index].get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  return i_previous;
}
HalfEdgeIndex Mesh::get_next_index(HalfEdgeIndex index) const {
  auto i_next = _mesh_edges[index].get_next();
  assertm(i_next != -1, "No next halfedge");
  return i_next;
}
HalfEdgeIndex Mesh::get_opposite_index(HalfEdgeIndex index) const {
  return _mesh_edges[index].get_opposite();
}

unordered_set<HalfEdge> Mesh::get_active_edges() const { return _active_edges; }

bool Mesh::has_active_edge(const HalfEdge &halfedge) const {
  return _active_edges.find(halfedge) != _active_edges.end();
}

bool Mesh::has_active_edge(const HalfEdgeIndex &index) const {
  assertm(index < _mesh_edges.size(), "Invalid halfedge index!");
  return _active_edges.find(_mesh_edges[index]) != _active_edges.end();
}

size_t Mesh::get_active_edges_size() const { return _active_edges.size(); }

void Mesh::_bound_consecutive(HalfEdge *previous, HalfEdgeIndex i_previous,
                              HalfEdge *next, HalfEdgeIndex i_next) const {
  previous->set_next(i_next);
  next->set_previous(i_previous);
}

void Mesh::_bound_opposite(HalfEdge *edge1, HalfEdgeIndex i_edge1,
                           HalfEdge *edge2, HalfEdgeIndex i_edge2) {
  edge1->set_opposite(i_edge2);
  edge2->set_opposite(i_edge1);
  edge1->set_inside();
  edge2->set_inside();
  _active_edges.erase(*edge1);
  _active_edges.erase(*edge2);
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

  HalfEdge &AB = _mesh_edges[i_AB];
  MeshPointIndex i_A = AB.get_A();
  MeshPointIndex i_B = AB.get_B();
  MeshPoint &A = _mesh_points[i_A];
  MeshPoint &B = _mesh_points[i_B];

  Edge e_BA(AB.get_point_B(), AB.get_point_A());
  HalfEdge BA(e_BA, i_B, i_A);  // i_edge
  HalfEdgeIndex i_BA = i_edge;
  _bound_opposite(&AB, i_AB, &BA, i_BA);

  std::cout << (AB.is_inside() ? "Correctly set" : "Incorrectly set")
            << std::endl;

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
  Face F(ABP, i_AP);
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

  std::cout << ((PB.is_active() && AP.is_active()) ? "correct" : "incorrect")
            << std::endl;
  _mesh_points.push_back(P);
  _mesh_triangles.push_back(F);
}

void Mesh::_assert_halfedge_index(const HalfEdgeIndex &index) const {
  assertm(index < _mesh_edges.size(), "Invalid halfedge index!");
}
void Mesh::_assert_meshpoint_index(const MeshPointIndex &index) const {
  assertm(index < _mesh_points.size(), "Invalid meshpoint index!");
}
void Mesh::_assert_face_index(const FaceIndex &index) const {
  assertm(index < _mesh_triangles.size(), "Invalid face index!");
}