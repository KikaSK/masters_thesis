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

  AB.set_active();
  BC.set_active();
  CA.set_active();
  
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
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in is_active function!");
  return _mesh_edges[index].is_active();
}
bool Mesh::is_checked(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in is_cheked function!");
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
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_halfedge function!");
  return _mesh_edges[index];
}
MeshPoint Mesh::get_meshpoint(MeshPointIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_meshpoint function!");
  return _mesh_points[index];
}
Face Mesh::get_face(FaceIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_face function!");
  return _mesh_triangles[index];
}

HalfEdge Mesh::get_previous_halfedge(const HalfEdge &halfedge) const {
  auto i_previous = halfedge.get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  assertm(i_previous < _mesh_edges.size() && i_previous >= 0,
          "Invalid previous halfedge index in get_previous_halfedge function!");
  return _mesh_edges[i_previous];
}
HalfEdge Mesh::get_previous_halfedge(HalfEdgeIndex index) const {
  auto i_previous = _mesh_edges[index].get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  assertm(i_previous < _mesh_edges.size() && i_previous >= 0,
          "Invalid previous halfedge index in get_previous_halfedge function!");
  return _mesh_edges[i_previous];
}
HalfEdge Mesh::get_next_halfedge(const HalfEdge &halfedge) const {
  auto i_next = halfedge.get_next();
  assertm(i_next != -1, "No next halfedge");
  assertm(i_next < _mesh_edges.size() && i_next >= 0,
          "Invalid next halfedge index in get_next_halfedge function!");
  return _mesh_edges[i_next];
}
HalfEdge Mesh::get_next_halfedge(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_next_halfedge function!");
  auto i_next = _mesh_edges[index].get_next();
  assertm(i_next != -1, "No next halfedge");
  assertm(i_next < _mesh_edges.size() && i_next >= 0,
          "Invalid next halfedge index in get_next_halfedge function!");
  return _mesh_edges[i_next];
}
std::optional<HalfEdge> Mesh::get_opposite_halfedge(
    const HalfEdge &halfedge) const {
  auto i_opposite = halfedge.get_opposite();
  if (i_opposite == -1) return std::nullopt;
  assertm(i_opposite < _mesh_edges.size() && i_opposite >= 0,
          "Invalid halfedge index in get_opposite_halfedge function!");
  return _mesh_edges[i_opposite];
}
std::optional<HalfEdge> Mesh::get_opposite_halfedge(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_opposite_halfedge function!");
  auto i_opposite = _mesh_edges[index].get_opposite();
  if (i_opposite == -1) return std::nullopt;
  assertm(i_opposite < _mesh_edges.size() && i_opposite >= 0,
          "Invalid halfedge index in get_opposite_halfedge function!");
  return _mesh_edges[i_opposite];
}

HalfEdgeIndex Mesh::get_previous_index(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_previous_index function!");
  auto i_previous = _mesh_edges[index].get_previous();
  assertm(i_previous != -1, "No previous halfedge");
  assertm(i_previous < _mesh_edges.size() && i_previous >= 0,
          "Invalid previous halfedge index in get_previous_index function!");
  return i_previous;
}
HalfEdgeIndex Mesh::get_next_index(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in get_next_index function!");
  auto i_next = _mesh_edges[index].get_next();
  assertm(i_next != -1, "No next halfedge");
  assertm(i_next < _mesh_edges.size() && i_next >= 0,
          "Invalid next halfedge index in get_next_index function!");
  return i_next;
}
HalfEdgeIndex Mesh::get_opposite_index(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size(),
          "Invalid halfedge index in get_opposite_index function!");
  return _mesh_edges[index].get_opposite();
}

unordered_set<HalfEdge> Mesh::get_active_edges() const { return _active_edges; }

bool Mesh::has_active_edge(const HalfEdge &halfedge) const {
  return _active_edges.find(halfedge) != _active_edges.end();
}

bool Mesh::has_active_edge(const HalfEdgeIndex &index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in has_active_edge function!");
  return _active_edges.find(_mesh_edges[index]) != _active_edges.end();
}

size_t Mesh::get_active_edges_size() const { return _active_edges.size(); }

void Mesh::_bound_consecutive(HalfEdge *previous,
                              const HalfEdgeIndex i_previous, HalfEdge *next,
                              const HalfEdgeIndex i_next) const {
  previous->set_next(i_next);
  next->set_previous(i_previous);
}

void Mesh::_bound_opposite(HalfEdge *edge1, const HalfEdgeIndex i_edge1,
                           HalfEdge *edge2, const HalfEdgeIndex i_edge2) {
  edge1->set_opposite(i_edge2);
  edge2->set_opposite(i_edge1);
  edge1->set_inside();
  edge2->set_inside();
  _active_edges.erase(*edge1);
  _active_edges.erase(*edge2);
}

void Mesh::_bound_face(const FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                       HalfEdge *edge3) const {
  edge1->set_incident(i_face);
  edge2->set_incident(i_face);
  edge3->set_incident(i_face);
}

// having edge BA, we are trying to bound it with edge AB if it exists, setting
// AB inside edge
void Mesh::_bound_opposite_outgoing(const MeshPoint &A,
                                    const MeshPointIndex i_B,
                                    const HalfEdgeIndex i_BA) {
  const std::vector<HalfEdgeIndex> A_outgoing =
      A.get_outgoing();  // all edges outgoing from A
  for (auto outgoing : A_outgoing) {
    // std::cout << outgoing << endl;
    assertm(outgoing < _mesh_edges.size(),
            "Invalid halfedge index in _bound_opposite_outgoing function!");
    if (_mesh_edges[outgoing].get_B() /*second point of outgoing edge*/ ==
        i_B)  // we found edge AB
    {
      _bound_opposite(&(_mesh_edges[outgoing]), outgoing, &(_mesh_edges[i_BA]),
                      i_BA);
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

  if (type != "new") assertm(i_P >= 0, "Invalid index of given meshpoint.");

  HalfEdgeIndex i_edge = _mesh_edges.size();
  MeshPointIndex i_point = _mesh_points.size();
  FaceIndex i_face = _mesh_triangles.size();

  assertm(i_AB < _mesh_edges.size(),
          "Invalid halfedge index in add_triangle function!");
  HalfEdge &AB = _mesh_edges[i_AB];
  MeshPointIndex i_A = AB.get_A();
  MeshPointIndex i_B = AB.get_B();
  MeshPoint &A = _mesh_points[i_A];
  MeshPoint &B = _mesh_points[i_B];

  Edge e_BA(AB.get_point_B(), AB.get_point_A());
  HalfEdge BA(e_BA, i_B, i_A);  // i_edge
  HalfEdgeIndex i_BA = i_edge;
  _bound_opposite(&AB, i_AB, &BA, i_BA);

  assertm(AB.is_inside() && BA.is_inside(),
          "The edge is not set to inside edge!");

  MeshPoint new_meshpoint(new_point);
  if (type == "new") {
    _mesh_points.push_back(new_meshpoint);
    i_P = i_point;
  }

  MeshPoint &P = _mesh_points[i_P];

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
  std::cout << "Adding outgoing " << i_AP << " to " << A << endl
            << i_PB << " to " << P << ", " << i_BA << " to " << B << endl;
  std::cout << "current state of " << P << endl;
  for (auto o : P.get_outgoing()) std::cout << o << endl;
  */
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

  _mesh_edges.push_back(BA);
  _mesh_edges.push_back(AP);
  _mesh_edges.push_back(PB);

  if (type == "fill") {  // the edges are inside and we need to bound them with
                         // opposite edges
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
  } else if (type == "previous") {           // only edge PB is active edge
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(_mesh_edges[i_PB]);
  } else if (type == "next") {               // only edge AP is active edge
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
    _mesh_edges[i_AP].set_active();
    _active_edges.insert(_mesh_edges[i_AP]);
  } else if (type == "new" ||
             type == "overlap") {  // both AP and PB are new active edges
    _mesh_edges[i_AP].set_active();
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(_mesh_edges[i_PB]);
    _active_edges.insert(_mesh_edges[i_AP]);
  }
  /*
    Insert the new faces to lists
  */
  _mesh_triangles.push_back(F);

  if (type == "new") _mesh_points.push_back(P);
}