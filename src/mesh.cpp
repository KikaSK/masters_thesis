#include "mesh.h"

Mesh::Mesh(const Triangle &T) {
  Face F(T);
  F.set_halfedge(0);
  _mesh_triangles.push_back(F);
  MeshPoint A(F.get_triangle().A(), 0);  // the outgoing halfedge will be AB
  MeshPoint B(F.get_triangle().B(), 1);  // the outgoing halfedge will be BC
  MeshPoint C(F.get_triangle().C(), 2);  // the outgoing halfedge will be CA
  A.set_index(0);
  B.set_index(1);
  C.set_index(2);

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

  AB.set_index(0);
  BC.set_index(1);
  CA.set_index(2);

  AB.set_active();
  BC.set_active();
  CA.set_active();

  _mesh_edges.push_back(AB);
  _mesh_edges.push_back(BC);
  _mesh_edges.push_back(CA);

  _active_edges.insert(0);
  _active_edges.insert(1);
  _active_edges.insert(2);

  edges_check("in mesh constructor:");
  // cout_mesh();
}

void Mesh::cout_triangles() const {
  for (auto T : _mesh_triangles) {
    std::cout << T.get_triangle() << endl;
  }
}
void Mesh::cout_triangles_number() const {
  std::cout << "Number of faces: " << _mesh_triangles.size() << endl;
}
void Mesh::cout_mesh() const {
  std::cout << "Triangles: " << endl;
  for (int i = 0; i < _mesh_triangles.size(); ++i) {
    std::cout << "Triangle " << i << ": " << endl
              << _mesh_triangles[i].get_triangle();
    std::cout << "  halfedge: " << _mesh_triangles[i].get_halfedge();
  }

  std::cout << endl << "Edges: " << endl;
  for (int i = 0; i < _mesh_edges.size(); ++i) {
    std::cout << "Edge " << i << ": " << endl << _mesh_edges[i].get_edge();
    std::cout << (_mesh_edges[i].is_active() ? "  is active" : "  is inside")
              << endl;
    std::cout << "  meshpoint1: " << _mesh_edges[i].get_A()
              << "meshpoint2: " << _mesh_edges[i].get_B() << endl;
    std::cout << "  opposite: " << _mesh_edges[i].get_opposite() << endl;
    std::cout << "  previous: " << _mesh_edges[i].get_previous() << endl;
    std::cout << "  next: " << _mesh_edges[i].get_next() << endl;
    std::cout << "  incident face: " << _mesh_edges[i].get_incident() << endl;
  }

  std::cout << "Points: " << endl;
  for (int i = 0; i < _mesh_points.size(); ++i) {
    std::cout << "Point " << i << endl
              << _mesh_points[i].get_point() << ":" << endl;
    std::cout << "  outgoing: " << endl;
    for (auto e : _mesh_points[i].get_outgoing()) {
      std::cout << "    " << e << endl;
    }
  }
  std::cout << endl;
}

bool Mesh::is_active(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in is_active function!");
  assertm(_mesh_edges[index].is_active() == has_active_edge(index),
          "Active edges inconsistent!");
  return _mesh_edges[index].is_active();
}
bool Mesh::is_checked(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in is_cheked function!");
  assertm(_mesh_edges[index].is_checked() == has_checked_edge(index),
          "Checked edges inconsistent!");
  return _mesh_edges[index].is_checked();
}
bool Mesh::is_bounding(HalfEdgeIndex index) const {
  assertm(index < _mesh_edges.size() && index >= 0,
          "Invalid halfedge index in is_bounding function!");
  assertm(!_mesh_edges[index].is_bounding(), "Found bounding edge!");
  return _mesh_edges[index].is_bounding();
}
bool Mesh::is_boundary(HalfEdgeIndex index) const {
  return is_active(index) || is_checked(index) || is_bounding(index);
}
bool Mesh::is_in_mesh(const Edge &edge) const {
  for (auto mesh_edge : _mesh_edges) {
    // comparing only endpoints in respective order
    if (edge == mesh_edge.get_edge()) {
      return true;
    }
  }
  return false;
}
HalfEdgeIndex Mesh::get_edge_index(const Edge &edge) const {
  for (HalfEdge mesh_edge : _mesh_edges) {
    // comparing only endpoints in respective order
    if (edge == mesh_edge.get_edge()) {
      return mesh_edge.get_index();
    }
  }
  return kInvalidEdgeIndex;
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

FaceIndex Mesh::get_incident_face(const HalfEdgeIndex &index) const {
  return _mesh_edges[index].get_incident();
}

unordered_set<HalfEdgeIndex> Mesh::get_active_edges() const {
  return _active_edges;
}
unordered_set<HalfEdgeIndex> Mesh::get_checked_edges() const {
  return _checked_edges;
}
unordered_set<HalfEdgeIndex> Mesh::get_bounding_edges() const {
  return _bounding_edges;
}
vector<MeshPoint> Mesh::get_mesh_points() const { return _mesh_points; }
vector<HalfEdge> Mesh::get_mesh_edges() const { return _mesh_edges; }
vector<Face> Mesh::get_mesh_faces() const { return _mesh_triangles; }
HalfEdgeIndex Mesh::get_active_edge() const {
  assertm(!active_edges_empty(), "Picking edge from empty active edges list!");
  return *_active_edges.begin();
}
HalfEdgeIndex Mesh::get_checked_edge() const {
  assertm(!checked_edges_empty(),
          "Picking edge from empty checked edges list!");
  return *_checked_edges.begin();
}
void Mesh::remove_active_edge(const HalfEdgeIndex &index) {
  assertm(_active_edges.find(index) != _active_edges.end(),
          "Removing non-existent active edge!");
  _active_edges.erase(index);
  _mesh_edges[index].set_inside();
}
void Mesh::add_edge_to_active(const HalfEdgeIndex &index) {
  assertm(_active_edges.find(index) == _active_edges.end(),
          "Inserting active edge twice!");
  _active_edges.insert(index);
  _mesh_edges[index].set_active();
}
void Mesh::remove_checked_edge(const HalfEdgeIndex &index) {
  assertm(_checked_edges.find(index) != _checked_edges.end(),
          "Removing non-existent checked edge!");
  _checked_edges.erase(index);
  _mesh_edges[index].set_inside();
  assertm(!is_checked(index), "not checked edge checked!");
  return;
}
void Mesh::add_edge_to_checked(const HalfEdgeIndex &index) {
  assertm(_checked_edges.find(index) == _checked_edges.end(),
          "Inserting checked edge twice!");
  _checked_edges.insert(index);
  _mesh_edges[index].set_checked();
  assertm(is_checked(index), "checked edge not checked!");
  return;
}
bool Mesh::has_active_edge(const HalfEdge &halfedge) const {
  return _active_edges.find(halfedge.get_index()) != _active_edges.end();
}
bool Mesh::has_active_edge(const HalfEdgeIndex &index) const {
  return _active_edges.find(index) != _active_edges.end();
}
bool Mesh::has_checked_edge(const HalfEdge &halfedge) const {
  return _checked_edges.find(halfedge.get_index()) != _checked_edges.end();
}
bool Mesh::has_checked_edge(const HalfEdgeIndex &index) const {
  return _checked_edges.find(index) != _checked_edges.end();
}
bool Mesh::has_outgoing_next(const HalfEdge &working_edge,
                             const Point &point) const {
  for (HalfEdgeIndex outgoing :
       _mesh_points[working_edge.get_B()].get_outgoing()) {
    if (_mesh_edges[outgoing].get_edge().B() == point) return true;
  }
  return false;
}
bool Mesh::has_incoming_prev(const HalfEdge &working_edge,
                             const Point &point) const {
  for (HalfEdgeIndex outgoing :
       _mesh_points[working_edge.get_A()].get_outgoing()) {
    HalfEdgeIndex incoming = _mesh_edges[outgoing].get_previous();
    if (_mesh_edges[incoming].get_edge().A() == point) return true;
  }
  return false;
}

size_t Mesh::get_active_edges_size() const { return _active_edges.size(); }
size_t Mesh::get_checked_edges_size() const { return _checked_edges.size(); }
bool Mesh::active_edges_empty() const { return (get_active_edges_size() == 0); }
bool Mesh::checked_edges_empty() const {
  return (get_checked_edges_size() == 0);
}
bool Mesh::is_boundary_point(const MeshPoint &P) const {
  for (auto halfedge : P.get_outgoing()) {
    if (_mesh_edges[halfedge].is_boundary()) return true;
  }
  for (auto halfedge : P.get_outgoing()) {
    HalfEdgeIndex incoming = get_previous_index(halfedge);
    if (_mesh_edges[incoming].is_boundary()) return true;
  }
  return false;
}

// Adding triangle ABP where P is new_point
void Mesh::add_triangle(HalfEdgeIndex i_AB, Point new_point,
                        std::string type /*new, next, previous, overlap, fill*/,
                        MeshPointIndex i_P) {
  assertm(type == "new" || type == "next" || type == "previous" ||
              type == "overlap" || type == "fill",
          "Invalid type!");

  if (type != "new") assertm(i_P >= 0, "Invalid index of given meshpoint.");

  HalfEdgeIndex i_edge = _mesh_edges.size();
  MeshPointIndex i_point = _mesh_points.size();
  FaceIndex i_face = _mesh_triangles.size();

  if (type == "new") {
    MeshPoint new_meshpoint(new_point);
    new_meshpoint.set_index(i_point);
    _mesh_points.push_back(new_meshpoint);
    i_P = i_point;
  }
  // std::cout << i_AB << endl;
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
  BA.set_index(i_BA);
  _bound_opposite(&AB, i_AB, &BA, i_BA);

  assertm(AB.is_inside() && BA.is_inside(),
          "The edge is not set to inside edge!");

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

  AP.set_index(i_AP);
  PB.set_index(i_PB);

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

  _mesh_edges.push_back(BA);
  _mesh_edges.push_back(AP);
  _mesh_edges.push_back(PB);

  if (type != "new") {
    type = _find_type(i_AB, P, i_P);
  }

  if (type == "fill") {  // the edges are inside and we need to bound them with
                         // opposite edges
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
  } else if (type == "previous") {           // only edge PB is active edge
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(i_PB);
  } else if (type == "next") {               // only edge AP is active edge
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
    _mesh_edges[i_AP].set_active();
    _active_edges.insert(i_AP);
  } else if (type == "new" ||
             type == "overlap") {  // both AP and PB are new active edges
    _mesh_edges[i_AP].set_active();
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(i_PB);
    _active_edges.insert(i_AP);
  }
  /*
    Insert the new faces to lists
  */
  _mesh_triangles.push_back(F);
}

// returns vector of points inside Delaunay sphere
vector<MeshPoint> Mesh::get_breakers(const Triangle &triangle) const {
  assertm(triangle.is_triangle(), "Getting breakers of non valid triangle!");
  Point circumcenter = triangle.get_circumcenter();
  numeric dist = Vector(circumcenter, triangle.A()).get_length();

  vector<MeshPoint> breakers;
  for (MeshPoint vertex : _mesh_points) {
    numeric dist1 = Vector(circumcenter, vertex.get_point()).get_length();
    if (dist1 < 1.1 * dist && vertex.get_point() != triangle.A() &&
        vertex.get_point() != triangle.B() &&
        vertex.get_point() != triangle.C()) {
      bool found = false;

      if (is_boundary_point(vertex)) breakers.push_back(vertex);
    }
  }
  return breakers;
}

// checks Delaunay constraint for triangle T
bool Mesh::check_Delaunay(const Mesh &mesh, const Triangle &new_triangle,
                          const HalfEdge &working_edge,
                          const Face &incident_face) const {
  MeshPoint other_point =
      mesh.get_meshpoint(mesh.get_next_halfedge(working_edge).get_B());

  assertm(new_triangle.is_triangle(),
          "Checking Delaunay of non valid triangle!");

  Point circumcenter = new_triangle.get_circumcenter();

  numeric distA = Vector(circumcenter, new_triangle.A()).get_length();
  numeric distB = Vector(circumcenter, new_triangle.B()).get_length();
  numeric distC = Vector(circumcenter, new_triangle.C()).get_length();
  assertm(abs(distA - distB) + abs(distB - distC) + abs(distC - distA) < kEps,
          "Wrong circumcenter in Delaunay!");
  numeric dist = std::min(std::min(distA, distB), distC);
  const Triangle &T = new_triangle;

  // gravity center condition
  for (Face face : mesh.get_mesh_faces()) {
    Point gravity_center = face.get_gravity_center();
    numeric gc_dist = Vector(circumcenter, gravity_center).get_length();

    if (gc_dist < dist - kEps) {
      if (!(T.AB() % face.AB() || T.AB() % face.BC() || T.AB() % face.CA() ||
            T.BC() % face.AB() || T.BC() % face.BC() || T.BC() % face.CA() ||
            T.CA() % face.AB() || T.CA() % face.BC() || T.CA() % face.CA() ||
            T.C() == face.A() || T.C() == face.B() || T.C() == face.C())) {
        //std::cout << "Delaunay returned FALSE0!" << endl;
        return false;
      }
    }

    /*
    Point Tr_circumcenter = Tr.get_circumcenter();
    numeric tr_radius = Vector(Tr_circumcenter, Tr.A()).get_length();
    numeric tr_dist = Vector(Tr_circumcenter, T.C()).get_length();
    if(tr_dist<tr_radius)
    {
      //breaking others triangle delaunay
      return false;
    }
    */

    if (face.A() != T.A() && face.A() != T.B() && face.A() != T.C() &&
        face.A() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.A()).get_length();
      if (dist1 < dist - 10 * kEps) {
        //std::cout << "Delaunay returned FALSE1!" << endl;
        return false;
      }
    }
    if (face.B() != T.A() && face.B() != T.B() && face.B() != T.C() &&
        face.B() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.B()).get_length();
      if (dist1 < dist - 10 * kEps) {
        //std::cout << "Delaunay returned FALSE2!" << endl;
        return false;
      }
    }
    if (face.C() != T.A() && face.C() != T.B() && face.C() != T.C() &&
        face.C() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.C()).get_length();
      if (dist1 < dist - 10 * kEps) {
        //std::cout << "Delaunay returned FALSE3!" << endl;
        return false;
      }
    }
  }

  return true;
}

// makes .obj file from mesh triangles
void Mesh::obj_format(const std::string &name) const {
  std::ofstream out(name + ".obj");
  for (size_t i = 0; i < _mesh_points.size(); ++i) {
    out << "v " << _mesh_points[i].x() << " " << _mesh_points[i].y() << " "
        << _mesh_points[i].z() << endl;
  }
  for (size_t i = 0; i < _mesh_triangles.size(); ++i) {
    HalfEdge halfedge = get_halfedge(_mesh_triangles[i].get_halfedge());
    HalfEdge next_halfedge =
        get_next_halfedge(_mesh_triangles[i].get_halfedge());
    out << "f " << halfedge.get_A() + 1 << " " << halfedge.get_B() + 1 << " "
        << next_halfedge.get_B() + 1 << endl;
  }
}

// private member functions
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
  _active_edges.erase(i_edge1);
  _active_edges.erase(i_edge2);
}

void Mesh::_bound_face(const FaceIndex i_face, HalfEdge *edge1, HalfEdge *edge2,
                       HalfEdge *edge3) const {
  edge1->set_incident(i_face);
  edge2->set_incident(i_face);
  edge3->set_incident(i_face);
}

// having edge BA, we are trying to bound it with edge AB if it exists,
// setting AB inside edge
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

void Mesh::edges_check(const std::string &message,
                       const HalfEdgeIndex working_edge) const {
  std::cout << "Edges check " << message << endl;
  std::cout << _mesh_edges.size() << endl;
  for (int i = 0; i < _mesh_edges.size(); ++i) {
    const HalfEdge &edge = _mesh_edges[i];
    if (edge.get_index() == working_edge) {
      continue;
    }

    const auto edge_index = edge.get_index();
    bool is_ok = true;
    if (edge.is_checked() != is_checked(edge_index)) {
      std::cout << "Inconsistent checked edges!" << endl;
      is_ok = false;
    }
    if (edge.is_active() != is_active(edge_index)) {
      std::cout << "Inconsistent active edges!" << endl;
      is_ok = false;
    }
    if (edge.is_bounding()) {
      std::cout << "Found bounding edge!" << endl;
      is_ok = false;
    }
    if (edge.is_inside() !=
        (!is_active(edge_index) && !is_checked(edge_index))) {
      std::cout << "Inconsistent inside edges!" << endl;
      is_ok = false;
    }
    if (edge.is_boundary() && edge.get_opposite() != kInvalidEdgeIndex) {
      std::cout << "Inconsistent boundary edges case 1!" << endl;
      is_ok = false;
    }
    if (edge.is_inside() && edge.get_opposite() == kInvalidEdgeIndex) {
      std::cout << "Inconsistent boundary edges case 2!" << endl;
      is_ok = false;
    }
    if (edge.is_boundary() &&
        is_in_mesh(Edge(edge.get_point_B(), edge.get_point_A()))) {
      std::cout << "Inconsistent opposite edges!" << endl;
      is_ok = false;
    }
    if (edge != _mesh_edges[edge.get_index()]) {
      std::cout << "Inconsistent indices!" << endl;
      is_ok = false;
    }

    if (!is_ok) assertm(false, "Edge check failed!");
  }

  return;
}
std::string Mesh::_find_type(const HalfEdgeIndex index_AB, const MeshPoint &P,
                             MeshPointIndex index_P) const {
  const MeshPoint &A = _mesh_points[_mesh_edges[index_AB].get_A()];
  const MeshPoint &B = _mesh_points[_mesh_edges[index_AB].get_B()];

  bool is_previous = false, is_next = false;

  for (HalfEdgeIndex outgoing : B.get_outgoing()) {
    const HalfEdge &outgoing_edge = _mesh_edges[outgoing];
    const MeshPoint &outgoing_point = _mesh_points[outgoing_edge.get_B()];
    if (P == outgoing_point) {
      is_next = true;
      break;
    }
  }
  for (HalfEdgeIndex outgoing : A.get_outgoing()) {
    const HalfEdge &outgoing_edge = _mesh_edges[outgoing];
    const HalfEdge &incoming_edge = _mesh_edges[outgoing_edge.get_previous()];
    const MeshPoint &incoming_point = _mesh_points[incoming_edge.get_A()];
    if (P == incoming_point) {
      is_previous = true;
      break;
    }
  }
  if (is_previous && is_next) return "fill";
  if (is_previous) return "previous";
  if (is_next) return "next";
  return "overlap";
}