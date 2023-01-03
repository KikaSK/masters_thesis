#include "mesh.h"

Mesh::Mesh(const Triangle &T, const BoundingBox &bounding_box) {
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

  _point_tree.insert(A);
  _point_tree.insert(B);
  _point_tree.insert(C);

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
void Mesh::set_active_edges(const unordered_set<HalfEdgeIndex> &edges) {
  _active_edges.clear();
  _active_edges = edges;
  return;
}
unordered_set<HalfEdgeIndex> Mesh::get_checked_edges() const {
  return _checked_edges;
}
void Mesh::set_checked_edges(const unordered_set<HalfEdgeIndex> &edges) {
  _active_edges.clear();
  _active_edges = edges;
  return;
}
unordered_set<HalfEdgeIndex> Mesh::get_bounding_edges() const {
  return _bounding_edges;
}
void Mesh::move_active_edges_to_checked() {
  assertm(!checked_edges_empty(), "Moving empty checked edges!");
  _active_edges.clear();
  while (!checked_edges_empty()) {
    HalfEdgeIndex edge = get_checked_edge();
    remove_checked_edge(edge);
    add_edge_to_active(edge);
  }
  assertm(!active_edges_empty() && checked_edges_empty(),
          "Problem in moving checked edges!");
  return;
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

vector<MeshPoint> Mesh::get_prev(const HalfEdge &halfedge) const {
  vector<MeshPoint> prev;
  const MeshPoint &A = _mesh_points[halfedge.get_A()];
  for (HalfEdgeIndex outgoing_index : A.get_outgoing()) {
    HalfEdgeIndex previous_halfedge_index =
        _mesh_edges[outgoing_index].get_previous();
    if (is_boundary(previous_halfedge_index)) {
      const HalfEdge &previous_halfedge = _mesh_edges[previous_halfedge_index];
      prev.push_back(_mesh_points[previous_halfedge.get_A()]);
    }
  }
  return prev;
}

vector<MeshPoint> Mesh::get_next(const HalfEdge &halfedge) const {
  vector<MeshPoint> next;
  // const MeshPoint &A = _mesh_points[halfedge.get_A()];
  const MeshPoint &B = _mesh_points[halfedge.get_B()];
  for (HalfEdgeIndex next_halfedge_index : B.get_outgoing()) {
    if (is_boundary(next_halfedge_index)) {
      const HalfEdge &next_halfedge = _mesh_edges[next_halfedge_index];
      next.push_back(_mesh_points[next_halfedge.get_B()]);
    }
  }
  return next;
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
  // edges_check("add triangle beg", i_AB);
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
  Triangle ABP(A.get_point(), B.get_point(), P.get_point());
  Face F(ABP, i_AP);
  _bound_face(i_face, &BA, &AP, &PB);
  _mesh_edges.push_back(BA);
  _mesh_edges.push_back(AP);
  _mesh_edges.push_back(PB);

  if (type != "new") {
    type = _find_type(i_AB, P, i_P);
  }
  if (type == "fill") {
    // the edges are inside and we need to bound them with opposite edges
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
  } else if (type == "previous") {           // only edge PB is border edge
    _bound_opposite_outgoing(P, i_A, i_AP);  // try to bound AP with PA
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(i_PB);
  } else if (type == "next") {               // only edge AP is border edge
    _bound_opposite_outgoing(B, i_P, i_PB);  // try to bound PB with BP
    _mesh_edges[i_AP].set_active();
    _active_edges.insert(i_AP);
  } else if (type == "new" ||
             type == "overlap") {  // both AP and PB are new active edges
    _mesh_edges[i_AP].set_active();
    _mesh_edges[i_PB].set_active();
    _active_edges.insert(i_PB);
    _active_edges.insert(i_AP);
    _point_tree.insert(P);
  }
  _mesh_triangles.push_back(F);
}

// returns vector of points inside Delaunay sphere
vector<MeshPoint> Mesh::_linear_breakers_getter(
    const Triangle &triangle) const {
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

vector<MeshPoint> Mesh::_tree_breakers_getter(const Triangle &triangle) const {
  assertm(triangle.is_triangle(), "Getting breakers of non valid triangle!");
  Point circumcenter = triangle.get_circumcenter();
  numeric dist = Vector(circumcenter, triangle.A()).get_length();

  vector<MeshPoint> potential_breakers = get_meshpoints_in_interval(
      circumcenter.x() - 1.1 * dist, circumcenter.x() + 1.1 * dist,
      circumcenter.y() - 1.1 * dist, circumcenter.y() + 1.1 * dist,
      circumcenter.z() - 1.1 * dist, circumcenter.z() + 1.1 * dist);

  vector<MeshPoint> breakers;
  for (MeshPoint meshpoint : potential_breakers) {
    if (meshpoint.get_point() != triangle.A() &&
        meshpoint.get_point() != triangle.B() &&
        meshpoint.get_point() != triangle.B()) {
      breakers.push_back(meshpoint);
    }
  }
  return breakers;
}

vector<MeshPoint> Mesh::get_breakers(const Triangle &triangle) const {
  return _tree_breakers_getter(triangle);
}

vector<MeshPoint> Mesh::get_meshpoints_in_interval(numeric min_x, numeric max_x,
                                                   numeric min_y, numeric max_y,
                                                   numeric min_z,
                                                   numeric max_z) const {
  const MeshPoint min_point = MeshPoint(min_x, min_y, min_z);
  const MeshPoint max_point = MeshPoint(max_x, max_y, max_z);

  return _point_tree.between(min_point, max_point);
}

// checks Delaunay constraint for triangle T
bool Mesh::check_Delaunay(const HalfEdge &working_edge,
                          const Point &new_point) const {
  MeshPoint other_point =
      get_meshpoint(get_next_halfedge(working_edge).get_B());

  const Triangle new_triangle(working_edge.get_point_B(),
                              working_edge.get_point_A(), new_point);

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
  for (Face face : get_mesh_faces()) {
    Point gravity_center = face.get_gravity_center();
    numeric gc_dist = Vector(circumcenter, gravity_center).get_length();

    if (gc_dist < dist - 10 * kEps) {
      if (!(T.AB() % face.AB() || T.AB() % face.BC() || T.AB() % face.CA() ||
            T.BC() % face.AB() || T.BC() % face.BC() || T.BC() % face.CA() ||
            T.CA() % face.AB() || T.CA() % face.BC() || T.CA() % face.CA() ||
            T.C() == face.A() || T.C() == face.B() || T.C() == face.C())) {
        // std::cout << "Delaunay returned FALSE0!" << endl;
        return false;
      }
    }

    Point face_circumcenter = face.get_circumcenter();
    numeric face_radius =
        std::min(std::min(Vector(face_circumcenter, face.A()).get_length(),
                          Vector(face_circumcenter, face.B()).get_length()),
                 Vector(face_circumcenter, face.C()).get_length());
    // distance of new point from face circumcenter
    numeric face_dist = Vector(face_circumcenter, T.C()).get_length();
    if (face_dist < face_radius - 10 * kEps) {
      // breaking others triangle delaunay
      return false;
    }

    if (face.A() != T.A() && face.A() != T.B() && face.A() != T.C() &&
        face.A() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.A()).get_length();
      if (dist1 < dist - 10 * kEps) {
        // std::cout << "Delaunay returned FALSE1!" << endl;
        return false;
      }
    }
    if (face.B() != T.A() && face.B() != T.B() && face.B() != T.C() &&
        face.B() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.B()).get_length();
      if (dist1 < dist - 10 * kEps) {
        // std::cout << "Delaunay returned FALSE2!" << endl;
        return false;
      }
    }
    if (face.C() != T.A() && face.C() != T.B() && face.C() != T.C() &&
        face.C() != other_point.get_point()) {
      numeric dist1 = Vector(circumcenter, face.C()).get_length();
      if (dist1 < dist - 10 * kEps) {
        // std::cout << "Delaunay returned FALSE3!" << endl;
        return false;
      }
    }
  }

  return true;
}

// makes .obj file from mesh triangles
void Mesh::obj_format(const std::string &name) const {
  // std::cout << name << endl;
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
  _checked_edges.erase(i_edge1);
  _checked_edges.erase(i_edge2);
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
    assertm(outgoing > -1 && outgoing < _mesh_edges.size(),
            "Invalid halfedge index in _bound_opposite_outgoing function!");
    if (_mesh_edges[outgoing].get_B() /*second point of outgoing edge*/ ==
        i_B)  // we found edge AB
    {
      // std::cout << "index inside bound opposite outgoing: " << outgoing <<
      // endl; std::cout << "bounding " << outgoing << " with " << i_BA << endl;
      _bound_opposite(&(_mesh_edges[outgoing]), outgoing, &(_mesh_edges[i_BA]),
                      i_BA);
    }
  }
}

void Mesh::edges_check(const std::string &message,
                       const HalfEdgeIndex working_edge) const {
  std::cout << "Edges check " << message << endl;
  for (int i = 0; i < _mesh_edges.size(); ++i) {
    const HalfEdge &edge = _mesh_edges[i];
    if (edge.get_index() == working_edge) {
      continue;
    }

    const auto edge_index = edge.get_index();
    bool is_ok = true;
    // std::cout << "a" << endl;
    // std::cout << "Edges count: " << _mesh_edges.size() << endl;
    // std::cout << "Edge index: " << edge_index << endl;
    if (edge.is_checked() != is_checked(edge_index)) {
      std::cout << "Inconsistent checked edges!" << endl;
      is_ok = false;
    }
    if (edge_index == 545) std::cout << "ok?" << is_ok << endl;
    // std::cout << "b" << endl;
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

// angle BAP in range (-Pi, Pi) with respect to neighbour triangle
std::optional<numeric> Mesh::angle(const HalfEdge &working_edge, const Point &P,
                                   const Face &incident_face,
                                   const bool clockwise) const {
  if (P == working_edge.get_point_A() || P == working_edge.get_point_B())
    return std::nullopt;

  Point A = get_meshpoint(working_edge.get_A()).get_point();
  Point B = get_meshpoint(working_edge.get_B()).get_point();

  // change the points for other direction
  if (clockwise) std::swap(A, B);

  Triangle T = Triangle(A, B, P);
  if (!T.is_triangle()) return std::nullopt;
  assertm(incident_face.is_triangle(), "Invalid incident triangle!");

  // third point in neighbour triangle
  Point C = get_meshpoint(get_next_halfedge(working_edge).get_B()).get_point();

  Vector AB = Vector(A, B);
  Vector AP = Vector(A, P);
  Vector AC = Vector(A, C);

  assertm(!AB.is_zero() && !AP.is_zero() && !AC.is_zero(), "Zero vectors!");

  numeric cos = (AB.unit() * AP.unit());
  assertm(abs(cos) <= numeric(1), "Wrong value of cosine!");

  numeric angle = acos(cos);

  assertm(angle >= 0, "Angle in the wrong range!");

  Vector ABcrossAC = (AB ^ AC);
  Vector ABcrossAP = (AB ^ AP);

  // if same_direction is > 0 the normals are pointing in aprox same direction
  // thus the points lie on the same side

  // TODO(kuska) check if we need this if I have good normals

  numeric same_direction = ABcrossAP * ABcrossAC;

  if (same_direction > 0) {
    angle = -angle;
  }
  assertm(abs(angle) <= GiNaC::ex_to<numeric>(GiNaC::Pi.evalf()),
          "Angle in the wrong range!");

  return angle;
}

// true if angle is between 0 and 3*pi/4 with respect to neighbour triangle
bool Mesh::good_orientation(const HalfEdge &working_edge, const Point &P,
                            const Face &incident_face) const {
  std::optional<numeric> angle1 = angle(working_edge, P, incident_face, false);
  std::optional<numeric> angle2 = angle(working_edge, P, incident_face, false);
  assertm(angle1.has_value() && angle2.has_value(),
          "Angle does not have value!");
  return angle1.value() > 0 &&
         angle2.value() < 3 * GiNaC::ex_to<numeric>(GiNaC::Pi.evalf()) / 4;
}

// https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d

// returns ditance between point and line segment given by working edge
numeric Mesh::line_point_dist(const HalfEdge &working_edge, const Point &P,
                              const Face &incident_face) const {
  assertm(
      !Vector(working_edge.get_point_A(), working_edge.get_point_B()).is_zero(),
      "Edge is zero vector!");
  Vector AB_unit =
      Vector(working_edge.get_point_A(), working_edge.get_point_B()).unit();
  Vector AP = Vector(working_edge.get_point_A(), P);
  numeric t = AP * AB_unit;
  Point p = Point(working_edge.get_point_A(), t * AB_unit);

  std::optional<numeric> angle1 = angle(working_edge, P, incident_face, true);
  std::optional<numeric> angle2 =
      angle(working_edge, P, incident_face, false);  // anticlockwise

  assertm(angle1.has_value() && angle2.has_value(),
          "Angles does not have value!");

  if (abs(angle1.value()) <= GiNaC::ex_to<numeric>(GiNaC::Pi.evalf()) / 2 &&
      abs(angle2.value()) <= GiNaC::ex_to<numeric>(GiNaC::Pi.evalf()) / 2) {
    return Vector(P, p).get_length();
  } else if (abs(angle1.value()) >
             GiNaC::ex_to<numeric>(GiNaC::Pi.evalf()) / 2) {
    return Vector(working_edge.get_point_A(), P).get_length();
  } else {
    return Vector(working_edge.get_point_B(), P).get_length();
  }
  assertm(false, "Wrong line point distance!");
  return 1000;
}
