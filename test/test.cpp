#include <memory>

#include "edge.h"
#include "face.h"
#include "gtest/gtest.h"
#include "half_edge.h"
#include "mesh.h"
#include "mesh_point.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"

#define eps 10e-6

namespace {

// TEST example

TEST(DIVISION, Zero) { EXPECT_EQ(5 / 1, 5); }

// TEST Point

TEST(POINT, Equality_True) {
  Point p1(1, 2, 3);
  Point q1(1, 2, 3);
  EXPECT_TRUE(p1 == q1);
  Point p2(0, 0, 0);
  Point q2(0, 0, 0);
  EXPECT_TRUE(p2 == q2);
  Point p3(0, 0, 0);
  Point q3(10e-9, 10e-9, 10e-9);
  EXPECT_TRUE(p3 == q3);
  Point p4(-1.984769347, 2.98764978, -6.56316);
  Point q4(-1.984769349, 2.98764978, -6.56316);
  EXPECT_TRUE(p4 == q4);
  Point p5(5, -4, -385683746);
  Point q5(5.000000002, -4.000000001, -385683746.000000006);
  EXPECT_TRUE(p5 == q5);
  Point p6(5, -4, -385683746);
  Point q6(5.00000000, -4.00000000, -385683746.000000006);
  EXPECT_TRUE(p6 == q6);
}

TEST(POINT, Equality_False) {
  Point p1(0, 0, 0);
  Point q1(10e-6, 10e-6, 10e-6);
  EXPECT_FALSE(p1 == q1);
  Point p2(0, 0, 0);
  Point q2(10e-6, 10e-6, 0);
  EXPECT_FALSE(p2 == q2);
  Point p3(0, 0, 0);
  Point q3(numeric(10e-5), 0, 0);
  EXPECT_FALSE(p3 == q3);
  Point p4(0, 0, 0);
  Point q4(numeric(0.00001), 0, 0);
  EXPECT_FALSE(p4 == q4);
  Point p5(1, 2, 3);
  Point q5(-1, 2, 3);
  EXPECT_FALSE(p5 == q5);
  Point p6(0, 0, 0);
  Point q6(10e-6, 0, 0);
  EXPECT_FALSE(p6 == q6);
}

TEST(POINT, Inquality_False) {
  Point p1(1, 2, 3);
  Point q1(1, 2, 3);
  EXPECT_FALSE(p1 != q1);
  Point p2(0, 0, 0);
  Point q2(0, 0, 0);
  EXPECT_FALSE(p2 != q2);
  Point p3(0, 0, 0);
  Point q3(10e-9, 10e-9, 10e-9);
  EXPECT_FALSE(p3 != q3);
  Point p4(-1.984769347, 2.98764978, -6.56316);
  Point q4(-1.984769349, 2.98764978, -6.56316);
  EXPECT_FALSE(p4 != q4);
  Point p5(5, -4, -385683746);
  Point q5(5.000000002, -4.000000001, -385683746.000000006);
  EXPECT_FALSE(p5 != q5);
  Point p6(5, -4, -385683746);
  Point q6(5.00000000, -4.00000000, -385683746.000000006);
  EXPECT_FALSE(p6 != q6);
}

TEST(POINT, Inequality_True) {
  Point p1(0, 0, 0);
  Point q1(10e-6, 10e-6, 10e-6);
  EXPECT_TRUE(p1 != q1);
  Point p2(0, 0, 0);
  Point q2(10e-6, 10e-6, 0);
  EXPECT_TRUE(p2 != q2);
  Point p3(0, 0, 0);
  Point q3(numeric(10e-5), 0, 0);
  EXPECT_TRUE(p3 != q3);
  Point p4(0, 0, 0);
  Point q4(numeric(0.00001), 0, 0);
  EXPECT_TRUE(p4 != q4);
  Point p5(1, 2, 3);
  Point q5(-1, 2, 3);
  EXPECT_TRUE(p5 != q5);
  Point p6(0, 0, 0);
  Point q6(10e-6, 0, 0);
  EXPECT_TRUE(p6 != q6);
}

TEST(POINT, Get) {
  Point p1(1, 2, 3);
  EXPECT_TRUE(p1.x() == 1);
  EXPECT_TRUE(p1.y() == 2);
  EXPECT_TRUE(p1.z() == 3);
  Point p2(1.098954325076408372837608, 2, 3);
  EXPECT_NEAR(p2.x().to_double(), 1.09895432507640837846406, eps);
  Point p4(1.098954325076408372837608, 2, 3);
  EXPECT_TRUE(p4.x() == 1.09895432507640837846406);
  Point p3(-1.38536296, 2.28375802, 3.3285740386264);
  EXPECT_TRUE(p3.x() == -1.38536296);
  EXPECT_TRUE(p3.y() == 2.28375802);
  EXPECT_TRUE(p3.z() == 3.3285740386264);
  EXPECT_FALSE(p3.x() == -1.385362961);
  EXPECT_FALSE(p3.y() == 2.28375802376);
  EXPECT_FALSE(p3.z() == 3.328574038626);
}

// TEST Vector

TEST(VECTOR, Length) {
  Vector v1(1, 2, 3);
  EXPECT_TRUE(v1.get_length() == sqrt(numeric(14)));
  Vector v2(Point(1, 2, 3), Point(0, 0, 0));
  EXPECT_TRUE(v2.get_length() == sqrt(numeric(14)));
}

TEST(VECTOR, Unit) {
  Vector v1(1, 2, 3);
  EXPECT_TRUE(v1.unit().get_length() - 1 < eps);
  EXPECT_TRUE(v1.unit() == Vector(1 / sqrt(numeric(14)), 2 / sqrt(numeric(14)),
                                  3 / sqrt(numeric(14))));
}

TEST(VECTOR, Perpendicular) {
  Vector v1(1, 2, 3);
  EXPECT_TRUE(v1.get_any_perpendicular() * v1 < eps);
  Vector v2(10e-6, 10e-6, 10e-6);
  EXPECT_TRUE(v2.get_any_perpendicular() * v2 < eps);
  Vector v3(-1.984769347, 2.98764978, -6.56316);
  EXPECT_TRUE(v3.get_any_perpendicular() * v3 < eps);
  Vector v4(5, -4, -385683746);
  EXPECT_TRUE(v4.get_any_perpendicular() * v4 < eps);
  Vector v5(5, -4, -385683746);
  EXPECT_TRUE(v5.get_any_perpendicular() * v5 < eps);
}

TEST(VECTOR, Cross_Product) {
  Vector v1(1, 2, 3);
  Vector q1(8, 9, 3);
  EXPECT_TRUE((v1 ^ q1) == Vector(-21, 21, -7));
  Vector v2(10e-3, 10e-3, 10e-3);
  EXPECT_TRUE((v2 ^ q1) ==
              Vector(-numeric(3) / 50, numeric(1) / 20, numeric(1) / 100));
  Vector v4(5, -4, -385683746);
  EXPECT_TRUE((v4 ^ q1) == Vector(3471153702, -3085469983, 77));
}

// TEST Edge

TEST(EDGE, LengthMidpoint) {
  GiNaC::Digits = 15;
  Point A = Point(1, 2, 3);
  Point B = Point(0, 0, 0);
  Edge e1(A, B);
  EXPECT_TRUE(B.x() == 0);
  EXPECT_TRUE(e1.get_length() == sqrt(numeric(14)));
  EXPECT_TRUE(e1.get_midpoint() ==
              Point(numeric(1) / 2, numeric(1), numeric(3) / 2));
}

TEST(EDGE, Equality) {}

// TEST MeshPoint

TEST(MESHPOINT, Equality) {
  MeshPoint p1(numeric(1), numeric(2), numeric(3));
  MeshPoint q1(1, 2, 3);
  EXPECT_TRUE(p1 == q1);
  MeshPoint p3(0, 0, 0);
  MeshPoint q3(10e-9, 10e-9, 10e-9);
  EXPECT_TRUE(p3 == q3);
  MeshPoint p4(0, 0, 0);
  MeshPoint q4(numeric(0.00001), 0, 0);
  EXPECT_FALSE(p4 == q4);
  MeshPoint p5(1, 2, 3);
  MeshPoint q5(-1, 2, 3);
  EXPECT_FALSE(p5 == q5);
  MeshPoint p6(0, 0, 0);
  MeshPoint q6(10e-6, 0, 0);
  EXPECT_FALSE(p6 == q6);
  MeshPoint p7(1, 2, 3);
  MeshPoint q7(1, 2, 3);
  EXPECT_FALSE(p7 != q7);
  MeshPoint p8(0, 0, 0);
  MeshPoint q8(10e-6, 10e-6, 10e-6);
  EXPECT_TRUE(p8 != q8);
}

TEST(MESHPOINT, Constructors) {
  Point A(1, 2, 3);
  Vector u(1, 2, 3);
  MeshPoint p1(1, 2, 3);
  MeshPoint p2(A);
  MeshPoint p3(A, u);
  MeshPoint p4(p1);
  MeshPoint p6(p1, u);

  EXPECT_TRUE(p1 == MeshPoint(1, 2, 3));
  EXPECT_TRUE(p3 == MeshPoint(2, 4, 6));
  EXPECT_TRUE(p4 == p1);
  EXPECT_TRUE(p6 == p3);
}

TEST(MESHPOINT, Get) {
  MeshPoint p1(1, 2, 3);
  EXPECT_TRUE(p1.x() == 1);
  EXPECT_TRUE(p1.y() == 2);
  EXPECT_TRUE(p1.z() == 3);
  MeshPoint p2(1.098954325076408372837608, 2, 3);
  EXPECT_NEAR(p2.x().to_double(), 1.09895432507640837846406, eps);
  MeshPoint p4(1.098954325076408372837608, 2, 3);
  EXPECT_TRUE(p4.x() == 1.09895432507640837846406);
  MeshPoint p3(-1.38536296, 2.28375802, 3.3285740386264);
  EXPECT_TRUE(p3.x() == -1.38536296);
  EXPECT_TRUE(p3.y() == 2.28375802);
  EXPECT_TRUE(p3.z() == 3.3285740386264);
  EXPECT_FALSE(p3.x() == -1.385362961);
  EXPECT_FALSE(p3.y() == 2.28375802376);
  EXPECT_FALSE(p3.z() == 3.328574038626);
}

TEST(MESHPOINT, Addition) {
  MeshPoint p1(1, 2, 3);
  EXPECT_TRUE(p1 + p1 == Point(2, 4, 6));
}

// TEST Mesh

TEST(MESH, CreateMesh) {
  Point P0(0, 0, 0);
  Point Px(1, 0, 0);
  Point Py(0, 1, 0);
  Point Pz(0, 0, 1);
  Edge Exy(Px, Py);
  Triangle T1(P0, Px, Py);
  Face F0xy(T1);
  Mesh M(F0xy);

  // Test all previous, next, opposite relationships of P0Px
  EXPECT_TRUE((M.get_previous_halfedge(0)).get_edge() ==
              Edge(Py, P0));  // Previous of P0Px should be PyP0
  EXPECT_TRUE((M.get_next_halfedge(0)).get_edge() ==
              Edge(Px, Py));                   // Next of P0Px should be PxPy
  EXPECT_TRUE(M.get_opposite_index(0) == -1);  // P0Px does not have opposite

  // Test all previous, next, opposite relationships of PxPy
  EXPECT_TRUE((M.get_previous_halfedge(1)).get_edge() ==
              Edge(P0, Px));  // Previous of PxPy should be P0Px
  EXPECT_TRUE((M.get_next_halfedge(1)).get_edge() ==
              Edge(Py, P0));                   // Next of PxPy should be PyP0
  EXPECT_TRUE(M.get_opposite_index(1) == -1);  // PxPy does not have opposite

  // Test all previous, next, opposite relationships of PyP0
  EXPECT_TRUE((M.get_previous_halfedge(2)).get_edge() ==
              Edge(Px, Py));  // Previous of PyP0 should be PxPy
  EXPECT_TRUE((M.get_next_halfedge(2)).get_edge() ==
              Edge(P0, Px));                   // Next of PyP0 should be P0Px
  EXPECT_TRUE(M.get_opposite_index(2) == -1);  // PyP0 does not have opposite

  // Test outgoing edges
  EXPECT_TRUE(M.get_halfedge(M.get_meshpoint(0).get_outgoing()[0]).get_edge() ==
              Edge(P0, Px));  // Outgoing from P0 is P0Px
  EXPECT_TRUE(M.get_halfedge(M.get_meshpoint(1).get_outgoing()[0]).get_edge() ==
              Edge(Px, Py));  // Outgoing from Px is PxPy
  EXPECT_TRUE(M.get_halfedge(M.get_meshpoint(2).get_outgoing()[0]).get_edge() ==
              Edge(Py, P0));  // Outgoing from Py is PyP0

  // Test edge points
  EXPECT_TRUE(M.get_halfedge(0).get_point_A() == P0);
  EXPECT_TRUE(M.get_halfedge(0).get_point_B() == Px);
  EXPECT_TRUE(M.get_halfedge(1).get_point_A() == Px);
  EXPECT_TRUE(M.get_halfedge(1).get_point_B() == Py);
  EXPECT_TRUE(M.get_halfedge(2).get_point_A() == Py);
  EXPECT_TRUE(M.get_halfedge(2).get_point_B() == P0);

  // Test incident face
  EXPECT_TRUE(M.get_face(0).get_triangle() == T1);
  EXPECT_TRUE(M.get_halfedge(0).get_incident() == 0);
  EXPECT_TRUE(M.get_halfedge(1).get_incident() == 0);
  EXPECT_TRUE(M.get_halfedge(2).get_incident() == 0);
  EXPECT_TRUE(M.get_face(0).get_halfedge() == 0 ||
              M.get_face(0).get_halfedge() == 1 ||
              M.get_face(0).get_halfedge() == 2);

  // Test active edges list
  EXPECT_TRUE(M.has_active_edge(0));
  EXPECT_TRUE(M.has_active_edge(1));
  EXPECT_TRUE(M.has_active_edge(2));
}

TEST(MESH, AddNewTriangle) {
  Point P0(0, 0, 0);
  Point Px(1, 0, 0);
  Point Py(0, 1, 0);
  Point Pz(0, 0, 1);
  Edge Exy(Px, Py);
  Triangle T1(P0, Px, Py);
  Face F0xy(T1);
  Mesh M(F0xy);

  // Add triangle P0PxPz to edge P0Px

  HalfEdgeIndex i_P0Px(0), i_PxPy(1), i_PyP0(2);
  M.add_triangle(i_P0Px, Pz, "new");
  HalfEdgeIndex i_PxP0(3), i_P0Pz(4), i_PzPx(5);

  // All edges
  const HalfEdge& P0Px = M.get_halfedge(i_P0Px);
  const HalfEdge& PxPy = M.get_halfedge(i_PxPy);
  const HalfEdge& PyP0 = M.get_halfedge(i_PyP0);
  const HalfEdge& PxP0 = M.get_halfedge(i_PxP0);
  const HalfEdge& P0Pz = M.get_halfedge(i_P0Pz);
  const HalfEdge& PzPx = M.get_halfedge(i_PzPx);

  // All points
  const MeshPoint& mp_P0 = M.get_meshpoint(P0Px.get_A());
  const MeshPoint& mp_Px = M.get_meshpoint(P0Px.get_B());
  const MeshPoint& mp_Py = M.get_meshpoint(PyP0.get_A());
  const MeshPoint& mp_Pz = M.get_meshpoint(P0Pz.get_B());

  // All faces
  const Face& F0 = M.get_face(0);
  const Face& F1 = M.get_face(1);

  // Test if P0Px and PxP0 are correctly set to inside edges
  EXPECT_TRUE(P0Px.is_inside());  // P0Px is inside edge
  EXPECT_TRUE(PxP0.is_inside());  // PxP0 is inside edge

  // Test if all other edges are active edges
  EXPECT_TRUE(PxPy.is_active());
  EXPECT_TRUE(PyP0.is_active());
  EXPECT_TRUE(P0Pz.is_active());
  EXPECT_TRUE(PzPx.is_active());

  // Test all previous/next relationships
  EXPECT_TRUE(P0Px.get_opposite() == i_PxP0);
  EXPECT_TRUE(PxP0.get_opposite() == i_P0Px);

  EXPECT_TRUE(P0Px.get_previous() == i_PyP0);
  EXPECT_TRUE(P0Px.get_next() == i_PxPy);
  EXPECT_TRUE(PxPy.get_previous() == i_P0Px);
  EXPECT_TRUE(PxPy.get_next() == i_PyP0);
  EXPECT_TRUE(PyP0.get_previous() == i_PxPy);
  EXPECT_TRUE(PyP0.get_next() == i_P0Px);

  EXPECT_TRUE(PxP0.get_previous() == i_PzPx);
  EXPECT_TRUE(PxP0.get_next() == i_P0Pz);
  EXPECT_TRUE(P0Pz.get_previous() == i_PxP0);
  EXPECT_TRUE(P0Pz.get_next() == i_PzPx);
  EXPECT_TRUE(PzPx.get_previous() == i_P0Pz);
  EXPECT_TRUE(PzPx.get_next() == i_PxP0);

  // Test outgoing relationships
  EXPECT_TRUE(mp_P0.has_outgoing(i_P0Px));
  EXPECT_TRUE(mp_P0.has_outgoing(i_P0Pz));
  EXPECT_TRUE(mp_Px.has_outgoing(i_PxPy));
  EXPECT_TRUE(mp_Py.has_outgoing(i_PyP0));
  EXPECT_TRUE(mp_Pz.has_outgoing(i_PzPx));

  // Test incident face
  EXPECT_TRUE(P0Px.get_incident() == 0);
  EXPECT_TRUE(PxPy.get_incident() == 0);
  EXPECT_TRUE(PyP0.get_incident() == 0);
  EXPECT_TRUE(PxP0.get_incident() == 1);
  EXPECT_TRUE(P0Pz.get_incident() == 1);
  EXPECT_TRUE(PzPx.get_incident() == 1);

  // Test face halfedges
  EXPECT_TRUE(F0.get_halfedge() == i_P0Px || F0.get_halfedge() == i_PxPy ||
              F0.get_halfedge() == i_PyP0);
  EXPECT_TRUE(F1.get_halfedge() == i_P0Pz || F1.get_halfedge() == i_PzPx ||
              F1.get_halfedge() == i_PxP0);

  // Test active edges list

  EXPECT_FALSE(M.has_active_edge(0));  // 0x
  EXPECT_TRUE(M.has_active_edge(1));   // xy
  EXPECT_TRUE(M.has_active_edge(2));   // y0
  EXPECT_FALSE(M.has_active_edge(3));  // x0
  EXPECT_TRUE(M.has_active_edge(4));   // 0z
  EXPECT_TRUE(M.has_active_edge(5));   // zx

  M.cout_triangles();

  /*
    std::cout << *M._active_edges.begin() << endl
              << *std::next(M._active_edges.begin()) << endl
              << *std::next(M._active_edges.begin(), 2);
  */
}
}  // namespace
