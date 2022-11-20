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
  MeshPoint p1(1, 2, 3);
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
  HalfEdge HExy(Exy, 1, 2);
  Face F0xy(-1, T1);
  Mesh M(F0xy);
  M.cout_triangles();
  EXPECT_TRUE(P0 == Point(0, 0, 0));
}
}  // namespace
