#include <ginac/ginac.h>

#include <iostream>

#include "edge.h"
#include "face.h"
#include "half_edge.h"
#include "mesh.h"
#include "mesh_point.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"

int main() {
  GiNaC::Digits = 15;

  MeshPoint p(1, 2, 3);

  Point P0(0, 0, 0);
  Point Px(1, 0, 0);
  Point Py(0, 1, 0);
  Point Pz(0, 0, 1);
  Edge Exy(Px, Py);
  Triangle T1(P0, Px, Py);
  HalfEdge HExy(Exy, 1, 2);
  Mesh M(T1);
  M.cout_triangles();

  std::cout << *M._active_edges.begin() << endl
            << *std::next(M._active_edges.begin()) << endl
            << *std::next(M._active_edges.begin(), 2);

  // add triangle
  M.add_triangle(0, Pz, "new");

  M.cout_triangles();

  std::cout << "SUCCESS" << std::endl;
}
