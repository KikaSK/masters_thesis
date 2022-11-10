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