#pragma once

#include <ginac/ginac.h>

using HalfEdgeIndex = int;
using MeshPointIndex = int;
using FaceIndex = int;

static GiNaC::numeric kEps = GiNaC::numeric(10e-7);
void set_epsylon(GiNaC::numeric edge_size);

constexpr HalfEdgeIndex kInvalidEdgeIndex = -1;
constexpr MeshPointIndex kInvalidPointIndex = -1;
constexpr FaceIndex kInvalidFaceIndex = -1;