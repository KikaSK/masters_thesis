#pragma once

#include <ginac/ginac.h>

using HalfEdgeIndex = int;
using MeshPointIndex = int;
using FaceIndex = int;

static const GiNaC::numeric kEps = GiNaC::numeric(10e-6);
constexpr HalfEdgeIndex kInvalidEdgeIndex = -1;
constexpr MeshPointIndex kInvalidPointIndex = -1;
constexpr FaceIndex kInvalidFaceIndex = -1;