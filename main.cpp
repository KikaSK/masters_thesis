#include <ginac/ginac.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <queue>
#include <vector>

#include "algorithms.h"
#include "assertm.h"
#include "basic_algorithm.h"
#include "bounding_box.h"
#include "constants.h"
#include "edge.h"
#include "function.h"
#include "mesh.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"

using std::cout;
using std::endl;
using std::queue;
using std::vector;
using namespace GiNaC;

// substitutes to function and returns numeric
numeric substitute(const Function F, GiNaC::ex il) {
  return ex_to<numeric>(F.get_function().subs(il).evalf());
}

// runs given input, creates output file in "/output" folder
void run_input(const int i, const string folder, const string index) {
  string index_str = index;

  std::ifstream input_file("./inputs" + folder + "/input" + std::to_string(i),
                           std::ifstream::in);
  assertm(input_file.is_open(), "Failed opening the file!");
  vector<string> parsed_input;
  for (string str; getline(input_file, str);) {
    string s;
    bool writing = false;
    for (char c : str) {
      if (writing && c != '"') {
        s.push_back(c);
      }
      if (c == '"') writing = !writing;
    }
    parsed_input.push_back(s);
  }

  realsymbol x("x"), y("y"), z("z");

  symtab table;
  table["x"] = x;
  table["y"] = y;
  table["z"] = z;

  parser reader(table);

  string name = index_str + "_" + parsed_input[0] + "_" + parsed_input[1];
  cout << "Triangulation of " << parsed_input[0] << " with edge size "
       << parsed_input[1] << endl;
  ex input_F = reader(parsed_input[2]);
  numeric min_x = ex_to<numeric>(stod(parsed_input[3]));
  numeric max_x = ex_to<numeric>(stod(parsed_input[4]));
  numeric min_y = ex_to<numeric>(stod(parsed_input[5]));
  numeric max_y = ex_to<numeric>(stod(parsed_input[6]));
  numeric min_z = ex_to<numeric>(stod(parsed_input[7]));
  numeric max_z = ex_to<numeric>(stod(parsed_input[8]));
  BoundingBox my_bounding_box(min_x, max_x, min_y, max_y, min_z, max_z);
  numeric e_size = ex_to<numeric>(stod(parsed_input[9]));
  numeric seed_point_x = ex_to<numeric>(stod(parsed_input[10]));
  numeric seed_point_y = ex_to<numeric>(stod(parsed_input[11]));
  numeric seed_point_z = ex_to<numeric>(stod(parsed_input[12]));
  Point seed_point(seed_point_x, seed_point_y, seed_point_z);

  set_epsylon(e_size);
  vector<ex> input_dF;
  input_dF.push_back(diff(input_F, x));
  input_dF.push_back(diff(input_F, y));
  input_dF.push_back(diff(input_F, z));

  Function F(x, y, z, input_F, input_dF);

  BasicAlgorithm alg("./outputs/" + name, F, seed_point, e_size, x, y, z,
                     my_bounding_box);
  cout << "Basic algorithm created, calling for calculate()!" << endl;
  // alg.my_mesh.obj_format("./outputs/" + name);

  /*assertm(seed_triangle.A() == sing_point,
          "Triangle not created with singular point.");
  assertm(F.eval_at_point(sing_point) < kEps,
          "Singular point not lying on surface!");
  assertm(false, "breakpoint");*/
  alg.calculate();
}

// runs multiple inputs, creates output files in "/output" folder
void run_all(const int beg, const int end, const string folder,
             const string index) {
  for (int i = beg; i <= end; ++i) run_input(i, folder, index);
}

void test_find_seed_triangle();

int main() {
  Digits = 15;

  // spustame prikazom "make"

  // spusti vstup "input0" v priecinku "inputs/finite_surfaces/sphere",
  // vystupny subor vlozi do priecinka "/outputs" a nazve ho s predponou
  // "my_run_input"
  // run_input(0, "/finite_surfaces/sphere", "my_run_input");
  // run_all(3, 4, "/finite_surfaces/sphere", "my_run_input");
  run_input(0, "/sing_surfaces/A1", "my_run_input");
  run_input(0, "/sing_surfaces/A2", "my_run_input");
  run_input(0, "/sing_surfaces/A3", "my_run_input");
  run_input(0, "/sing_surfaces/A4", "my_run_input");
  run_input(0, "/sing_surfaces/A5", "my_run_input");
  // run_all(0, 3, "/finite_surfaces/sphere", "my_run_input");
  //   run_all(0, 1, "/finite_surfaces/blobby", "my_run_input");
  //   run_all(0, 1, "/finite_surfaces/torus", "my_run_input");
  //   run_all(3, 3, "/finite_surfaces/cubed_sphere", "my_run_input");
  //   run_input(0, "/cut_surfaces", "crop_to_box");
  //     run_all(3, 3, "/finite_surfaces/cubed_sphere", "my_run_input");
  //     run_all(0, 0, "/finite_surfaces/ellipsoid", "my_run_input");
  //     run_all(1, 3, "/finite_surfaces/tetrahedron", "my_run_input");
  //     run_all(0, 3, "/finite_surfaces/joined_spheres", "my_run_input");
  //     run_all(0, 3, "/finite_surfaces/genus", "my_run_input");

  // spusti vstupy "input0", "input1", "input2" v priecinku
  // "inputs/infinite_surfaces/hyperboloid", vystupny subor vlozi do priecinka
  // "/outputs" a nazve ho s predponou "my_run_all"
  // run_all(0, 2, "/infinite_surfaces/hyperboloid", "my_run_all");
}

// testing
/*
void test_find_seed_triangle() {
  realsymbol x("x"), y("y"), z("z");
  numeric e_size = 0.01;
  {
    BoundingBox test_bounding_box =
        BoundingBox(-1.5, 1.5, -1.5, 1.5, -1.5, 1.5);
    ex input_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 1;
    vector<ex> input_dF;
    input_dF.push_back(2 * x);
    input_dF.push_back(2 * y);
    input_dF.push_back(2 * z);

    Function F(x, y, z, input_F, input_dF);
    Point seed(1, 0, 0);

    Triangle t = find_seed_triangle(F, seed, e_size, test_bounding_box);

    auto total_length =
        t.AB().get_length() + t.BC().get_length() + t.CA().get_length();
    assertm(abs(total_length - e_size * 3) < 10e-2,
            "Triangle not of specified length.");

    numeric value1 =
        F.substitute(lst{x == t.A().x(), y == t.A().y(), z == t.A().z()});
    numeric value2 =
        F.substitute(lst{x == t.B().x(), y == t.B().y(), z == t.B().z()});
    numeric value3 =
        F.substitute(lst{x == t.C().x(), y == t.C().y(), z == t.C().z()});

    assertm(value1 < 10e-6 && value2 < 10e-6 && value3 < 10e-6, "Bad test!");
  }
  {
    BoundingBox test_bounding_box =
        BoundingBox(-1.5, 1.5, -1.5, 1.5, -1.5, 1.5);
    ex input_G =
        pow((x - 1) / 2, 2) + pow((y - 1) / 3, 2) + pow((z - 1), 2) - 1;
    vector<ex> input_dG;
    input_dG.push_back((x - 1) / 2);
    input_dG.push_back(2 * (y - 1) / 3);
    input_dG.push_back(2 * (z - 1));

    Function G(x, y, z, input_G, input_dG);
    Point seed(1, 1, 2);

    Triangle t = find_seed_triangle(G, seed, e_size, test_bounding_box);
    auto total_length =
        t.AB().get_length() + t.BC().get_length() + t.CA().get_length();
    assertm(abs(total_length - e_size * 3) < 10e-2,
            "Triangle not of specified length.");

    numeric value1 =
        G.substitute(lst{x == t.A().x(), y == t.A().y(), z == t.A().z()});
    numeric value2 =
        G.substitute(lst{x == t.B().x(), y == t.B().y(), z == t.B().z()});
    numeric value3 =
        G.substitute(lst{x == t.C().x(), y == t.C().y(), z == t.C().z()});

    assertm(value1 < 10e-6 && value2 < 10e-6 && value3 < 10e-6, "Bad test!");
  }
}
*/
