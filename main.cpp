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
  numeric num_of_singular = ex_to<numeric>(stod(parsed_input[13]));
  vector<Point> singular_points;
  vector<vector<Vector>> singular_directions;
  vector<int> types;  // not used
  for (int i = 0; i < num_of_singular - 0.5; ++i) {
    numeric p_x = ex_to<numeric>(stod(parsed_input[14 + 3 * i]));
    numeric p_y = ex_to<numeric>(stod(parsed_input[15 + 3 * i]));
    numeric p_z = ex_to<numeric>(stod(parsed_input[16 + 3 * i]));
    singular_points.push_back(Point(p_x, p_y, p_z));
    singular_directions.push_back(vector<Vector>());
    numeric num_of_branches = ex_to<numeric>(stod(parsed_input[17 + 3 * i]));
    for (int j = 0; j < num_of_branches - 0.5; ++j) {
      numeric d_x = ex_to<numeric>(stod(parsed_input[18 + 3 * i + 3 * j]));
      numeric d_y = ex_to<numeric>(stod(parsed_input[19 + 3 * i + 3 * j]));
      numeric d_z = ex_to<numeric>(stod(parsed_input[20 + 3 * i + 3 * j]));
      singular_directions[singular_directions.size() - 1].push_back(
          Vector(d_x, d_y, d_z));
    }
  }
  /*
  for (auto p : singular_points) std::cout << "point: " << p << endl;
  for (auto b : singular_directions) {
    for (auto d : b) std::cout << "direction: " << d << endl;
    std::cout << endl;
  }
  */
  set_epsylon(e_size);
  vector<ex> input_dF;
  input_dF.push_back(diff(input_F, x));
  input_dF.push_back(diff(input_F, y));
  input_dF.push_back(diff(input_F, z));

  Function F(x, y, z, input_F, input_dF);
  std::cout << F.get_function() << endl;

  BasicAlgorithm alg("./outputs/" + name, F, seed_point, e_size, x, y, z,
                     my_bounding_box, singular_points, singular_directions);
  cout << "Basic algorithm created, calling for calculate()!" << endl;
  // alg.my_mesh.obj_format("./outputs/" + name);

  /*assertm(seed_triangle.A() == sing_point,
          "Triangle not created with singular point.");
  assertm(F.eval_at_point(sing_point) < kEps,
          "Singular point not lying on surface!");
  assertm(false, "breakpoint");*/
  alg.calculate();
}
// runs given input, creates output file in "/output" folder
void run_input_plane(const int i, const string folder, const string index) {
  string index_str = index;

  std::ifstream input_file(
      "./inputs_plane" + folder + "/input" + std::to_string(i),
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
    std::cout << s << endl;
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
  numeric e_size = ex_to<numeric>(stod(parsed_input[1]));
  // ex input_F = reader(parsed_input[2]);
  numeric min_x = ex_to<numeric>(stod(parsed_input[2]));
  numeric max_x = ex_to<numeric>(stod(parsed_input[3]));
  numeric min_y = ex_to<numeric>(stod(parsed_input[4]));
  numeric max_y = ex_to<numeric>(stod(parsed_input[5]));
  numeric min_z = ex_to<numeric>(stod(parsed_input[6]));
  numeric max_z = ex_to<numeric>(stod(parsed_input[7]));
  BoundingBox my_bounding_box(min_x, max_x, min_y, max_y, min_z, max_z);
  numeric num_of_singular = ex_to<numeric>(stod(parsed_input[8]));
  vector<Point> singular_points;
  vector<ex> singular_equations;
  vector<numeric> heights;
  vector<vector<Vector>> singular_directions;
  vector<int> types;
  for (int i = 0; i < num_of_singular - 0.5; ++i) {
    string str_n = parsed_input[9 + 4 * i];
    int n = stoi(str_n);
    types.push_back(n);
    double p_y = stod(parsed_input[10 + 4 * i]);
    double p_z = stod(parsed_input[11 + 4 * i]);
    double h = stod(parsed_input[12 + 4 * i]);
    double k = 3.1415 / (2 * sqrt(std::pow(h, n + 1)));
    double q = 4 * h / (3.1415 * (n + 1));
    double p_x = -h - q;
    singular_points.push_back(Point(numeric(p_x), numeric(p_y), numeric(p_z)));
    // equation of A_n singularity translated to point p_x, p_y, p_z
    ex sing_equation = pow(x - p_x, n + 1) - pow(y - p_y, 2) - pow(z - p_z, 2);
    ex plane_q = -x - q;
    ex inter_sing_plane =
        sing_equation + plane_q - sqrt(pow(sing_equation, 2) + pow(plane_q, 2));
    ex cosine = x + q + q * cos(k * sqrt(pow(y - p_y, 2) + pow(z - p_z, 2)));
    ex cylinder = 4 * pow(h, n + 1) - pow(z - p_z, 2) - pow(y - p_y, 2);
    ex inter_cylinder_cosine =
        cosine + cylinder - sqrt(pow(cosine, 2) + pow(cylinder, 2));
    ex plane_0 = x;
    ex union_plane_0_cosine =
        inter_cylinder_cosine + plane_0 +
        sqrt(pow(inter_cylinder_cosine, 2) + pow(plane_0, 2));

    ex inter_plane_q_cosine =
        -plane_q + union_plane_0_cosine -
        sqrt(pow(plane_q, 2) + pow(union_plane_0_cosine, 2));

    ex union_res_equation =
        inter_plane_q_cosine + inter_sing_plane +
        sqrt(pow(inter_plane_q_cosine, 2) + pow(inter_sing_plane, 2));

    singular_equations.push_back(union_res_equation);
    std::vector<Vector> sing_dir = {Vector(1, 0, 0)};
    singular_directions.push_back(sing_dir);
  }
  ex final_eq = singular_equations[0];
  for (int i = 1; i < singular_equations.size(); ++i) {
    final_eq = final_eq + singular_equations[i] +
               sqrt(pow(final_eq, 2) + pow(singular_equations[i], 2));
  }
  std::cout << final_eq << endl;
  set_epsylon(e_size);

  vector<ex> input_dF;
  input_dF.push_back(diff(final_eq, x));
  input_dF.push_back(diff(final_eq, y));
  input_dF.push_back(diff(final_eq, z));

  Function F(x, y, z, final_eq, input_dF);
  BasicAlgorithm alg("./outputs/" + name, F, singular_points[0], e_size, x, y,
                     z, my_bounding_box, singular_points, singular_directions,
                     types);
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
void run_all_plane(const int beg, const int end, const string folder,
                   const string index) {
  for (int i = beg; i <= end; ++i) run_input_plane(i, folder, index);
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
  // run_input(0, "/sing_surfaces/A1", "my_run_input");
  // run_input(0, "/sing_surfaces/A2", "my_run_input");
  // run_input(0, "/sing_surfaces/A3", "my_run_input");
  // run_input(0, "/sing_surfaces/A4", "my_run_input");
  // run_input(0, "/sing_surfaces/A1", "my_run_input");
  run_all_plane(2, 2, "/1_singularity/A2", "test");
  // run_all(2, 2, "/sing_surfaces/A2-plane",
  //"test-plane-basic-height-case3-adaptive");
  //    run_all(0, 1, "/finite_surfaces/blobby", "my_run_input");
  //    run_all(0, 1, "/finite_surfaces/torus", "my_run_input");
  //    run_all(3, 3, "/finite_surfaces/cubed_sphere", "my_run_input");
  //    run_input(0, "/cut_surfaces", "crop_to_box");
  //      run_all(3, 3, "/finite_surfaces/cubed_sphere", "my_run_input");
  //      run_all(0, 0, "/finite_surfaces/ellipsoid", "my_run_input");
  //      run_all(1, 3, "/finite_surfaces/tetrahedron", "my_run_input");
  //      run_all(0, 3, "/finite_surfaces/joined_spheres", "my_run_input");
  //      run_all(0, 3, "/finite_surfaces/genus", "my_run_input");

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
