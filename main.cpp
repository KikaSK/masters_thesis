#include <ginac/ginac.h>

#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <queue>
#include <regex>
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
#include "singularity.h"
#include "triangle.h"
#include "vector.h"

using std::cout;
using std::endl;
using std::queue;
using std::vector;
using namespace GiNaC;
using namespace std::chrono;

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
  // regex for An-- singularities
  std::regex Amm_regex("A[0-9]+--.*");
  std::regex Apm_regex("A[0-9]+[+]-.*");
  std::regex Dmm_regex("D[0-9]+--.*");
  std::regex Dpm_regex("D[0-9]+[+]-.*");
  std::regex Epm_regex("E[0-9]+[+]-.*");
  std::regex Epp_regex("E[0-9]+[+][+].*");

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
  vector<int> types;
  vector<Singularity> singularities;
  for (int i = 0; i < num_of_singular - 0.5; ++i) {
    Singularity sing;
    if (std::regex_match(parsed_input[0], Amm_regex)) {
      sing.set_type(SingularityType::Amm);
    } else if (std::regex_match(parsed_input[0], Apm_regex)) {
      sing.set_type(SingularityType::Apm);
    } else if (std::regex_match(parsed_input[0], Dmm_regex)) {
      sing.set_type(SingularityType::Dmm);
    } else if (std::regex_match(parsed_input[0], Dpm_regex)) {
      sing.set_type(SingularityType::Dpm);
    } else if (std::regex_match(parsed_input[0], Epm_regex)) {
      sing.set_type(SingularityType::Epm);
    } else if (std::regex_match(parsed_input[0], Epp_regex)) {
      sing.set_type(SingularityType::Epp);
    } else {
      assertm(false, "Wrong input!");
    }
    string str_n;
    for (int i = 1; i < 4; ++i) {
      if (parsed_input[0][i] != '-' && parsed_input[0][i] != '+')
        str_n.push_back(parsed_input[0][i]);
      else
        break;
    }
    types.push_back(stoi(str_n));
    sing.set_n(stoi(str_n));
    numeric p_x = ex_to<numeric>(stod(parsed_input[14 + 3 * i]));
    numeric p_y = ex_to<numeric>(stod(parsed_input[15 + 3 * i]));
    numeric p_z = ex_to<numeric>(stod(parsed_input[16 + 3 * i]));
    singular_points.push_back(Point(p_x, p_y, p_z));
    sing.set_location(Point(p_x, p_y, p_z));
    singular_directions.push_back(vector<Vector>());
    numeric num_of_branches = ex_to<numeric>(stod(parsed_input[17 + 3 * i]));
    for (int j = 0; j < num_of_branches - 0.5; ++j) {
      numeric d_x = ex_to<numeric>(stod(parsed_input[18 + 3 * i + 3 * j]));
      numeric d_y = ex_to<numeric>(stod(parsed_input[19 + 3 * i + 3 * j]));
      numeric d_z = ex_to<numeric>(stod(parsed_input[20 + 3 * i + 3 * j]));
      singular_directions[singular_directions.size() - 1].push_back(
          Vector(d_x, d_y, d_z));
      sing.add_direction(Vector(d_x, d_y, d_z));
    }
    singularities.push_back(sing);
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

  // regular algorithm
  if (num_of_singular == 0) {
    Triangle seed_triangle =
        find_seed_triangle(F, seed_point, e_size, my_bounding_box);

    assertm(seed_triangle.AB() != seed_triangle.BC() &&
                seed_triangle.AB() != seed_triangle.CA() &&
                seed_triangle.BC() != seed_triangle.CA(),
            "Seed triangle contains duplicit edges!");

    BasicAlgorithm alg("./outputs/" + name, F, seed_triangle, e_size, x, y, z,
                       my_bounding_box);
    cout << "Basic algorithm created, calling for calculate()!" << endl;
    // alg.my_mesh.obj_format("./outputs/" + name);
    auto start = high_resolution_clock::now();
    alg.calculate();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Edge size: " << e_size << endl;
    std::cout << "Duration: " << duration.count() << endl;

  }
  // algorithm starting in singular points
  else {
    BasicAlgorithm alg("./outputs/" + name, F, seed_point, e_size, x, y, z,
                       my_bounding_box, singularities);
    cout << "Basic algorithm created, calling for calculate()!" << endl;
    // alg.my_mesh.obj_format("./outputs/" + name);

    alg.calculate();
  }
  return;
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
  vector<Singularity> singularities;
  for (int i = 0; i < num_of_singular - 0.5; ++i) {
    Singularity sing;
    string str_n = parsed_input[9 + 4 * i];
    int n = stoi(str_n);
    types.push_back(n);
    sing.set_n(n);
    sing.set_type(SingularityType::Amm);
    double p_y = stod(parsed_input[10 + 4 * i]);
    double p_z = stod(parsed_input[11 + 4 * i]);
    double h = stod(parsed_input[12 + 4 * i]);
    double k = 3.1415 / (2 * sqrt(std::pow(h, n + 1)));
    double q = 4 * h / (3.1415 * (n + 1));
    double p_x = -h - q;
    sing.set_location(Point(numeric(p_x), numeric(p_y), numeric(p_z)));
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
    sing.set_directions(sing_dir);
    singular_directions.push_back(sing_dir);
    singularities.push_back(sing);
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

  Point seed_point = Point(0, 0, 0);
  BasicAlgorithm alg("./outputs/plane/" + name, F, seed_point, e_size, x, y, z,
                     my_bounding_box, singularities);
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

vector<Point> get_polyline(const Function &F, const Function &G,
                           const Point &starting_point, const numeric &e_size,
                           const BoundingBox &bounding_box) {
  std::optional<Point> proj_starting_point = starting_point;
  assertm(proj_starting_point.has_value(), "No value!");
  int iter = 0;
  while ((!F.is_on(proj_starting_point.value()) ||
          !G.is_on(proj_starting_point.value())) &&
         iter < 30) {
    const Vector grad_F = F.get_gradient_at_point(proj_starting_point.value());
    const Vector grad_G = G.get_gradient_at_point(proj_starting_point.value());
    assertm(!grad_F.is_zero() && !grad_G.is_zero(), "Zero gradient!");
    proj_starting_point =
        project(proj_starting_point.value(), grad_F.unit(), F, e_size);
    assertm(proj_starting_point.has_value(), "No value!");
    proj_starting_point =
        project(proj_starting_point.value(), grad_G.unit(), G, e_size);
    assertm(proj_starting_point.has_value(), "No value!");
    iter++;
  }
  assertm(iter != 30, "Too many iterations!");

  // first branch of points

  vector<Point> points_1;
  Point last_point = proj_starting_point.value();
  points_1.push_back(last_point);
  Vector grad_F = F.get_gradient_at_point(last_point);
  Vector grad_G = G.get_gradient_at_point(last_point);
  Vector tangent_vector = (grad_F ^ grad_G).unit();

  int iter_2 = 0;
  while (bounding_box.is_inside(last_point) &&
         !bounding_box.is_on(last_point) && iter_2 < 10000) {
    std::optional<Point> new_point =
        Point(last_point.x() + e_size * tangent_vector.x(),
              last_point.y() + e_size * tangent_vector.y(),
              last_point.z() + e_size * tangent_vector.z());
    assertm(new_point.has_value(), "No value!");
    iter = 0;
    while ((!F.is_on(new_point.value()) || !G.is_on(new_point.value())) &&
           iter < 30) {
      const Vector grad_F = F.get_gradient_at_point(new_point.value());
      const Vector grad_G = G.get_gradient_at_point(new_point.value());
      assertm(!grad_F.is_zero() && !grad_G.is_zero(), "Zero gradient!");
      new_point = project(new_point.value(), grad_F.unit(), F, e_size);
      new_point = project(new_point.value(), grad_G.unit(), G, e_size);
      assertm(new_point.has_value(), "No value!");
      iter++;
    }
    assertm(iter != 30, "Too many iterations!");

    new_point =
        bounding_box.crop_to_box(last_point, new_point.value(), e_size, F);
    // std::cout << iter_2 << " new_point: " << new_point.value() << endl;

    // special case - closed curve
    if (iter_2 > 3 &&
        Vector(new_point.value(), proj_starting_point.value()).get_length() <
            3.0 * e_size / 4.0) {
      new_point = proj_starting_point;
      points_1.push_back(new_point.value());
      return points_1;
    }

    last_point = new_point.value();
    grad_F = F.get_gradient_at_point(last_point);
    grad_G = G.get_gradient_at_point(last_point);
    tangent_vector = (grad_F ^ grad_G).unit();
    points_1.push_back(last_point);
    iter_2++;
  }

  // second branch of points

  vector<Point> points_2;
  last_point = proj_starting_point.value();
  points_2.push_back(last_point);
  grad_F = F.get_gradient_at_point(last_point);
  grad_G = G.get_gradient_at_point(last_point);
  tangent_vector = (grad_F ^ grad_G).unit();

  iter_2 = 0;
  while (bounding_box.is_inside(last_point) &&
         !bounding_box.is_on(last_point) && iter_2 < 100) {
    std::optional<Point> new_point =
        Point(last_point.x() - e_size * tangent_vector.x(),
              last_point.y() - e_size * tangent_vector.y(),
              last_point.z() - e_size * tangent_vector.z());

    assertm(new_point.has_value(), "No value!");
    iter = 0;
    while ((!F.is_on(new_point.value()) || !G.is_on(new_point.value())) &&
           iter < 30) {
      const Vector grad_F = F.get_gradient_at_point(new_point.value());
      const Vector grad_G = G.get_gradient_at_point(new_point.value());
      assertm(!grad_F.is_zero() && !grad_G.is_zero(), "Zero gradient!");
      new_point = project(new_point.value(), grad_F.unit(), F, e_size);
      new_point = project(new_point.value(), grad_G.unit(), G, e_size);
      assertm(new_point.has_value(), "No value!");
      iter++;
    }
    assertm(iter != 30, "Too many iterations!");

    new_point =
        bounding_box.crop_to_box(last_point, new_point.value(), e_size, F);
    last_point = new_point.value();
    grad_F = F.get_gradient_at_point(last_point);
    grad_G = G.get_gradient_at_point(last_point);
    tangent_vector = (grad_F ^ grad_G).unit();
    points_2.push_back(last_point);
    iter_2++;
  }
  if (points_1.empty()) return points_2;
  if (points_2.empty()) return points_1;

  reverse(points_1.begin(), points_1.end());
  points_1.pop_back();
  points_1.insert(points_1.end(), points_2.begin(), points_2.end());
  std::cout << "finished polyline" << endl;
  return points_1;
}

void run_polyline(const string function_F, const string function_G,
                  const vector<Point> &seeds, const numeric &e_size,
                  const string type, const string name = "default_name") {
  realsymbol x("x"), y("y"), z("z");

  symtab table;
  table["x"] = x;
  table["y"] = y;
  table["z"] = z;

  parser reader(table);
  ex ex_F = reader(function_F);
  // ex ex_F = pow(x, 2) + pow(y, 2) + pow(z, 2) - 4;
  vector<ex> input_dF;
  input_dF.push_back(diff(ex_F, x));
  input_dF.push_back(diff(ex_F, y));
  input_dF.push_back(diff(ex_F, z));

  Function F(x, y, z, ex_F, input_dF);
  ex ex_G = reader(function_G);
  // ex ex_G = pow(x, 2) - pow(y, 2) + pow(z, 2) - 1;
  vector<ex> input_dG;
  input_dG.push_back(diff(ex_G, x));
  input_dG.push_back(diff(ex_G, y));
  input_dG.push_back(diff(ex_G, z));

  Function G(x, y, z, ex_G, input_dG);

  numeric min_x = -5;
  numeric max_x = 5;
  numeric min_y = -3.5;
  numeric max_y = 3.5;
  numeric min_z = -4;
  numeric max_z = 4;
  BoundingBox my_bounding_box(min_x, max_x, min_y, max_y, min_z, max_z);

  vector<vector<Point>> polylines;
  for (Point seed : seeds) {
    // const Point starting_point(sqrt(5.0 / 2.0), sqrt(3.0 / 2.0), 0);
    vector<Point> polyline = get_polyline(F, G, seed, e_size, my_bounding_box);
    for (auto P : polyline) {
      std::cout << P << endl;
    }
    std::cout << endl;
    polylines.push_back(polyline);
  }
  ex my_function;
  if (type == "intersection") {
    my_function = F.get_function() + G.get_function() +
                  sqrt(pow(F.get_function(), 2) + pow(G.get_function(), 2));
  } else if (type == "union") {
    my_function = F.get_function() + G.get_function() -
                  sqrt(pow(F.get_function(), 2) + pow(G.get_function(), 2));
  } else if (type == "difference") {
    my_function = F.get_function() - G.get_function() +
                  sqrt(pow(F.get_function(), 2) + pow(G.get_function(), 2));
  }
  vector<ex> d_my_function;
  d_my_function.push_back(diff(my_function, x));
  d_my_function.push_back(diff(my_function, y));
  d_my_function.push_back(diff(my_function, z));

  Function my_function_FG(x, y, z, my_function, d_my_function);

  BasicAlgorithm alg("./outputs/curves/" + name, my_function_FG, F, G,
                     polylines, e_size, x, y, z, my_bounding_box, type);
  alg.calculate();
  return;
}

int main() {
  Digits = 15;

  // spustame prikazom "make"

  // spusti vstup "input0" v priecinku "inputs/finite_surfaces/sphere",
  // vystupny subor vlozi do priecinka "/outputs" a nazve ho s predponou
  // "my_run_input"
  // run_all(4, 4, "/finite_surfaces/genus", "my_run_input");
  assertm(false, "TADAAAAA");
  run_all(8, 15, "/finite_surfaces/sphere_timer", "timer");
  // run_all_plane(3, 3, "/1_singularity/A2", "test");
  /*
    run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
                 0.45, "difference", "sphere_sphere_difference_0.45");
    run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
                 0.30, "difference", "sphere_sphere_difference_0.30");
    run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
                 0.15, "difference", "sphere_sphere_difference_0.15");
  */
  /*
   run_polyline("x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8", "x^2+y^2+z^2-5.5",
                {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
                 Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
                 Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
                0.15, "difference", "small_sphere_tanglecube_difference_0.15");
                */
  /*
run_polyline("x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8", "x^2+y^2+z^2-5.5",
  {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
   Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
   Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
  0.10, "difference", "small_sphere_tanglecube_difference_0.10");
run_polyline("x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8", "x^2+y^2+z^2-5.5",
  {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
   Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
   Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
  0.05, "difference", "small_sphere_tanglecube_difference_0.05");*/
  /*
  run_polyline("x^2+y^2+z^2-5.5", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
               {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
                Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
                Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
               0.30, "difference", "tanglecube_small_sphere_difference_0.30");
  run_polyline("x^2+y^2+z^2-5.5", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
               {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
                Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
                Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
               0.20, "difference", "tanglecube_small_sphere_difference_0.20");
  run_polyline("x^2+y^2+z^2-5.5", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
               {Point(-1.28903, 0.5, 1.89431), Point(-1.28903, 0.5, -1.89431),
                Point(1.28903, -1.89431, 0.5), Point(-1.89431, 1.28903, 0.5),
                Point(1.89431, -1.28903, 0.5), Point(-1.28903, 1.89431, 0.5)},
               0.10, "difference", "tanglecube_small_sphere_difference_0.10");*/
  /*
      run_polyline(
          "sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)"
          "+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1",
          "y", {Point(1.2038, 0, 0)}, 0.15, "blobby_plane_0.15");
      run_polyline(
          "sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)"
          "+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1",
          "y", {Point(1.2038, 0, 0)}, 0.05, "blobby_plane_0.05");
          */
  /*
run_polyline(
   "sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)"
   "+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1",
   "z", {Point(1.2038, 0, 0)}, 0.15, "blobby_plane_z_0.15");
run_polyline(
   "sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)"
   "+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1",
   "z", {Point(1.2038, 0, 0)}, 0.10, "blobby_plane_z_0.10");
run_polyline(
   "sqrt((x-1)*(x-1)+y*y+z*z)*sqrt((x+1)*(x+1)+y*y+z*z)*sqrt(x*x+(y-1)*(y-1)"
   "+z*z)*sqrt(x*x+(y+1)*(y+1)+z*z)-1.1",
   "z", {Point(1.2038, 0, 0)}, 0.05, "blobby_plane_z_0.05");
   */
  /*
    run_polyline("x^2+y^2+z^2-4", "x^2-y^2+z^2-1",
                 {Point(sqrt(5.0 / 2.0), sqrt(3.0 / 2.0), 0),
                  Point(sqrt(5.0 / 2.0), -sqrt(3.0 / 2.0), 0)},
                 0.45, "union", "sphere_hyperboloid_union_0.45");
    run_polyline("x^2+y^2+z^2-4", "x^2-y^2+z^2-1",
                 {Point(sqrt(5.0 / 2.0), sqrt(3.0 / 2.0), 0),
                  Point(sqrt(5.0 / 2.0), -sqrt(3.0 / 2.0), 0)},
                 0.30, "union", "sphere_hyperboloid_union_0.30");
    run_polyline("x^2+y^2+z^2-4", "x^2-y^2+z^2-1",
                 {Point(sqrt(5.0 / 2.0), sqrt(3.0 / 2.0), 0),
                  Point(sqrt(5.0 / 2.0), -sqrt(3.0 / 2.0), 0)},
                 0.15, "union", "sphere_hyperboloid_union_0.15");
    */
  /*
 run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
              0.45, "union", "sphere_sphere_union_0.45");
 run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
              0.30, "union", "sphere_sphere_union_0.30");
 run_polyline("(x-2)^2+y^2+z^2-4", "x^2+y^2+z^2-4", {Point(1, sqrt(3.0), 0)},
              0.15, "union", "sphere_sphere_union_0.15");
              */
  /*
run_polyline("x^2+y^2+z^2-9", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
  {Point(sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), -sqrt(1.75))},
  0.45, "union", "tanglecube_sphere__union_0.45");
run_polyline("x^2+y^2+z^2-9", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
  {Point(sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), -sqrt(1.75))},
  0.30, "union", "tanglecube_sphere__union_0.30");
run_polyline("x^2+y^2+z^2-9", "x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8",
  {Point(sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), sqrt(1.75)),
   Point(sqrt(4.5), -sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), sqrt(1.75), -sqrt(1.75)),
   Point(-sqrt(4.5), -sqrt(1.75), -sqrt(1.75))},
  0.15, "union", "tanglecube_sphere__union_0.15");*/
  /*
    run_polyline("x+1", "x^2+y^2+z^2-4", {Point(-1, sqrt(3), 0)}, 0.05,
                 "sphere_plane");*/
  //  run_all(3, 4, "/finite_surfaces/sphere", "my_run_input");
  // run_input(8, "/sing_surfaces/A2", "measure");
  // run_input(5, "/sing_surfaces/A3", "measure");
  // run_input(7, "/sing_surfaces/A4", "measure");
  // run_input(6, "/sing_surfaces/D4", "measure");
  // run_input(7, "/sing_surfaces/D4", "measure");
  // run_all(11, 11, "/sing_surfaces/D5", "measure");

  // run_all(9, 9, "/sing_surfaces/E6", "measure");

  // run_all(4, 4, "/sing_surfaces/E7", "measure");
  // run_all(4, 4, "/sing_surfaces/E8", "measure");
  //   run_all(3, 3, "/sing_surfaces/D4", "measure");

  // run_polyline(0.35);
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
