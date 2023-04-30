
#include <cassert>
#include <ginac/ginac.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "src/assertm.h"
#include "src/point.h"
#include "src/vector.h"
#include "src/edge.h"
#include "src/triangle.h"
#include "src/function.h"
#include "src/algorithms.h"
#include "src/bounding_box.h"

using namespace std;
using namespace GiNaC;

int is_in_points(const Point P, const vector<pair<Point, vector<int> > > &points){
    for(size_t i = 0; i < points.size(); ++i){
        if(points[i].first == P)
            return i;
    }
    return -1;
}

int is_in_edges(const Edge & E, const vector<pair<Edge, vector<Triangle> > > &edges){
    for (size_t i = 0; i < edges.size(); ++i){
        if(edges[i].first == E || edges[i].first == Edge(E.B(), E.A())){
            return i;
        }
    }
    return -1;
}
bool is_in_pointers(const int p, const vector<int>&pointers){
    for(auto i : pointers){
        if(i==p){
            return true;
        }
    }
    return false;
}

bool is_in_triangles(const Triangle T, const vector<Triangle> &triangles){
    for(auto t:triangles){
        if(T == t){
            return true;
        }
    }
    return false;
}

vector<numeric> average_gc_distance(const vector<Triangle> &mesh_triangles, 
                                    const Function &F, const numeric e_size) {
  int n = mesh_triangles.size();
  vector<numeric>result;
  numeric sum = 0;
  numeric min = 1000*e_size;
  numeric max = -1;
  for(auto T: mesh_triangles){
    Point gc = T.get_gravity_center();
    Vector direction = F.get_gradient_at_point(gc).unit();
    std::optional<Point> gc_proj = project(gc, direction, F, e_size);
    assertm (gc_proj.has_value(), "No value!");
    numeric dist = Vector(gc, gc_proj.value()).get_length();
    sum += dist;
    if(dist<min) min = dist;
    if(dist>max) max = dist;
  }
  result.push_back(sum/n);
  result.push_back(max);
  result.push_back(min);
  return result;
}

numeric average_side_length(const vector<pair<Edge, vector<Triangle>>> &mesh_edges) {
  int n = mesh_edges.size();
  numeric sum = 0;
  for(auto E: mesh_edges){
    sum += E.first.get_length();
  }
  return sum/n;
}
numeric average_side_ratio(const vector<Triangle> &mesh_triangles) {
  int n = mesh_triangles.size();
  numeric sum = 0;
  for(auto T: mesh_triangles){
    numeric a = T.BC().get_length();
    numeric b = T.CA().get_length();
    numeric c = T.AB().get_length();

    numeric biggest = std::max(a, std::max(b, c));
    numeric smallest = std::min(a, std::min(b, c));
    numeric ratio = biggest/smallest;
    sum = sum+ratio;
  }
  return sum/n;
}

pair<numeric, numeric> min_max_normal_angle(const vector<pair<Edge, vector<Triangle>>> &mesh_edges, const Function &F, const numeric e_size){
    numeric min_angle = 10;
    numeric max_angle =-1;
    bool bounding_edge = false;
    for(auto E : mesh_edges)
    {
        Edge e = E.first;
        vector<Triangle> T = E.second;
        if(T.size() == 2){
            Vector n1 = F.outside_normal(T[0], e_size);
            Vector n2 = F.outside_normal(T[1], e_size);
            numeric angle = acos(n1*n2);
            if(angle>max_angle) max_angle = angle;
            if(angle<min_angle) min_angle = angle;
            assertm(min_angle<=max_angle, "Min angle greater than max angle!");
        }
        else if(T.size() == 1)
            bounding_edge = true;
        else{
            cout<<"Weird number of neigghbour triangles: " << T.size() << endl;
            //assertm(false, "Wrong number of neighbours in mesh!");
        }
        
    }
    if(bounding_edge) cout<< "Bounding edge in a model!"<<endl;
    return pair(min_angle, max_angle);
}

pair<numeric, numeric> mean_std(const vector< pair<Point, vector<int> > > & points){
    numeric sum_avg = 0;
    numeric sum_std = 0;
    for (auto p : points){
        if(p.second.size() == 0) continue;
        Point P = p.first;
        numeric sum = 0;
        for(auto index : p.second){
            numeric dist = Vector(P, points[index].first).get_length();
            sum += dist;
        }
        numeric avg = sum/p.second.size();
        sum = 0;
        for (auto index : p.second){
            numeric clen = (Vector(P, points[index].first).get_length() - avg);
            clen *= clen; 
            sum += clen;
        }
        numeric std = sqrt(sum/p.second.size());
        sum_avg += avg;
        sum_std += std;
    }
    numeric avg_avg = sum_avg/points.size();
    numeric std_avg = sum_std/points.size();
    return pair(avg_avg, std_avg);
}




void measure (const vector<pair<Point, vector<int> > > &mesh_points, 
            const vector<Triangle> &mesh_triangles, 
            const vector<pair<Edge, vector<Triangle>> > &mesh_edges, const Function &F, 
            const numeric e_size, const string folder, const string &name){
    string output_dir = "./measure/measure_data/" + folder + "/" + name + ".out";
    std::cout<<"Output directory: " << output_dir << endl;

    std::ofstream out(output_dir);
    cout<<"Started measuring..."<<endl;
    double avg_side_length = to_double(average_side_length(mesh_edges));
    double avg_max_side_ratio = to_double(average_side_ratio(mesh_triangles));

    auto avg_max_min_gc_dist = average_gc_distance(mesh_triangles, F, e_size);
    double avg_gc_dist = to_double(avg_max_min_gc_dist[0]);
    double max_gc_dist = to_double(avg_max_min_gc_dist[1]);
    //double min_gc_dist = to_double(avg_max_min_gc_dist[2]);
    //numeric avg_gc_dist = average_gc_distance(mesh_triangles, F, e_size);
    //numeric ideal_triangle_area = e_size*e_size*sqrt(numeric(3))/4;
    auto avg_mean_std = mean_std(mesh_points);
    double avg_mean_neighbour_points_dist = to_double(avg_mean_std.first);
    double avg_std_neighbour_points_dist = to_double(avg_mean_std.second);
    //auto min_max_angle = min_max_normal_angle(mesh_edges, F, e_size);

    //double min_normal_angle = to_double(min_max_angle.first*100);
    //double max_normal_angle = to_double(min_max_angle.second);
    double d_e_size = to_double(e_size);

    cout<<"Measure done."<<endl;

    Digits = 15;


    out << std::fixed;
    out << std::setprecision(3) << avg_max_side_ratio << " & " << max_gc_dist << " & " <<
    avg_gc_dist << " & " << 
    avg_mean_neighbour_points_dist << " & " << 
    avg_std_neighbour_points_dist << "\\\\" << endl << endl;


    out << std::fixed;
    out << std::setprecision(3) << avg_max_side_ratio << " & " << max_gc_dist << " & " <<
    avg_gc_dist << " & " << 
    avg_std_neighbour_points_dist << "\\\\" << endl << endl;
/*
    out << std::setprecision(3) << d_e_size << " & " << avg_side_length/d_e_size << " & " << avg_gc_dist/d_e_size << " & " <<
    avg_max_side_ratio << " & " << max_gc_dist/d_e_size << " & " << 
    avg_mean_neighbour_points_dist/d_e_size << " & " << 
    avg_std_neighbour_points_dist/d_e_size << "\\\\" << endl << endl;
*/
    out<<avg_gc_dist/avg_side_length<<endl<<endl;

    out << "DATA:" << endl;
    out << "Edge size: " << e_size << endl;
    //out << "Area of triangle with edge size " << e_size << " is: " << ideal_triangle_area << endl;
    out << "Number of triangles in mesh: " << mesh_triangles.size() << endl;
    out << "Number of points in mesh: " << mesh_points.size() << endl << endl;
    
    out << "AVERAGES: " << endl;
    out << "Average side length: " << avg_side_length << endl;
    out << "Average maximum triangle side ratio: " << avg_max_side_ratio << endl;
    out << "Average gravity center distance: " << avg_gc_dist << endl;
    out << "Average of mean of nieghbour points distance: " << avg_mean_neighbour_points_dist << endl;
    out << "Average of standard deviation of nieghbour points distance: " << avg_std_neighbour_points_dist << endl << endl;
//   out << "Average area of triangles: " << avg_tri_area << endl << endl;

    out << "MIN and MAX: " << endl;
    out << "Hausdorff distance - Max gc distance: " << max_gc_dist << endl;
    //out << "Min neighbour triangles normals angle: " << min_normal_angle << endl;
    //out << "Max neighbour triangles normals angle: " << max_normal_angle << endl << endl;

    out << "RATIOS:" << endl;
//   out << "Average trinagle area to ideal triangle area: " << avg_tri_area/ideal_triangle_area << endl;
    out << "Average gravity center distance to edge size: " << avg_gc_dist/e_size << endl; 
    out << "Avg of mean of neighbour points distance to edge size: " << avg_mean_neighbour_points_dist/e_size << endl;
    out << "Avg of std of nieghbour points distance to edge size: " << avg_std_neighbour_points_dist/e_size << endl;
    out << "Average edge size to wished edge size: " << avg_side_length/e_size << endl;
    out << "Max gravity center distance to edge size: " << max_gc_dist/e_size << endl;

    std::cout<<"Output done."<<endl;
    return;
}

void load_data_old_obj(std::ifstream& obj_file, vector<pair<Point,vector<int>> >&mesh_points, vector<Triangle>&mesh_triangles, vector<pair<Edge, vector<Triangle>> >&mesh_edges, const BoundingBox& my_bounding_box) {
    vector<Edge>bounding_edges;
    string type;
    int round = 0;
    vector<Point>last_points;

    for(string type; obj_file >> type;){
        if(type == "v")
        {
            double a, b, c;
            obj_file >> a >> b >> c;
            numeric n_a, n_b, n_c;
            n_a = ex_to<numeric>(a);
            n_b = ex_to<numeric>(b);
            n_c = ex_to<numeric>(c);

            Point P(n_a, n_b, n_c);
            
            last_points.push_back(P);
            if(round == 2)
            {
                Triangle T(last_points[0], last_points[1], P);
                last_points.clear();
                mesh_triangles.push_back(T);
                int edge_index0 = is_in_edges(T.AB(), mesh_edges);
                int edge_index1 = is_in_edges(T.BC(), mesh_edges);
                int edge_index2 = is_in_edges(T.CA(), mesh_edges);
                if(edge_index0 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.AB(), Tvec));
                }
                else{
                        mesh_edges[edge_index0].second.push_back(T);
                    }
                
                if(edge_index1 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.BC(), Tvec));
                }
                else{
                    mesh_edges[edge_index1].second.push_back(T);
                }
                if(edge_index2 == -1){
                    vector<Triangle>Tvec;
                    Tvec.push_back(T);
                    mesh_edges.push_back(pair(T.CA(), Tvec));
                }
                else{
                    mesh_edges[edge_index2].second.push_back(T);
                }
                if(my_bounding_box.is_new_bounding_edge(T.AB())){
                    bounding_edges.push_back(T.AB());
                }
                if(my_bounding_box.is_new_bounding_edge(T.BC())){
                    bounding_edges.push_back(T.BC());
                }
                if(my_bounding_box.is_new_bounding_edge(T.CA())){
                    bounding_edges.push_back(T.CA());
                }
            }

            int ind = is_in_points(P, mesh_points);
            if(ind == -1){
                pair<Point, vector<int> > par = pair(P, vector<int>());
                mesh_points.push_back(par);
            }
        round = (round+1)%3;
        }
    }
    for(Triangle T : mesh_triangles){
        int indA = is_in_points(T.A(), mesh_points);
        int indB = is_in_points(T.B(), mesh_points);
        int indC = is_in_points(T.C(), mesh_points);
        
        if(!is_in_pointers(indB, mesh_points[indA].second))
            mesh_points[indA].second.push_back(indB);
        if(!is_in_pointers(indC, mesh_points[indA].second))
            mesh_points[indA].second.push_back(indC);
        
        if(!is_in_pointers(indA, mesh_points[indB].second))
            mesh_points[indB].second.push_back(indA);
        if(!is_in_pointers(indC, mesh_points[indB].second))
            mesh_points[indB].second.push_back(indC);

        if(!is_in_pointers(indA, mesh_points[indC].second))
            mesh_points[indC].second.push_back(indA);
        if(!is_in_pointers(indB, mesh_points[indC].second))
            mesh_points[indC].second.push_back(indB);
    }
}

void load_data_new_obj(std::ifstream& obj_file, vector<pair<Point,vector<int>> >&mesh_points, vector<Triangle>&mesh_triangles, vector<pair<Edge, vector<Triangle>> >&mesh_edges, const BoundingBox& my_bounding_box) {
    // vector<pair<Point,vector<int>> > contains all points and indices to all neighbour points
    // vector<Triangle>mesh_triangles contains all mesh triangles
    // vector<pair<Edge, vector<Triangle>> >mesh_edges contains each mesh edge once and both neighbour triangles
    //int counter = 0;
    for(string type; obj_file >> type;){
        //std::cout<<counter<<endl;
        //counter++;
        if(type == "v")
        {
            double a, b, c;
            obj_file >> a >> b >> c;
            numeric n_a, n_b, n_c;
            n_a = ex_to<numeric>(a);
            n_b = ex_to<numeric>(b);
            n_c = ex_to<numeric>(c);

            Point P(n_a, n_b, n_c); // meshpoint

            pair<Point, vector<int> > par = pair(P, vector<int>());
            mesh_points.push_back(par);
        }
        // face
        else if (type == "f") {
            int i0, i1, i2;
            obj_file >> i0 >> i1 >> i2;
            Point A = mesh_points[i0-1].first;
            Point B = mesh_points[i1-1].first;
            Point C = mesh_points[i2-1].first;
            if (A==B || B==C || A==C) continue;
            Triangle T(A, B, C);
            if(!T.is_triangle()) continue;
            mesh_triangles.push_back(T);
            int indA = i0-1;
            int indB = i1-1;
            int indC = i2-1;
            
            if(!is_in_pointers(indB, mesh_points[indA].second))
                mesh_points[indA].second.push_back(indB);
            if(!is_in_pointers(indC, mesh_points[indA].second))
                mesh_points[indA].second.push_back(indC);
            
            if(!is_in_pointers(indA, mesh_points[indB].second))
                mesh_points[indB].second.push_back(indA);
            if(!is_in_pointers(indC, mesh_points[indB].second))
                mesh_points[indB].second.push_back(indC);

            if(!is_in_pointers(indA, mesh_points[indC].second))
                mesh_points[indC].second.push_back(indA);
            if(!is_in_pointers(indB, mesh_points[indC].second))
                mesh_points[indC].second.push_back(indB);
            int edge_index0 = is_in_edges(T.AB(), mesh_edges);
            int edge_index1 = is_in_edges(T.BC(), mesh_edges);
            int edge_index2 = is_in_edges(T.CA(), mesh_edges);
            if(edge_index0 == -1){
                vector<Triangle>Tvec;
                Tvec.push_back(T);
                mesh_edges.push_back(pair(T.AB(), Tvec));
            }
            else {
                mesh_edges[edge_index0].second.push_back(T);
            }
            
            if(edge_index1 == -1) {
                vector<Triangle>Tvec;
                Tvec.push_back(T);
                mesh_edges.push_back(pair(T.BC(), Tvec));
            }
            else {
                mesh_edges[edge_index1].second.push_back(T);
            } 
            if(edge_index2 == -1) {
                vector<Triangle>Tvec;
                Tvec.push_back(T);
                mesh_edges.push_back(pair(T.CA(), Tvec));
            }
            else {
                mesh_edges[edge_index2].second.push_back(T);
            }
        }
    }
    return;
}

//measures given input, creates output in corresponding folder inside "/measure_data" folder
void run_input(const int i, const string folder, const string index){

    //string index_str = index;

    string file_name0 = "./measure/inputs" + folder + "/input" + to_string(i);
    std::ifstream input_file (file_name0, std::ifstream::in);
    assertm(input_file.is_open(), "Failed opening the first input file!");
    vector<string>parsed_input;
    for(string str; getline(input_file, str);){
        string s;
        bool writing = false;
        for(char c : str){
        if(writing && c!='"'){
            s.push_back(c);
        }
        if(c=='"') writing = !writing;
        }
        parsed_input.push_back(s);
    }
    realsymbol x("x"), y("y"), z("z");

    symtab table;
    table["x"] = x;
    table["y"] = y;
    table["z"] = z;
    
    parser reader(table);
    string name;
    if(index == "singsurf")
        name = parsed_input[0];
    else if(index == "uniform" || index == "adaptive")
        name = parsed_input[0] + "_" + parsed_input[1];

    ex input_F = reader(parsed_input[2]);
    
    numeric min_x = ex_to<numeric>(stod(parsed_input[3]));
    numeric max_x = ex_to<numeric>(stod(parsed_input[4]));
    numeric min_y = ex_to<numeric>(stod(parsed_input[5]));
    numeric max_y = ex_to<numeric>(stod(parsed_input[6]));
    numeric min_z = ex_to<numeric>(stod(parsed_input[7]));
    numeric max_z = ex_to<numeric>(stod(parsed_input[8]));
    BoundingBox my_bounding_box(min_x, max_x, min_y, max_y, min_z, max_z);

    numeric e_size = ex_to<numeric>(stod(parsed_input[9]));

    vector<ex> input_dF;
    input_dF.push_back(diff(input_F, x));
    input_dF.push_back(diff(input_F, y));
    input_dF.push_back(diff(input_F, z));

    Function F(x, y, z, input_F, input_dF);
    cout<<"Name:" << name << endl;
    string file_name;
    if (index == "uniform")
        file_name = "/uniform/" + name + "_uniform" + ".obj";
    else if(index == "adaptive")
        file_name = "/adaptive/" + name + "_adaptive" + ".obj";
    else if(index == "singsurf")
        file_name = "/singsurf/" + name + ".obj";

    string file_path = "./measure/outputs" + file_name;
    std::cout<< "Path: "<< file_path<<endl;

    std::ifstream obj_file (file_path, std::ifstream::in);
    assertm(obj_file.is_open(), "Failed opening the second input file!");

    vector<pair<Point, vector<int> > >mesh_points;
    vector<Triangle>mesh_triangles;
    vector<pair<Edge, vector<Triangle>> >mesh_edges;
    std::cout<<"Started loading data!"<<endl;
    //load_data_old_obj(obj_file, mesh_points, mesh_triangles, mesh_edges, my_bounding_box);
    load_data_new_obj(obj_file, mesh_points, mesh_triangles, mesh_edges, my_bounding_box);
    cout<<"Succesfuly loaded data!" << endl;

    measure(mesh_points, mesh_triangles, mesh_edges, F, e_size, index, name);
    
}

//measures multiple inputs, creates output in corresponding folder inside "/measure_data" folder
void run_all(const int beg, const int end, const string folder, const string name){
    for(int i = beg; i<=end; ++i){
        run_input(i, folder, name);
    }
}

int main(){

    // spustame prikazom "make run_measure" 

    // zmeria model vyprodukovany vstupnym suborom "input0" v priecinku "/measure/inputs/finite_surfaces/sphere",
    // ocakava, ze model sa nachadza v priecinku "/measure/outputs/finite_surfaces/sphere" a ze 
    // bol vytvoreny s predponou "measure" (predpona, ktoru si volime pri spustani algoritmu v main.cpp)
    // vystupny subor vlozi do priecinka "/measure/measure_data/finite_surfaces/sphere" a nazve ho s 
    // predponou "measure"
    run_input(4, "/sing_surfaces/A1", "adaptive");
    /*run_all(3, 3, "/sing_surfaces/A1", "singsurf");
    run_all(3, 3, "/sing_surfaces/A1", "me");
    run_all(4, 5, "/sing_surfaces/A2", "singsurf");
    run_all(4, 5, "/sing_surfaces/A2", "me");
    
    run_all(2, 3, "/sing_surfaces/A3", "singsurf");
    run_all(2, 3, "/sing_surfaces/A3", "me");
    run_all(2, 3, "/sing_surfaces/A4", "singsurf");
    run_all(2, 3, "/sing_surfaces/A4", "me");
    */
    
    //run_all(2, 2, "/sing_surfaces/D4", "singsurf");
    //run_all(2, 2, "/sing_surfaces/D4", "me");
    //run_all(6, 7, "/sing_surfaces/D5", "singsurf");
    //run_all(6, 7, "/sing_surfaces/D5", "me");
    /*
    run_all(4, 5, "/sing_surfaces/E6", "singsurf");
    run_all(4, 5, "/sing_surfaces/E6", "me");
    run_all(2, 2, "/sing_surfaces/E7", "singsurf");
    run_all(2, 2, "/sing_surfaces/E7", "me");
    run_all(2, 2, "/sing_surfaces/E8", "singsurf");
    run_all(2, 2, "/sing_surfaces/E8", "me");
    */
    // zmeria modely vyprodukovane vstupnymi suboromi "input1", "input2", "input3" v priecinku 
    // "/measure/inputs/adaptive_height/neighbour_influence_adapt", ocakava, ze model sa nachadza 
    // v priecinku "/measure/outputs/adaptive_height/neighbour_influence_adapt" a ze 
    // bol vytvoreny s predponou "adaptive" (predpona, ktoru si volime pri spustani algoritmu v main.cpp)
    // vystupny subor vlozi do priecinka "/measure/measure_data/adaptive_height/neighbour_influence_adapt" 
    // a nazve ho s predponou "adaptive"
    //run_all(1, 3, "/adaptive_height/neighbour_influence_adapt", "adaptive");
}