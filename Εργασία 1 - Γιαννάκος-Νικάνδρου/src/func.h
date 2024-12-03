#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <sstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
// #include <boost/json.hpp>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Kernel/global_functions.h>

using namespace boost::property_tree;

using Kernel = CGAL::Simple_cartesian<double>;
using Point_2 = Kernel::Point_2;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

// using DelTr = CGAL::Delaunay_triangulation_2<Kernel>; // Probably not needed
using CDT = CGAL::Constrained_Delaunay_triangulation_2<Kernel>;

// contstraint delaunay triangulation
typedef CGAL::Triangulation_vertex_base_with_info_2<int, Kernel> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Edge Edge;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

using namespace std;

const int MAX_NO__STEINER_POINTS = 1;

// export.cpp
void export_to_svg(const CDT &cdt, const std::string &filename);

// trianglulation.cpp
CDT triangulation(vector<Point_2> &points, vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints);
bool add_steiner_point_on_edge(CDT &cdt, const CDT::Edge &edge, const vector<pair<Point_2, Point_2>> &constraints);
bool attempt_to_flip(CDT &cdt, CDT::Finite_faces_iterator face_it, CDT::Edge edge);

// common.cpp
double angle_between_points(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3);
void check_cdt_validity(const CDT &cdt);
bool is_point_inside_constraints(const Point_2 &point, const vector<pair<Point_2, Point_2>> &constraints);
void analyze_obtuse_angles(const CDT &cdt);
Point_2 mean_point_of_adjacent_triangles(CDT &cdt, CDT::Face_handle face, const vector<pair<Point_2, Point_2>> &constraints);

// io.c
bool read_json_file(const string &file_path, string &instance_uid, vector<Point_2> &points, vector<int> &region_boundary, int &num_constraints, vector<pair<int, int>> &additional_constraints,
                    string &method, ptree &parameters, bool &delaunay);
void create_json_output(const CDT &cdt, const std::string &filename);

#endif