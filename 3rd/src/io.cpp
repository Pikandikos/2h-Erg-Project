#include "./func.h"
#include <iostream>

bool read_json_file(const string &file_path, string &instance_uid, vector<Point_2> &points, vector<int> &region_boundary, int &num_constraints, vector<pair<int, int>> &additional_constraints, string &method, ptree &parameters, bool &delaunay)
{
    ptree pt;
    try
    {
        read_json(file_path, pt);
    }
    catch (const json_parser_error &err)
    {
        cerr << "Error parsing JSON file: " << err.what() << endl;
        return false;
    }

    instance_uid = pt.get<string>("instance_uid");

    ptree points_x = pt.get_child("points_x");
    ptree points_y = pt.get_child("points_y");

    auto x_it = points_x.begin();
    auto y_it = points_y.begin();

    while (x_it != points_x.end() && y_it != points_y.end())
    {
        points.emplace_back(x_it->second.get_value<int>(), y_it->second.get_value<int>());
        ++x_it;
        ++y_it;
    }

    for (const auto &boundary : pt.get_child("region_boundary"))
    {
        region_boundary.push_back(boundary.second.get_value<int>());
    }

    for (const auto &constraint : pt.get_child("additional_constraints"))
    {
        auto first = constraint.second.front().second.get_value<int>();
        auto second = constraint.second.back().second.get_value<int>();
        additional_constraints.emplace_back(first, second);
    }

    num_constraints = pt.get<int>("num_constraints");
    method = pt.get<string>("method");

    // Load parameters specific to the chosen method
    string parameters_key = "parameters_" + method;
    if (pt.find(parameters_key) != pt.not_found())
    {
        parameters = pt.get_child(parameters_key);
    }
    else
    {
        cerr << "Error: Parameters for method \"" << method << "\" not found in JSON file." << endl;
        return false;
    }

    delaunay = pt.get<bool>("delaunay");

    return true;
}

// Function to create JSON output from the CDT and save it to a file
#include <sstream> // For stringstream
void create_json_output(const CDT &cdt, const std::string &filename,
                        std::string instance_uid,
                        int obtuse_count, const std::string &method,
                        const std::vector<double> &method_parameters)
{
    using boost::property_tree::ptree;
    ptree json_output;

    // Static fields
    json_output.put("content_type", "CG_SHOP_2025_Solution");
    json_output.put("instance_uid", instance_uid);

    // Create arrays for steiner_points_x and steiner_points_y
    ptree steiner_points_x, steiner_points_y;
    for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); vit++)
    {
        ptree x, y;
        std::stringstream ss_x, ss_y;
        ss_x << vit->point().x();
        ss_y << vit->point().y();

        x.put("", ss_x.str());
        y.put("", ss_y.str());

        steiner_points_x.push_back(std::make_pair("", x));
        steiner_points_y.push_back(std::make_pair("", y));
    }

    json_output.add_child("steiner_points_x", steiner_points_x);
    json_output.add_child("steiner_points_y", steiner_points_y);

    // Create edges array
    ptree edges;
    for (CDT::Finite_edges_iterator edge_it = cdt.finite_edges_begin(); edge_it != cdt.finite_edges_end(); edge_it++)
    {
        CDT::Vertex_handle v1 = edge_it->first->vertex(cdt.ccw(edge_it->second));
        CDT::Vertex_handle v2 = edge_it->first->vertex(cdt.cw(edge_it->second));

        ptree edge;
        edge.push_back(std::make_pair("", ptree(std::to_string(v1->point().x()))));
        edge.push_back(std::make_pair("", ptree(std::to_string(v2->point().x()))));
        edges.push_back(std::make_pair("", edge));
    }
    json_output.add_child("edges", edges);

    // Add obtuse_count
    json_output.put("obtuse_count", obtuse_count);

    // Add method
    json_output.put("method", method);

    // Add parameters based on the method
    ptree params;
    if (method == "local" && method_parameters.size() == 1)
    {
        params.put("L", method_parameters[0]);
    }
    else if (method == "sa" && method_parameters.size() == 3)
    {
        params.put("alpha", method_parameters[0]);
        params.put("beta", method_parameters[1]);
        params.put("L", method_parameters[2]);
    }
    else if (method == "ant" && method_parameters.size() == 7)
    {
        params.put("alpha", method_parameters[0]);
        params.put("beta", method_parameters[1]);
        params.put("xi", method_parameters[2]);
        params.put("psi", method_parameters[3]);
        params.put("lambda", method_parameters[4]);
        params.put("kappa", method_parameters[5]);
        params.put("L", method_parameters[6]);
    }
    else
    {
        std::cerr << "Error: Incorrect number of parameters for the method " << method << std::endl;
        return;
    }
    json_output.add_child("parameters", params);

    // Write JSON output to a file
    std::ofstream output_file(filename);
    if (output_file.is_open())
    {
        write_json(output_file, json_output);
        output_file.close();
    }
    else
    {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
}
