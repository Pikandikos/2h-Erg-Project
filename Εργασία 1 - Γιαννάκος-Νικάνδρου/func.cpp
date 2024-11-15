#include "./func.h"
#include <iostream>

// Reads the json folder given in the file path and distributes it to the vectors
bool read_json_file(const string &file_path, string &instance_uid, vector<Point_2> &points, vector<int> &region_boundary, int &num_constraints, vector<pair<int, int>> &additional_constraints)
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

    // Updated handling of "additional_constraints" as array
    for (const auto &constraint : pt.get_child("additional_constraints"))
    {
        num_constraints++;
        auto first = constraint.second.front().second.get_value<int>();
        auto second = constraint.second.back().second.get_value<int>();
        additional_constraints.emplace_back(first, second);
    }

    return true;
}

// γωνια που σχηματιζουν p1, p2, p3, ελεγχος για αμβλυγωνιο (επιστρεφει τιμη >90º)
double angle_between_points(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    Kernel::Vector_2 v1 = p2 - p1;                        // διανυσμα που ξεκιναει απο p1 προς p2
    Kernel::Vector_2 v2 = p3 - p1;                        // διανυσμα που ξεκιναει απο p1 προς p3
    double dot_product = v1 * v2;                         // εσωτερικο γινομενο
    double magnitude_v1 = std::sqrt(v1.squared_length()); // μέτρο διανυσματος
    double magnitude_v2 = std::sqrt(v2.squared_length()); // μέτρο διανυσματος

    // Ελεγχος για μηδενικο διανυσμα για αποφύγης διαίρεσης με το μηδέν
    if (magnitude_v1 == 0 || magnitude_v2 == 0)
    {
        return 0; // Μηδενικη γωνια αν και τα 2 σημεια ιδια
    }
    // Clamp the dot_product / (magnitude_v1 * magnitude_v2) to the range [-1, 1]
    double cosine_angle = dot_product / (magnitude_v1 * magnitude_v2);
    cosine_angle = std::max(-1.0, std::min(1.0, cosine_angle)); // Ensure valid range

    return std::acos(dot_product / (magnitude_v1 * magnitude_v2)) * 180.0 / M_PI; // υπολογισμος γωνιας
}

void check_cdt_validity(const CDT &cdt)
{
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        if (!fit->is_valid())
        {
            std::cerr << "Invalid face detected!" << endl;
        }
    }
    // cout << "All Faces valid for now..." << endl;
    return;
}

bool add_steiner_point_on_edge(CDT &cdt, const CDT::Edge &edge)
{
    // Check if the edge is valid
    if (!edge.first->is_valid())
    {
        cerr << "Invalid edge detected, skipping Steiner point insertion" << endl;
        return false;
    }

    CDT::Vertex_handle vh1 = edge.first->vertex((edge.second + 1) % 3);
    CDT::Vertex_handle vh2 = edge.first->vertex((edge.second + 2) % 3);

    // Check if the vertex handles are valid
    if (!vh1->is_valid() || !vh2->is_valid())
    {
        cerr << "Invalid vertex handle detected, skipping Steiner point insertion" << endl;
        return false; // Early exit to avoid inserting into an invalid edge
    }

    // Attempt to use the circumcenter as the Steiner point
    Point_2 circumcenter = CGAL::circumcenter(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point());
    bool circumcenter_exists = false;

    // Check if the circumcenter already exists as a vertex
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->point() == circumcenter)
        {
            circumcenter_exists = true;
            break;
        }
    }

    if (!circumcenter_exists)
    {
        // Try to insert the circumcenter
        CDT::Vertex_handle new_vertex = cdt.insert(circumcenter);
        if (new_vertex != nullptr)
        {
            cout << "Steiner Point (circumcenter) added: (" << circumcenter.x() << ", " << circumcenter.y() << ")" << endl;
            return true;
        }
    }

    // If the circumcenter could not be inserted or already exists, fallback to midpoint
    Point_2 midpoint = CGAL::midpoint(vh1->point(), vh2->point());
    bool midpoint_exists = false;

    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->point() == midpoint)
        {
            midpoint_exists = true;
            break;
        }
    }

    if (!midpoint_exists)
    {
        CDT::Vertex_handle new_vertex = cdt.insert(midpoint);
        if (new_vertex != nullptr)
        {
            cout << "Steiner Point (midpoint) added: (" << midpoint.x() << ", " << midpoint.y() << ")" << endl;
            return true;
        }
    }

    std::cerr << "Both circumcenter and midpoint Steiner points already exist, skipping insertion." << std::endl;
    return false;
}

bool attempt_to_flip(CDT &cdt, CDT::Finite_faces_iterator face_it, CDT::Edge edge)
{
    CDT copied_cdt = cdt; // Create a copy of the CDT

    // Ensure the edge has two distinct faces
    CDT::Face_handle face0 = edge.first;
    CDT::Face_handle face1 = face0->neighbor(edge.second);
    if (face1 == nullptr || face0 == nullptr || face0 == face1)
    {
        std::cerr << "Error: Edge does not have two distinct valid faces for flipping... Extiting attempt_to_flip" << std::endl;
        return false;
    }

    if (cdt.is_infinite(face0) || cdt.is_infinite(face1))
    {
        std::cerr << "One of the faces is infinite." << std::endl;
        return false; // Handle this case appropriately
    }

    // Perform the edge flip in the copied CDT
    copied_cdt.flip(face0, edge.second);
    check_cdt_validity(cdt);

    // Check angles in the copied CDT
    bool all_acute = true;
    for (CDT::Finite_faces_iterator fit = copied_cdt.finite_faces_begin(); fit != copied_cdt.finite_faces_end(); ++fit)
    {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        double angle1 = angle_between_points(p1, p2, p3);
        double angle2 = angle_between_points(p2, p1, p3);
        double angle3 = angle_between_points(p3, p1, p2);

        if (angle1 > 90 || angle2 > 90 || angle3 > 90)
        {
            all_acute = false;
            break;
        }
    }

    if (all_acute)
    {
        cout << "Flipping Edge!!!" << endl;
        cdt = copied_cdt; // Apply the flip if all angles are acute
        return true;
    }
    else
    {
        std::cerr << "Flip resulted in non-acute angles, reverting changes." << std::endl;
        return false; // Revert to the original CDT
    }
}

CDT triangulation(vector<Point_2> &points, vector<int> &region_boundary)
{
    // τριγωνοποίηση Delaunay
    CDT cdt;

    // προσθήκη σημείων από τον vector points
    for (const auto &point : points)
    {
        cdt.insert(point); // εισαγωγή σημείου στην τριγωνοποίηση
        check_cdt_validity(cdt);
    }

    // προσθήκη περιορισμένων ακμών (PSLG)
    for (std::size_t i = 0; i < region_boundary.size() - 1; i++ /*i += 2*/)
    { // ζεύγος δεικτών αναπαριστουν ακμή, ορια περιοχης που θα τριγωνοποιηθει
        cdt.insert_constraint(points[region_boundary[i]], points[region_boundary[i + 1]]);
    }

    check_cdt_validity(cdt);

    // επανάληψη για προσθήκη σημείων Steiner αν υπάρχουν αμβλυγώνια τρίγωνα
    bool all_acute = false;
    bool steiner_point_inserted;
    bool steiner_point_added_this_rotation;
    int no_of_steiner_points_added = 0;

    while (!all_acute)
    {
        all_acute = true; // υποθετω ολα τα τριγωνα οξυγωνια
        cout << endl
             << "Starting to check angles... again" << endl;

        int face_count = cdt.number_of_faces();
        cout << "Number of faces: " << face_count << endl;

        for (CDT::Finite_faces_iterator face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); face_it++)
        {
            steiner_point_added_this_rotation = false;
            if (!face_it->is_valid())
            {
                std::cerr << "Invalid face detected, skipping." << endl;
                continue; // Skip invalid faces
            }

            Point_2 p1 = face_it->vertex(0)->point();
            Point_2 p2 = face_it->vertex(1)->point();
            Point_2 p3 = face_it->vertex(2)->point();

            cout << "P1:" << p1 << "  P2:" << p2 << "  P3:" << p3 << endl;

            // Υπολογισμός των γωνιών του τριγώνου
            double angle1 = angle_between_points(p1, p2, p3);
            double angle2 = angle_between_points(p2, p1, p3);
            double angle3 = angle_between_points(p3, p1, p2);

            cout << "Angle1:" << angle1 << "  Angle2:" << angle2 << "  Angle3:" << angle3 << endl;

            if (angle1 == 0 || angle2 == 0 || angle3 == 0)
            {
                cout << "Angle = 0 Found !!!!" << endl;
                continue;
            }
            // if (angle1 > 90 || angle2 > 90 || angle3 > 90)
            // {
            //     all_acute = false;
            //     bool flipped = false;

            //     // Try to flip the edges in the order of angles, fallback to Steiner point if flipping fails
            //     if (angle1 > 90)
            //         flipped = attempt_to_flip(cdt, face_it, CDT::Edge(face_it, 2));

            //     if (!flipped && angle2 > 90)
            //         flipped = attempt_to_flip(cdt, face_it, CDT::Edge(face_it, 0));

            //     if (!flipped && angle3 > 90)
            //         flipped = attempt_to_flip(cdt, face_it, CDT::Edge(face_it, 1));

            //     if (!flipped) // If flipping fails, add a Steiner point
            //     {
            if (angle1 > 90)
            {
                all_acute = false;
                bool flipped = false;
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 2));
                if (steiner_point_inserted) // steiner point hasn't been skipped
                {
                    no_of_steiner_points_added++;
                    cout << "No. of Steiner Points: " << no_of_steiner_points_added << endl;
                    face_it = cdt.finite_faces_begin();
                    steiner_point_added_this_rotation = true;
                    break; // Steiner Point has been added... Now Restart the checking of angles
                }
            }
            else if (angle2 > 90)
            {
                all_acute = false;
                bool flipped = false;
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 0));
                if (steiner_point_inserted) // steiner point hasn't been skipped
                {
                    no_of_steiner_points_added++;
                    cout << "No. of Steiner Points: " << no_of_steiner_points_added << endl;
                    face_it = cdt.finite_faces_begin();
                    steiner_point_added_this_rotation = true;
                    break;
                }
            }
            else if (angle3 > 90)
            {
                all_acute = false;
                bool flipped = false;
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 1));
                if (steiner_point_inserted) // steiner point hasn't been skipped
                {
                    no_of_steiner_points_added++;
                    cout << "No. of Steiner Points: " << no_of_steiner_points_added << endl;
                    face_it = cdt.finite_faces_begin();
                    steiner_point_added_this_rotation = true;
                    break;
                }
            }
            // }
            // }
        }
        if (steiner_point_added_this_rotation == false) // Will need to add flip if it was working
        {
            cout << "All faces/triangles are acute" << endl;
            all_acute = true; // If the loop ever ends it means that all faces are acute
        }
    }
    return cdt;
}

// Function to create JSON output from the CDT and save it to a file
void create_json_output(const CDT &cdt, const std::string &filename)
{
    using boost::property_tree::ptree;
    ptree json_output;

    // Set static fields
    json_output.put("content_type", "CG_SHOP_2025_Solution");
    json_output.put("instance_uid", "unique_instance_id");

    // Create arrays for steiner_points_x and steiner_points_y
    ptree steiner_points_x, steiner_points_y;
    for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); vit++)
    {
        ptree x, y;
        x.put("", vit->point().x());
        y.put("", vit->point().y());

        steiner_points_x.push_back(std::make_pair("", x));
        steiner_points_y.push_back(std::make_pair("", y));
    }

    json_output.add_child("steiner_points_x", steiner_points_x);
    json_output.add_child("steiner_points_y", steiner_points_y);

    // Create edges array
    ptree edges;
    for (CDT::Finite_edges_iterator edge_it = cdt.finite_edges_begin(); edge_it != cdt.finite_edges_end(); edge_it++)
    {
        ptree edge_array;

        // CDT::Edge e1 = edges.first->vertex(cdt.ccw(edge_it->second)); // First vertex of the edge
        // CDT::Edge e2 = edges.first->vertex(cdt.cw(edge_it->second));  // Second vertex of the edge
        // edge_array.push_back(std::make_pair("", )); // Assume 'info' provides a unique index
        // edge_array.push_back(std::make_pair("", ));

        // Get the vertices of the edge
        CDT::Vertex_handle v1 = edge_it->first->vertex(cdt.ccw(edge_it->second));
        CDT::Vertex_handle v2 = edge_it->first->vertex(cdt.cw(edge_it->second));

        // Access the points of the vertices
        Point_2 p1 = v1->point(); // Get the point of the first vertex
        Point_2 p2 = v2->point(); // Get the point of the second vertex

        // Store vertex coordinates in the edge array
        edge_array.put("vertex1.x", p1.x());
        edge_array.put("vertex1.y", p1.y());
        edge_array.put("vertex2.x", p2.x());
        edge_array.put("vertex2.y", p2.y());

        edges.push_back(std::make_pair("", edge_array));
    }
    json_output.add_child("edges", edges);

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