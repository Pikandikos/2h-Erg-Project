#include "./func.h"
#include <iostream>

bool add_steiner_point_on_edge(CDT &cdt, const CDT::Edge &edge, const vector<pair<Point_2, Point_2>> &constraints)
{
    // Check if the edge is valid
    if (!edge.first->is_valid())
    {
        cerr << "Invalid edge detected, skipping Steiner point insertion" << endl;
        return false;
    }

    CDT::Vertex_handle vh1 = edge.first->vertex((edge.second + 1) % 3);
    CDT::Vertex_handle vh2 = edge.first->vertex((edge.second + 2) % 3);

    // Validate the face handle and vertices
    if (!edge.first->is_valid() || vh1 == nullptr || vh2 == nullptr)
    {
        cerr << "Invalid edge or vertices detected, skipping Steiner point insertion" << endl;
        return false; // Early exit
    }

    // // Compute the Steiner point from the neighboring faces
    // Point_2 polygon_st = compute_steiner_point_from_neighbors(cdt, *edge.first);

    // // Check if the computed Steiner point already exists in the triangulation
    // bool polygon_st_exists = false;
    // for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    // {
    //     if (vit->point() == polygon_st)
    //     {
    //         polygon_st_exists = true;
    //         break;
    //     }
    // }

    // if (!polygon_st_exists)
    // {
    //     // Insert the Steiner point from neighboring faces
    //     CDT::Vertex_handle new_vertex = cdt.insert(polygon_st);
    //     if (new_vertex != nullptr)
    //     {
    //         cout << "Steiner Point (from neighbors) added: (" << polygon_st.x() << ", " << polygon_st.y() << ")" << endl;
    //         return true;
    //     }
    // }

    // Check for degeneracy before circumcenter calculation (εκφυλισμένη κορυφή)
    if (CGAL::collinear(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point()))
    {
        cerr << "Degenerate triangle detected, skipping circumcenter calculation" << endl;
        return false; // Early exit to avoid inserting into an invalid edge
    }

    // Attempt to use the circumcenter as the Steiner point
    Point_2 circumcenter = CGAL::circumcenter(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point());
    bool circumcenter_not_inside = false;
    bool circumcenter_exists = false;

    // Validate if the circumcenter is within constraints
    if (!is_point_inside_constraints(circumcenter, constraints))
    {
        std::cerr << "Circumcenter is outside the constraints, skipping insertion." << endl;
        circumcenter_not_inside = false;
    }

    // Check if the circumcenter already exists as a vertex
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->point() == circumcenter)
        {
            circumcenter_exists = true;
            break;
        }
    }

    if (!(circumcenter_exists || circumcenter_not_inside))
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
    bool midpoint_not_inside = false;
    bool midpoint_exists = false;

    // Validate if the circumcenter is within constraints
    if (!is_point_inside_constraints(midpoint, constraints))
    {
        std::cerr << "Midpoint is outside the constraints too, skipping insertion." << endl;
        midpoint_not_inside = true;
    }

    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->point() == midpoint)
        {
            midpoint_exists = true;
            break;
        }
    }

    if (!(midpoint_exists || midpoint_not_inside))
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

CDT triangulation(vector<Point_2> &points, vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints)
{
    // τριγωνοποίηση Delaunay
    CDT cdt;

    // Manually store the constraints as pairs of points
    vector<std::pair<Point_2, Point_2>> constraints;

    // προσθήκη περιορισμένων ακμών (PSLG)
    // Ensure region boundary constraints are added in CCW order
    for (std::size_t i = 0; i < region_boundary.size() - 1; ++i)
    {
        Point_2 p1 = points[region_boundary[i]];
        Point_2 p2 = points[region_boundary[i + 1]];
        cdt.insert_constraint(p1, p2);
        constraints.push_back({p1, p2}); // Store the region boundary constraint
    }

    // Close the loop for region_boundary if needed
    if (region_boundary.size() > 2)
    {
        Point_2 p1 = points[region_boundary.back()];
        Point_2 p2 = points[region_boundary.front()];
        cdt.insert_constraint(p1, p2);
        constraints.push_back({p1, p2}); // Store the closing constraint
    }

    // Add additional constraints
    for (const auto &constraint : additional_constraints)
    {
        Point_2 p1 = points[constraint.first];
        Point_2 p2 = points[constraint.second];
        cdt.insert_constraint(p1, p2);
        constraints.push_back({p1, p2}); // Store the additional constraint
    }

    check_cdt_validity(cdt);

    // προσθήκη σημείων από τον vector points με έλεγχο των constraints
    for (const auto &point : points)
    {
        // Check if the point is inside the constraints
        if (is_point_inside_constraints(point, constraints))
        {
            cdt.insert(point); // εισαγωγή σημείου στην τριγωνοποίηση
            check_cdt_validity(cdt);
        }
        else
        {
            std::cerr << "Point (" << point.x() << ", " << point.y() << ") is outside constraints and will be skipped.\n";
        }
    }

    // επανάληψη για προσθήκη σημείων Steiner αν υπάρχουν αμβλυγώνια τρίγωνα
    bool all_acute = false;
    bool steiner_point_inserted;
    bool steiner_point_added_this_rotation;
    int no_of_steiner_points_added = 0;

    while (!all_acute && no_of_steiner_points_added < MAX_NO__STEINER_POINTS)
    {
        // Visualize the triangulation
        // export_to_svg(cdt, "output.svg");
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
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 0), constraints);
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
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 1), constraints);
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
                steiner_point_inserted = add_steiner_point_on_edge(cdt, CDT::Edge(face_it, 2), constraints);
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

double calculate_energy(const CDT& cdt, double alpha, double beta) {
    int obtuse_count = 0;
    int steiner_count = 0;

    // Count obtuse triangles
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        // Compute the squared edge lengths
        double a2 = CGAL::squared_distance(p2, p3);
        double b2 = CGAL::squared_distance(p1, p3);
        double c2 = CGAL::squared_distance(p1, p2);

        //ελεγχος για αμβλυα γωνια
        if (a2>b2+c2 || b2>a2+c2 || c2>a2+b2)
            obtuse_count++;

    }

    //υπολογισμος Steiner points
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
        if (vit->info().is_steiner) {
            steiner_count++;
        }
    }

    //τπος υπολογισμου ενεργειας απο δοαφανειες μαθηματος
    double energy = alpha * obtuse_count + beta * steiner_count;
    return energy;
}