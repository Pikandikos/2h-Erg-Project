// Find the Circumcenter
Point_2 circumcenter = CGAL::circumcenter(vh1->point(), vh2->point());

// Check if circumcenter exists
if (cdt.locate(circumcenter) != nullptr)
{
    cerr << "Circumcenter already exists, skipping insertion." << endl;
    cerr << "Skipping insertion for edge with vertices: ("
         << vh1->point() << ") and (" << vh2->point() << ") ";
    cerr << "Triangle vertices: (" << edge.first->vertex(0)->point() << "), "
         << "(" << edge.first->vertex(1)->point() << "), "
         << "(" << edge.first->vertex(2)->point() << ")" << endl;
}
else
{
    CDT::Vertex_handle new_vertex = cdt.insert(circumcenter);
    if (new_vertex != nullptr)
    {
        cout << "Steiner Point (circumcenter) added: (" << circumcenter.x() << ", " << circumcenter.y() << ")" << endl;
        return true;
    }
}

// Find the Midpoint
Point_2 midpoint = CGAL::midpoint(vh1->point(), vh2->point());

// Check if midpoint exists
if (cdt.locate(midpoint) != nullptr)
{
    cerr << "Midpoint already exists, skipping insertion." << endl;
    cerr << "Skipping insertion for edge with vertices: ("
         << vh1->point() << ") and (" << vh2->point() << ") ";
    cerr << "Triangle vertices: (" << edge.first->vertex(0)->point() << "), "
         << "(" << edge.first->vertex(1)->point() << "), "
         << "(" << edge.first->vertex(2)->point() << ")" << endl;
}
else
{
    CDT::Vertex_handle new_vertex = cdt.insert(midpoint);
    if (new_vertex != nullptr)
    {
        cout << "Steiner Point (midpoint) added: (" << midpoint.x() << ", " << midpoint.y() << ")" << endl;
        return true;
    }
}

if (cdt.locate(circumcenter) != nullptr && cdt.locate(midpoint) != nullptr)
{
    Point_2 perturbed_point = Point_2(circumcenter.x() + 0.01, circumcenter.y() + 0.01);
    CDT::Vertex_handle new_vertex = cdt.insert(perturbed_point);
    if (new_vertex != nullptr)
    {
        cout << "Steiner Point (perturbed circumcenter) added: (" << perturbed_point.x() << ", " << perturbed_point.y() << ")" << endl;
        return true;
    }
    cerr << "Perturbed point could not be added, skipping." << endl;
}
// ------------------------------------------------------------------------------------------------------- //

4. Handle Unresolvable Cases
    If acute angles cannot be achieved after a certain number of iterations,
    output a warning and gracefully terminate : int max_steiner_points = 100; // Set a limit to avoid infinite loops
if (no_of_steiner_points_added >= max_steiner_points)
{
    cerr << "Warning: Maximum Steiner points added. Acute triangulation may not be possible." << endl;
    break;
}

5. Final Validation and Debugging
    Add a function to validate the triangulation after the process : void
                                                                     validate_acute_triangulation(const CDT &cdt)
{
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        if (!fit->is_valid())
            continue;

        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        double angle1 = angle_between_points(p1, p2, p3);
        double angle2 = angle_between_points(p2, p1, p3);
        double angle3 = angle_between_points(p3, p1, p2);

        if (angle1 >= 90 || angle2 >= 90 || angle3 >= 90)
        {
            cerr << "Non-acute triangle detected: " << endl;
            cerr << "Vertices: (" << p1 << "), (" << p2 << "), (" << p3 << ")" << endl;
            cerr << "Angles: " << angle1 << ", " << angle2 << ", " << angle3 << endl;
        }
    }
}

// ------------------------------------------------------------------------------------------------------- //

cout << "Edge vertices: (" << vh1->point() << "), (" << vh2->point() << ")" << endl;
cout << "Triangle vertices: ("
     << edge.first->vertex(0)->point() << "), "
     << edge.first->vertex(1)->point() << "), "
     << edge.first->vertex(2)->point() << ")" << endl;

// ------------------------------------------------------------------------------------------------------- //

bool add_steiner_point_on_edge(CDT &cdt, const CDT::Edge &edge, const vector<pair<Point_2, Point_2>> &constraints)
{
    vector<Point_2> candidates;

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