#include "./func.h"
#include <iostream>

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

bool is_obtuse_triangle(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    double angle1 = angle_between_points(p1, p2, p3);
    double angle2 = angle_between_points(p2, p1, p3);
    double angle3 = angle_between_points(p3, p1, p2);
    return (angle1 > 90 || angle2 > 90 || angle3 > 90);
}

bool check_cdt_validity(const CDT &cdt)
{
    for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        if (!fit->is_valid())
        {
            std::cerr << "Invalid face detected!" << endl;
            return false;
        }
    }
    // cout << "All Faces valid for now..." << endl;
    return true;
}

bool is_point_inside_constraints(const Point_2 &point, const std::vector<Point_2> &region_points)
{
    Polygon_2 region_polygon(region_points.begin(), region_points.end());

    // Check the bounded side
    CGAL::Bounded_side side = region_polygon.bounded_side(point);

    // ON_BOUNDED_SIDE: Point is inside the polygon
    // ON_BOUNDARY: The point is exactly on one of the edges or vertices of the polygon
    return (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY);
}

void insert_points_into_cdt(vector<Point_2> &points, vector<int> &region_boundary, const vector<pair<int, int>> &additional_constraints, CDT &cdt)
{
    // Vector to store points in the region boundary
    vector<Point_2> region_points;

    // Insert all points into the CDT
    for (const auto &point : points)
    {
        cdt.insert(point);
    }
    cout << "Points inserted into CDT" << endl;

    // Add region boundary as constraints
    cout << "Adding region boundary constraints..." << endl;
    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int start_idx = region_boundary[i];
        // Wrap around for the last edge
        int end_idx = region_boundary[(i + 1) % region_boundary.size()];
        cdt.insert_constraint(points[start_idx], points[end_idx]);
    }

    // Add additional constraints
    cout << "Adding additional constraints..." << endl;
    for (const auto &constraint : additional_constraints)
    {
        int start_idx = constraint.first;
        int end_idx = constraint.second;
        cdt.insert_constraint(points[start_idx], points[end_idx]);
    }

    cout << "Constraints added successfully" << endl;
}

vector<Point_2> create_region_points(vector<Point_2> &points, vector<int> &region_boundary)
{
    // Vector to store points in the region boundary
    vector<Point_2> region_points;

    for (size_t i = 0; i < region_boundary.size(); ++i)
    {
        int idx = region_boundary[i];
        // Ensure unique points are added to the region_points vector
        if (region_points.empty() || region_points.back() != points[idx])
        {
            region_points.push_back(points[idx]);
        }
    }
    return region_points;
}

// Function to find and print obtuse angles in the CDT
int analyze_obtuse_angles(const CDT &cdt)
{
    int obtuse_count = 0; // Counter for obtuse angles

    for (CDT::Finite_faces_iterator face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); ++face_it)
    {
        // Get the points of the triangle
        Point_2 p1 = face_it->vertex(0)->point();
        Point_2 p2 = face_it->vertex(1)->point();
        Point_2 p3 = face_it->vertex(2)->point();

        // Calculate angles
        double angle1 = angle_between_points(p1, p2, p3);
        double angle2 = angle_between_points(p2, p1, p3);
        double angle3 = angle_between_points(p3, p1, p2);

        // Check for obtuse angles
        if (angle1 > 90 || angle2 > 90 || angle3 > 90)
        {
            obtuse_count++;
            // cout << "Obtuse triangle found:\n";
            // cout << "  Vertices: (" << p1 << "), (" << p2 << "), (" << p3 << ")\n";
            // cout << "  Angles: " << angle1 << "°, " << angle2 << "°, " << angle3 << "°\n";
        }
    }
    return obtuse_count;
}

Point_2 mean_point_of_adjacent_triangles(CDT &cdt, CDT::Face_handle face, const std::vector<Point_2> &region_points)
{
    std::vector<Point_2> centroids;

    Point_2 mean(0, 0);
    int count = 0;

    for (int i = 0; i < 3; ++i)
    {
        CDT::Face_handle neighbour = face->neighbor(i);
        if (cdt.is_infinite(neighbour))
            continue;

        Point_2 p1 = neighbour->vertex(0)->point();
        Point_2 p2 = neighbour->vertex(1)->point();
        Point_2 p3 = neighbour->vertex(2)->point();

        if (is_obtuse_triangle(p1, p2, p3))
        {
            centroids.push_back(CGAL::centroid(p1, p2, p3));
        }
    }

    if (!centroids.empty())
    {
        // Compute the centroid of the collected points
        return CGAL::centroid(centroids.begin(), centroids.end());
    }
    else
    {
        // Return the midpoint of the edge as a fallback
        return CGAL::midpoint(face->vertex(0)->point(), face->vertex(1)->point());
    }
}

Point_2 project_point_on_segment(const Point_2 &p, const Segment_2 &s)
{
    CGAL::Line_2<Kernel> line(s.source(), s.target());
    return line.projection(p);
}

// Function to check if a point already exists in the CDT
bool point_exists_in_cdt(const Point_2 &point, const CDT &cdt)
{
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->point() == point)
        {
            return true; // Point exists
        }
    }
    return false; // Point does not exist
}

string get_steiner_point_method(int i)
{
    switch (i)
    {
    case 0:
        return "Circumcenter";
    case 1:
        return "Midpoint";
    case 2:
        return "Projection";
    case 3:
        return "Centroid";
    case 4:
        return "Mean Point";
    default:
        return "Unknown Method"; // Fallback for invalid indices
    }
}

void populateVector(const std::string &method, const boost::property_tree::ptree &parameters, std::vector<double> &params)
{
    params.clear(); // Clear the vector before populating

    std::cout << "Parameters for " << method << ":" << std::endl;

    // Iterate over the parameters and store the values in the vector
    for (const auto &param : parameters)
    {
        try
        {
            // Convert the value to double and add to the vector
            params.push_back(std::stod(param.second.data()));
        }
        catch (const std::exception &e)
        {
            std::cerr << "Error converting parameter \"" << param.first
                      << "\" to double: " << e.what() << std::endl;
        }
    }
}

bool is_circumcenter_inside_neighbor_triangle(const CDT &cdt, const CDT::Face_handle &neighbor_face, const Point_2 &circumcenter)
{
    // Collect the vertices of the triangle (neighboring face)
    std::vector<Point_2> triangle_points;
    triangle_points.push_back(neighbor_face->vertex(0)->point());
    triangle_points.push_back(neighbor_face->vertex(1)->point());
    triangle_points.push_back(neighbor_face->vertex(2)->point());

    // Ensure correct usage of bounded_side_2, use the proper point type
    CGAL::Bounded_side neighbor_side = CGAL::bounded_side_2(triangle_points.begin(), triangle_points.end(), circumcenter, Kernel());

    // If the circumcenter is inside the triangle, the result will be ON_BOUNDED_SIDE
    return neighbor_side == CGAL::ON_BOUNDED_SIDE;
}