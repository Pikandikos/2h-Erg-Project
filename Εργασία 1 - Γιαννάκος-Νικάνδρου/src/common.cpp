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

bool is_point_inside_constraints(const Point_2 &point, const std::vector<std::pair<Point_2, Point_2>> &constraints)
{
    // Create a set of unique points from the constraints (ensures no duplicate points)
    std::set<Point_2> unique_points;
    for (const auto &constraint : constraints)
    {
        unique_points.insert(constraint.first);
        unique_points.insert(constraint.second);
    }

    // If there are not enough points to form a polygon, return false (degenerate case)
    if (unique_points.size() < 3)
    {
        std::cerr << "Constraints form a degenerate polygon (less than 3 unique points)." << std::endl;
        return false;
    }

    // Create a vector of the unique points to form the polygon
    std::vector<Point_2> polygon_points(unique_points.begin(), unique_points.end());

    // Check for degenerate cases where points might be collinear
    // We will check if the points form a valid simple polygon
    Polygon_2 constraint_polygon(polygon_points.begin(), polygon_points.end());

    // Ensure the polygon is simple (non-self-intersecting)
    if (!constraint_polygon.is_simple())
    {
        std::cerr << "Constraints form a non-simple polygon!" << std::endl;
        return false;
    }

    // Check if the point lies inside the polygon
    return constraint_polygon.bounded_side(point) != CGAL::ON_UNBOUNDED_SIDE;
}

#include <tuple>

// Function to find and print obtuse angles in the CDT
void analyze_obtuse_angles(const CDT &cdt)
{
    int obtuse_count = 0; // Counter for obtuse angles

    cout << "Analyzing triangles for obtuse angles...\n";

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
    cout << "Total obtuse angles found: " << obtuse_count << "\n";
}

Point_2 mean_point_of_adjacent_triangles(CDT &cdt, CDT::Face_handle face, const vector<pair<Point_2, Point_2>> &constraints)
{
    Point_2 mean(0, 0);
    int count = 0;

    for (int i = 0; i < 3; ++i)
    {
        CDT::Face_handle neighbor = face->neighbor(i);
        if (cdt.is_infinite(neighbor))
            continue;

        Point_2 p1 = neighbor->vertex(0)->point();
        Point_2 p2 = neighbor->vertex(1)->point();
        Point_2 p3 = neighbor->vertex(2)->point();

        if (is_obtuse_triangle(p1, p2, p3))
        {
            mean = mean + CGAL::centroid(p1, p2, p3);
            count++;
        }
    }

    return (count > 0) ? mean / count : face->vertex(0)->point(); // Fallback
}
