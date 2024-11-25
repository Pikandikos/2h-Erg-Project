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

// Point_2 compute_steiner_point_from_neighbors(CDT &cdt, CDT::Face_iterator face_it)
// {
//     std::vector<Point_2> neighbors_polygon;

//     // Iterate over neighboring faces and collect their vertices
//     for (auto edge : cdt.finite_edges_of(face_it))
//     {
//         CDT::Face_iterator neighbor_face = edge.second; // Get neighboring face
//         if (neighbor_face != face_it)
//         {
//             // Collect vertices that form the boundary of the polygon
//             Point_2 neighbor_vertex = edge.first->get_opposite_vertex(face_it);
//             neighbors_polygon.push_back(neighbor_vertex);
//         }
//     }

//     // Find a suitable Steiner point inside the polygon (e.g., centroid)
//     Point_2 steiner_point = compute_centroid(neighbors_polygon);

//     return steiner_point;
// }

// bool add_steiner_point_from_neighbors(CDT &cdt, CDT::Face_iterator face_it)
// {
//     // Get a Steiner point using the neighboring faces
//     Point_2 steiner_point = compute_steiner_point_from_neighbors(cdt, face_it);

//     // Check if the point is valid (e.g., inside the convex hull of neighbors)
//     if (is_valid_steiner_point(cdt, steiner_point))
//     {
//         // Insert the Steiner point into the triangulation
//         cdt.insert(steiner_point);
//         return true;
//     }
//     return false;
// }
