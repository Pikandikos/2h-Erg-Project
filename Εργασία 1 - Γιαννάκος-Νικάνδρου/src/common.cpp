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