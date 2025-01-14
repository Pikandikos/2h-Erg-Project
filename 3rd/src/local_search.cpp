#include "./func.h"
#include "insert_no_flip.h"
#include <iostream>

bool add_steiner_point_local_search(CDT &cdt, const CDT::Edge &edge, const std::vector<Point_2> &region_points)
{
    // Valid Steiner point candidates
    std::vector<Point_2> candidate_points;
    std::vector<contender> st_contenders;

    // Check if the edge is valid
    if (!edge.first->is_valid())
    {
        cerr << "Invalid edge detected, skipping Steiner point insertion" << endl;
        return false;
    }

    CDT::Vertex_handle vh1 = edge.first->vertex((edge.second + 1) % 3);
    CDT::Vertex_handle vh2 = edge.first->vertex((edge.second + 2) % 3);

    // Validate the face handle and vertices
    if (vh1 == nullptr || vh2 == nullptr)
    {
        std::cerr << "Invalid vertices detected, skipping Steiner point insertion." << endl;
        return false;
    }

    // Check for degeneracy before circumcenter calculation (εκφυλισμένη κορυφή)
    if (CGAL::collinear(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point()))
    {
        cerr << "Degenerate triangle detected, skipping Steiner point calculation" << endl;
        return false; // To avoid inserting into an invalid edge
    }

    // Attempt to use the circumcenter as the Steiner point
    Point_2 circumcenter = CGAL::circumcenter(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point());
    candidate_points.push_back(circumcenter);

    Point_2 midpoint = CGAL::midpoint(vh1->point(), vh2->point());
    candidate_points.push_back(midpoint);

    Point_2 projection = project_point_on_segment(edge.first->vertex(edge.second)->point(), CGAL::Segment_2<Kernel>(vh1->point(), vh2->point()));
    candidate_points.push_back(projection);

    // Centroid of the triangle
    Point_2 centroid = CGAL::centroid(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point());
    candidate_points.push_back(centroid);

    // Mean point of adjacent obtuse triangles
    Point_2 mean_point = mean_point_of_adjacent_triangles(cdt, edge.first, region_points);
    candidate_points.push_back(mean_point);

    for (int i = 0; i < candidate_points.size(); i++)
    {
        // Validate if the point is within constraints or already exists in the CDT
        if (point_exists_in_cdt(candidate_points[i], cdt))
        {
            std::cerr << "Candidate " << get_steiner_point_method(i) << " failed: duplicate point\n";
            continue;
        }

        if (!is_point_inside_constraints(candidate_points[i], region_points))
        {
            std::cerr << "Candidate " << get_steiner_point_method(i) << " failed: outside constraints\n";
            continue;
        }

        contender ct;
        ct.st_point = candidate_points[i];
        ct.method = get_steiner_point_method(i);
        ct.copy_cdt = cdt;

        // Before inserting, ensure the triangulation is in a valid state
        if (!cdt.is_valid())
        {
            std::cerr << "Triangulation is not valid before insertion!" << std::endl;
            return false;
        }

        // Check if the circumcenter is inside a neighboring triangle
        if (ct.method == "Circumcenter" && is_circumcenter_inside_neighbor_triangle(ct.copy_cdt, edge.first, ct.st_point))
        {
            // Remove the diagonal edge (if necessary) and insert Steiner point without flips
            // Before inserting, ensure the triangulation is in a valid state
            if (!remove_constraint_no_flip(ct.copy_cdt, edge))
            {
                std::cerr << "Failed to remove constrained edge!" << std::endl;
                continue;
            }
            insert_no_flip(ct.copy_cdt, candidate_points[i]);
            std::cout << "Circumcenter inside neighboring triangle, constraint removed, and Steiner point inserted." << std::endl;
        }
        ct.copy_cdt.insert(candidate_points[i]);
        export_to_svg(ct.copy_cdt, "output.svg");

        if (!check_cdt_validity(ct.copy_cdt))
        {
            std::cerr << "Invalid face detected during" << ct.method << " , skipping." << endl;
            continue; // Skip invalid faces
        }

        // Analyze obtuse angles in the new triangulation
        for (CDT::Finite_faces_iterator face_it = ct.copy_cdt.finite_faces_begin(); face_it != ct.copy_cdt.finite_faces_end(); ++face_it)
        {
            Point_2 p1 = face_it->vertex(0)->point();
            Point_2 p2 = face_it->vertex(1)->point();
            Point_2 p3 = face_it->vertex(2)->point();

            if (is_obtuse_triangle(p1, p2, p3))
            {
                ct.no_obtuse_faces++;
                double angle1 = angle_between_points(p1, p2, p3);
                double angle2 = angle_between_points(p2, p1, p3);
                double angle3 = angle_between_points(p3, p1, p2);

                // Sum of obtuse angles
                if (angle1 > 90)
                    ct.total_obtuse_angle_sum += angle1;
                if (angle2 > 90)
                    ct.total_obtuse_angle_sum += angle2;
                if (angle3 > 90)
                    ct.total_obtuse_angle_sum += angle3;

                // Track maximum angle
                ct.max_angle = std::max(ct.max_angle, std::max({angle1, angle2, angle3}));
            }
        }

        // Calculate current CDT's penalty score
        ct.cdt_penalty_score = (weight_obtuse_faces * ct.no_obtuse_faces) +
                               (weight_max_angle * ct.max_angle) +
                               (weight_total_obtuse_sum * ct.total_obtuse_angle_sum);

        st_contenders.push_back(ct);
    }

    // Compare the contenders based on custom metrics
    if (!st_contenders.empty())
    {
        int best_contender_index = 0;
        for (int i = 0; i < st_contenders.size(); i++)
        {
            if (st_contenders[best_contender_index].cdt_penalty_score > st_contenders[i].cdt_penalty_score)
                best_contender_index = i;
        }

        contender best_contender = st_contenders[best_contender_index];

        if (!is_point_inside_constraints(best_contender.st_point, region_points))
        {
            std::cerr << "Selected Steiner point failed constraints check, skipping insertion.\n";
            return false;
        }
        // cdt.insert(best_contender.st_point);

        cdt = best_contender.copy_cdt;

        cout << "Best Steiner point (" << best_contender.method << ") added at: ("
             << best_contender.st_point.x() << ", "
             << best_contender.st_point.y() << ") with penalty score: "
             << best_contender.cdt_penalty_score << "\n";
        return true;
    }

    std::cerr << "No Steiner points found that were valid, skipping insertion." << endl;
    return false;
}

void local_search(CDT &cdt, const std::vector<Point_2> &region_points, vector<double> method_parameters)
{
    cout << "\nStarting Local Search part of the Algorithm";

    int max_no_of_iterations = method_parameters[0]; // L

    // επανάληψη για προσθήκη σημείων Steiner αν υπάρχουν αμβλυγώνια τρίγωνα
    bool all_acute = false;
    bool steiner_point_inserted;
    int no_of_steiner_points_added = 0;

    // Sets up a check where it interrupts the loop if the no_obtuse_angles before
    // is smaller than after the addition of a steiner point
    int no_obtuse_angles_old = 100000;
    int no_obtuse_angles_curr = analyze_obtuse_angles(cdt);

    while (!all_acute && no_of_steiner_points_added < max_no_of_iterations)
    {
        export_to_svg(cdt, "output.svg");
        all_acute = true; // υποθετω ολα τα τριγωνα οξυγωνια
        cout << "\nStarting to check angles... again" << endl;

        cout << "Number of obtuse faces: " << analyze_obtuse_angles(cdt) << "/" << cdt.number_of_faces() << endl;
        no_obtuse_angles_old = analyze_obtuse_angles(cdt);

        for (CDT::Finite_faces_iterator face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); face_it++)
        {
            steiner_point_inserted = false;
            if (!face_it->is_valid())
            {
                std::cerr << "Invalid face detected, skipping." << endl;
                continue; // Skip invalid faces
            }

            Point_2 p1 = face_it->vertex(0)->point();
            Point_2 p2 = face_it->vertex(1)->point();
            Point_2 p3 = face_it->vertex(2)->point();

            cout << "P1:" << p1 << "  P2:" << p2 << "  P3:" << p3 << endl;

            double angles[3];
            // Υπολογισμός των γωνιών του τριγώνου
            angles[0] = angle_between_points(p1, p2, p3);
            angles[1] = angle_between_points(p2, p1, p3);
            angles[2] = angle_between_points(p3, p1, p2);

            cout << "Angle1:" << angles[0] << "  Angle2:" << angles[1] << "  Angle3:" << angles[2] << endl;

            if (angles[0] == 0 || angles[1] == 0 || angles[2] == 0)
            {
                cout << "Angle = 0 Found !!!!" << endl;
                continue;
            }

            int obtuse_index = -1;
            for (int i = 0; i < 3; ++i)
            {
                if (angles[i] > 90)
                {
                    obtuse_index = i;
                    break;
                }
            }

            // An obtuse angle exists
            if (obtuse_index != -1)
            {
                all_acute = false;
                bool flipped = false;
                steiner_point_inserted = add_steiner_point_local_search(cdt, CDT::Edge(face_it, obtuse_index), region_points);
                if (steiner_point_inserted) // steiner point hasn't been skipped
                {
                    no_of_steiner_points_added++;
                    cout << "No. of Steiner Points: " << no_of_steiner_points_added << endl;
                    face_it = cdt.finite_faces_begin();
                    break; // Steiner Point has been added... Now Restart the checking of angles
                }
            }
        }

        // In case triangles start increasing meaning algorithm starts getting less effective
        no_obtuse_angles_curr = analyze_obtuse_angles(cdt);
        if (no_obtuse_angles_old < no_obtuse_angles_curr)
        {
            cout << "Obtuse faces/triangles increasing... loop breaks" << endl;
            break;
        }
        no_obtuse_angles_old = no_obtuse_angles_curr;

        if (steiner_point_inserted == false && all_acute == false) // Will need to add flip if it was working
        {
            cout << "All faces/triangles are acute... loop finishes" << endl;
            all_acute = true; // If the loop ever ends it means that all faces are acute
        }
    }
}