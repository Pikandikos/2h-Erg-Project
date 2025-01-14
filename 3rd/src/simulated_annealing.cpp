#include "./func.h"
#include <iostream>

double calculate_energy(const CDT &cdt, double alpha, double beta)
{
    int obtuse_count = 0;
    int steiner_count = 0;

    // Count obtuse triangles
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
    {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        // Compute the squared edge lengths
        double a2 = CGAL::squared_distance(p2, p3);
        double b2 = CGAL::squared_distance(p1, p3);
        double c2 = CGAL::squared_distance(p1, p2);

        // Check for obtuse angle
        if (a2 > b2 + c2 || b2 > a2 + c2 || c2 > a2 + b2)
            obtuse_count++;
    }

    // Calculate Steiner points
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (point_exists_in_cdt(vit->point(), cdt))
            steiner_count++;
    }

    // Calculate energy
    double energy = alpha * obtuse_count + beta * steiner_count;
    return energy;
}

// Simulated Annealing
bool simulated_annealing(CDT &cdt, const std::vector<Point_2> &region_points, std::vector<double> method_parameters)
{
    cout << "\nStarting Simulated Annealing part of the Algorithm";

    double alpha = method_parameters[0];
    double beta = method_parameters[1];
    double T = 1.0;               // Iterations per temperature level
    int L = method_parameters[2]; // Number of iterations
    int total_iterations = 0;     // Track the number of iterations

    // Random generator setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    while (T > 1e-5) // Stop when temperature is close to zero
    {
        for (auto face_it = cdt.finite_faces_begin(); face_it != cdt.finite_faces_end(); ++face_it)
        {
            if (!face_it->is_valid())
            {
                std::cerr << "Invalid face detected, skipping." << std::endl;
                continue;
            }

            // Triangle vertices
            Point_2 a = face_it->vertex(0)->point();
            Point_2 b = face_it->vertex(1)->point();
            Point_2 c = face_it->vertex(2)->point();

            // Compute edge lengths squared
            double a2 = CGAL::squared_distance(b, c);
            double b2 = CGAL::squared_distance(a, c);
            double c2 = CGAL::squared_distance(a, b);

            // Compute angles of the triangle
            double angles[3];
            angles[0] = angle_between_points(a, b, c);
            angles[1] = angle_between_points(b, a, c);
            angles[2] = angle_between_points(c, a, b);

            // Find the obtuse angle
            int obtuse_index = -1;
            for (int i = 0; i < 3; ++i)
            {
                if (angles[i] > 90)
                {
                    obtuse_index = i;
                    break;
                }
            }

            // Process obtuse triangles only
            if (obtuse_index != -1)
            {
                // Backup the current CDT
                CDT backup = cdt;

                // Calculate initial energy
                double initial_energy = calculate_energy(cdt, alpha, beta);

                if (add_steiner_point_local_search(cdt, CDT::Edge(face_it, obtuse_index), region_points))
                {
                    // Calculate new energy
                    double new_energy = calculate_energy(cdt, alpha, beta);
                    double delta_energy = new_energy - initial_energy;

                    // Accept or reject the new configuration
                    if (!(delta_energy < 0 || exp(-delta_energy / T) > dist(gen)))
                    {
                        // Reject and restore backup
                        cdt = backup;
                        cout << "Reject and restore backup\n";
                    }
                }
            }

            total_iterations++;
            if (total_iterations >= L)
                break;
        }
        // Exponential cooling
        T *= 0.95;
    }
    return true;
}
