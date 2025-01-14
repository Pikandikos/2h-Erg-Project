#include "./func.h"

// Calculate the radius-to-height ratio of a triangle
double calculate_radius_to_height_ratio(const Point_2 &p1, const Point_2 &p2, const Point_2 &p3)
{
    if (!is_obtuse_triangle(p1, p2, p3))
    {
        return 0.0; // Face Not obtuse
    }

    // Calculate circumradius (R)
    double a = CGAL::squared_distance(p1, p2);
    double b = CGAL::squared_distance(p2, p3);
    double c = CGAL::squared_distance(p3, p1);

    double R = std::sqrt((a * b * c) / (a + b + c)); // Simplified formula for circumradius

    // Calculate height from the longest side
    double height = std::sqrt(CGAL::to_double(CGAL::area(p1, p2, p3)) * 4 / R);

    return R / height;
}

// Calculate the heuristic value for a triangle
double calculate_heuristic(const CDT::Face_handle &face)
{
    Point_2 p1 = face->vertex(0)->point();
    Point_2 p2 = face->vertex(1)->point();
    Point_2 p3 = face->vertex(2)->point();

    double rho = calculate_radius_to_height_ratio(p1, p2, p3);

    if (rho > 2.0)
    {
        return std::max(0.0, (rho - 1) / rho); // Vertex projection
    }
    else if (rho >= 1.0 && rho <= 2.0)
    {
        return rho / (2 + rho); // Circumcenter
    }
    else if (rho < 1.0)
    {
        return std::max(0.0, (3 - 2 * rho) / 3); // Midpoint of the longest edge
    }
    else
    {
        return 0.0; // Default fallback
    }
}

Point_2 select_next_point_ant(const CDT &cdt, const std::vector<Point_2> &region_points, const CDT::Edge &edge, std::map<Point_2, double> &pheromone, double xi, double psi)
{
    CDT::Vertex_handle vh1 = edge.first->vertex((edge.second + 1) % 3);
    CDT::Vertex_handle vh2 = edge.first->vertex((edge.second + 2) % 3);

    // Valid Steiner point candidates
    std::vector<Point_2> candidate_points;

    // Calculate Steiner points
    Point_2 circumcenter = CGAL::circumcenter(vh1->point(), vh2->point(), edge.first->vertex(edge.second)->point());
    candidate_points.push_back(circumcenter);

    Point_2 midpoint = CGAL::midpoint(vh1->point(), vh2->point());
    candidate_points.push_back(midpoint);

    Point_2 projection = project_point_on_segment(edge.first->vertex(edge.second)->point(), CGAL::Segment_2<Kernel>(vh1->point(), vh2->point()));
    candidate_points.push_back(projection);

    Point_2 mean_point = mean_point_of_adjacent_triangles(cdt, edge.first, region_points);
    candidate_points.push_back(mean_point);

    std::vector<double> probabilities(candidate_points.size());
    double total_probability = 0.0;

    // Calculate heuristic values and probabilities for each candidate
    for (int i = 0; i < candidate_points.size(); ++i)
    {
        Point_2 candidate = candidate_points[i];

        if (!point_exists_in_cdt(candidate, cdt))
        {
            CDT::Face_handle face = cdt.locate(candidate);
            if (cdt.is_infinite(face))
            {
                probabilities[i] = 0.0;
                continue;
            }

            double heuristic = calculate_heuristic(face);

            // Calculate the probability of this candidate based on the heuristic and pheromone
            double pheromone_value = pheromone[candidate];
            double probability = std::pow(pheromone_value, xi) * std::pow(heuristic, psi);
            probabilities[i] = probability;
            total_probability += probability;
        }
    }

    // Normalize probabilities
    for (double &probability : probabilities)
    {
        probability /= total_probability;
    }

    // Select a point based on the calculated probabilities
    double random_value = (double)rand() / RAND_MAX;
    double cumulative_probability = 0.0;

    for (int i = 0; i < probabilities.size(); ++i)
    {
        cumulative_probability += probabilities[i];
        if (random_value <= cumulative_probability)
        {
            cout << "Steiner Point to be tested: " << candidate_points[i] << endl;
            return candidate_points[i];
        }
    }

    cout << "Steiner Point to be tested: " << candidate_points[0] << endl;
    return candidate_points[0]; // Default fallback
}

Point_2 select_next_point_with_dynamic_pheromone(CDT &local_triangulation, const std::vector<Point_2> &region_points,
                                                 std::map<Point_2, double> &pheromone, double xi, double psi)
{
    // Weighing by pheromone levels: more pheromone means higher likelihood of selection
    std::vector<CDT::Finite_faces_iterator> obtuse_faces;

    // Collect obtuse faces from the triangulation
    for (CDT::Finite_faces_iterator face_it = local_triangulation.finite_faces_begin();
         face_it != local_triangulation.finite_faces_end(); ++face_it)
    {
        Point_2 p1 = face_it->vertex(0)->point();
        Point_2 p2 = face_it->vertex(1)->point();
        Point_2 p3 = face_it->vertex(2)->point();

        if (!local_triangulation.is_infinite(face_it) && is_obtuse_triangle(p1, p2, p3))
        {
            obtuse_faces.push_back(face_it);
        }
    }

    // Calculate pheromone-based weights for each obtuse face
    std::vector<double> weights;
    double total_weight = 0.0;

    for (auto &face_it : obtuse_faces)
    {
        Point_2 p1 = face_it->vertex(0)->point();
        Point_2 p2 = face_it->vertex(1)->point();
        Point_2 p3 = face_it->vertex(2)->point();

        // Calculate the total pheromone influence for the current obtuse face (sum of pheromone on vertices)
        double pheromone_weight = pheromone[p1] + pheromone[p2] + pheromone[p3];
        weights.push_back(pheromone_weight);
        total_weight += pheromone_weight;
    }

    // Calculate the probabilistic selection for each obtuse face based on pheromone trail
    std::vector<double> probabilities(obtuse_faces.size(), 0.0);

    for (int i = 0; i < obtuse_faces.size(); ++i)
    {
        probabilities[i] = weights[i] / total_weight; // Normalize to get a probability distribution
    }

    // Generate a random number and select a face based on the probability distribution
    double random_value = rand() / (double)RAND_MAX;
    double cumulative_prob = 0.0;
    int selected_face_idx = 0;

    for (int i = 0; i < obtuse_faces.size(); ++i)
    {
        cumulative_prob += probabilities[i];
        if (random_value <= cumulative_prob)
        {
            selected_face_idx = i;
            break;
        }
    }

    // Retrieve the selected obtuse triangle and calculate the next Steiner point
    CDT::Finite_faces_iterator selected_face = obtuse_faces[selected_face_idx];
    return select_next_point_ant(local_triangulation, region_points, CDT::Edge(selected_face, 0), pheromone, xi, psi);
}

// Main Ant Colony Optimization
bool ant_colony_optimization(CDT &cdt, const std::vector<Point_2> &region_points, std::vector<double> &method_parameters)
{
    double alpha = method_parameters[0]; // a,b weights for energy calculation
    double beta = method_parameters[1];
    double xi = method_parameters[2]; // xi, psi for calculating the possibilities of selecting a point
    double psi = method_parameters[3];
    double lambda = method_parameters[4]; // lambda for evaporation rate
    int kappa = method_parameters[5];     // No of ants
    int L = method_parameters[6];         // No of iterations

    // Initialize pheromone trails
    std::map<Point_2, double> pheromone;
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        pheromone[vit->point()] = 1.0; // Initialize pheromone trail
    }

    CDT best_triangulation = cdt;
    int st_aco_counter = 0;

    // Optimization cycles
    for (int i = 0; i < L; ++i)
    {
        CDT local_triangulation = best_triangulation;
        std::vector<Point_2> ants_selected_points;

        for (int ant = 0; ant < kappa; ++ant)
        {
            CDT::Finite_faces_iterator face_it = local_triangulation.finite_faces_begin();
            while (face_it != local_triangulation.finite_faces_end())
            {
                Point_2 p1 = face_it->vertex(0)->point();
                Point_2 p2 = face_it->vertex(1)->point();
                Point_2 p3 = face_it->vertex(2)->point();

                if (!local_triangulation.is_infinite(face_it) && is_obtuse_triangle(p1, p2, p3))
                {
                    // Use dynamic pheromone-based selection
                    Point_2 new_steiner_point = select_next_point_with_dynamic_pheromone(local_triangulation, region_points, pheromone, xi, psi);

                    // Check if the point already exists in the triangulation
                    if (!point_exists_in_cdt(new_steiner_point, local_triangulation))
                    {
                        ants_selected_points.push_back(new_steiner_point);
                        break; // Stop after selecting one valid Steiner point
                    }
                }
                ++face_it;
            }
        }

        // After the cycle, insert the selected Steiner points into the triangulation
        for (const auto &steiner_point : ants_selected_points)
        {
            local_triangulation.insert(steiner_point);
        }

        // Calculate the energy of the triangulation after inserting the selected points
        double energy = calculate_energy(local_triangulation, alpha, beta);
        cout << "Energy calculated: " << energy << endl;

        // Update the best triangulation if the current one has lower energy
        if (energy < calculate_energy(best_triangulation, alpha, beta))
        {
            best_triangulation = local_triangulation;
            cout << "Steiner Point(s) added using ACO: ";
            for (const auto &point : ants_selected_points)
            {
                cout << point << " ";
            }
            cout << endl;
            st_aco_counter++;
        }

        // Update pheromones based on the triangulation's faces
        for (auto &entry : pheromone)
        {
            entry.second = (1.0 - lambda) * entry.second; // Evaporation
        }

        for (CDT::Finite_faces_iterator face_it = best_triangulation.finite_faces_begin();
             face_it != best_triangulation.finite_faces_end(); ++face_it)
        {
            Point_2 p1 = face_it->vertex(0)->point();
            Point_2 p2 = face_it->vertex(1)->point();
            Point_2 p3 = face_it->vertex(2)->point();
            if (!best_triangulation.is_infinite(face_it) && is_obtuse_triangle(p1, p2, p3))
            {
                // Increase pheromone level for the vertices of obtuse triangles
                pheromone[face_it->vertex(0)->point()] += 1.0;
                pheromone[face_it->vertex(1)->point()] += 1.0;
                pheromone[face_it->vertex(2)->point()] += 1.0;
            }
        }
    }

    cout << "Number of Steiner Points added using ACO: " << st_aco_counter << endl;

    // Final best triangulation after ACO optimization
    cdt = best_triangulation;
    return true;
}
