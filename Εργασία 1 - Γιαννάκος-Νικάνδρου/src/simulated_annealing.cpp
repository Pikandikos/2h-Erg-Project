#include "./func.h"
#include <iostream>

// double calculate_energy(const CDT &cdt, double alpha, double beta)
// {
//     int obtuse_count = 0;
//     int steiner_count = 0;

//     // Count obtuse triangles
//     for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
//     {
//         Point_2 p1 = fit->vertex(0)->point();
//         Point_2 p2 = fit->vertex(1)->point();
//         Point_2 p3 = fit->vertex(2)->point();

//         // Compute the squared edge lengths
//         double a2 = CGAL::squared_distance(p2, p3);
//         double b2 = CGAL::squared_distance(p1, p3);
//         double c2 = CGAL::squared_distance(p1, p2);

//         // ελεγχος για αμβλυα γωνια
//         if (a2 > b2 + c2 || b2 > a2 + c2 || c2 > a2 + b2)
//             obtuse_count++;
//     }

//     // υπολογισμος Steiner points
//     for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
//     {
//         if (vit->info().is_steiner)
//         {
//             steiner_count++;
//         }
//     }

//     // τπος υπολογισμου ενεργειας απο δοαφανειες μαθηματος
//     double energy = alpha * obtuse_count + beta * steiner_count;
//     return energy;
// }

// // προσομοιωμενη ανόπτηση
// void simulated_annealing(Triangulation &triangulation, double alpha, double beta, int max_iterations)
// {
//     double T = 1.0;         // αρχικη θερμοκρασια
//     int L = max_iterations; // αριθμος επαναληψεων

//     // Define random number generation for Steiner options
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_int_distribution<int> dist(0, 4); // 5 Steiner point options

//     while (T > 0)
//     {
//         for (auto &triangle : triangulation.triangles)
//         {
//             // Check if the triangle is obtuse dynamically
//             auto vertices = triangle.get_vertices();
//             Point a = vertices[0], b = vertices[1], c = vertices[2];
//             double a2 = CGAL::squared_distance(b, c);
//             double b2 = CGAL::squared_distance(a, c);
//             double c2 = CGAL::squared_distance(a, b);

//             // Determine if the triangle is obtuse
//             if (a2 > b2 + c2 || b2 > a2 + c2 || c2 > a2 + b2)
//             {
//                 // Generate a random Steiner point location (pseudo-code)
//                 Point steiner_option = {/* Random point */};

//                 // Store the current triangulation as a backup
//                 Triangulation backup = triangulation;

//                 // Store the current energy before insertion
//                 double initial_energy = calculate_energy(triangulation, alpha, beta);

//                 // Insert Steiner point
//                 // add_steiner_point_local_search(triangulation, steiner_option);
//                 bool insertion_succesful = add_steiner_point_local_search(triangulation, steiner_option);

//                 // Calculate new energy
//                 if (insertion_succesful)
//                 {
//                     double new_energy = calculate_energy(triangulation, alpha, beta);
//                     double delta_energy = new_energy - initial_energy;

//                     // Accept or reject the new configuration
//                     if (delta_energy < 0 || exp(-delta_energy / T) > ((double)rand() / RAND_MAX))
//                     {
//                         // Accept new configuration
//                         continue;
//                     }
//                     else
//                     {
//                         // Reject and restore the previous state
//                         triangulation = backup;
//                     }
//                 }
//             }
//         }
//         // Decrease temperature
//         T -= 1.0 / L;
//     }
// }