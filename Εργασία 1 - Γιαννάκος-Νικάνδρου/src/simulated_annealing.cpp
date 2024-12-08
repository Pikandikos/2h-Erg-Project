#include "./func.h"
#include <iostream>

double calculate_energy(const CDT &cdt, std::vector<double> method_parameters)
{
    double alpha = method_parameters[0];
    double beta = method_parameters[1];
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

        // ελεγχος για αμβλυα γωνια
        if (a2 > b2 + c2 || b2 > a2 + c2 || c2 > a2 + b2)
            obtuse_count++;
    }

    // υπολογισμος Steiner points
    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        if (vit->info().is_steiner)
        {
            steiner_count++;
        }
    }

    // τπος υπολογισμου ενεργειας απο δοαφανειες μαθηματος
    double energy = alpha * obtuse_count + beta * steiner_count;
    return energy;
}

Point_2 generate_random_steiner_point(const Point_2 &a, const Point_2 &b, const Point_2 &c, std::mt19937 &gen, std::uniform_int_distribution<int> &dist)
{
    // δημιουργια συντεταγμενων
    double r1 = dist(gen) / 4.0; // τυχαια τιμη [0, 1]
    double r2 = dist(gen) / 4.0;
    if (r1 + r2 > 1.0)
    {
        r1 = 1.0 - r1;
        r2 = 1.0 - r2;
    }
    double r3 = 1.0 - r1 - r2;

    // τυχαιο σημειο μεσα στο τριγωνο
    return Point_2(r1 * a.x() + r2 * b.x() + r3 * c.x(),
                   r1 * a.y() + r2 * b.y() + r3 * c.y());
}

// προσομοιωμενη ανόπτηση
void simulated_annealing(const CDT &cdt, std::vector<double> method_parameters)
{
    double alpha = method_parameters[0];
    double beta = method_parameters[1];
    double T = 1.0;               // αρχικη θερμοκρασια
    int L = method_parameters[2]; // αριθμος επαναληψεων

    // για τα steiner
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, 4); // 5 επιλογες σημεια steiner

    while (T > 0)
    {
        for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
        {
            // κορυφες τριγωνου
            Point_2 a = fit->vertex(0)->point();
            Point_2 b = fit->vertex(1)->point();
            Point_2 c = fit->vertex(2)->point();

            // υπολογισμως μηκων (τετραγωνικων) των ακρων
            double a2 = CGAL::squared_distance(b, c);
            double b2 = CGAL::squared_distance(a, c);
            double c2 = CGAL::squared_distance(a, b);

            // αμβλυγωνιο τριγωνο
            if (a2 > b2 + c2 || b2 > a2 + c2 || c2 > a2 + b2)
            {
                // τυχαιο σημειο steiner
                Point_2 steiner_option = generate_random_steiner_point(a, b, c, gen, dist);

                // backup τριγωνοποιηση
                CDT backup = cdt;

                // υπολογισμος αρχικς ενεργειας
                double initial_energy = calculate_energy(cdt, method_parameters);

                // προσπαθεια πορσθηκης Steiner point
                bool insertion_successful = add_steiner_point(cdt, steiner_option);

                if (insertion_successful)
                {
                    // Calculate new energy
                    double new_energy = calculate_energy(cdt, method_parameters);
                    double delta_energy = new_energy - initial_energy;

                    // Accept or reject the new configuration
                    if (delta_energy >= 0 && exp(-delta_energy / T) <= ((double)rand() / RAND_MAX))
                    {
                        // Reject and restore backup
                        cdt = backup;
                    }
                }
            }
        }
        // Decrease temperature
        T -= 1.0 / L;
    }
}