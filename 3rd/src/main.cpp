#include <iostream>
#include "./func.h"

int main()
{
    string file_path = "../test_instances/instance_test_6_2.json";
    string instance_uid;
    int num__constraints = 0;
    vector<Point_2> points;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;
    string method;
    ptree parameters;
    vector<double> method_parameters;
    bool delaunay;

    if (read_json_file(file_path, instance_uid, points, region_boundary, num__constraints, additional_constraints, method, parameters, delaunay))
    {
        cout << "Instance UID: " << instance_uid << endl;
        cout << "Method: " << method << endl;
        cout << "Delaunay: " << (delaunay ? "true" : "false") << endl;

        cout << "Points:" << endl;
        for (const auto &point : points)
        {
            cout << "(" << point.x() << ", " << point.y() << ")" << endl;
        }

        cout << "Region Boundary:" << endl;
        for (int index : region_boundary)
        {
            cout << index << " ";
        }
        cout << endl;

        cout << "Additional Constraints:" << endl;
        for (const auto &constraint : additional_constraints)
        {
            cout << "(" << constraint.first << ", " << constraint.second << ")" << endl;
        }
        cout << "Num_constraints: " << num__constraints << endl;

        cout << "Parameters for " << method << ":" << endl;
        for (const auto &param : parameters)
        {
            cout << "  " << param.first << ": " << param.second.data() << endl;
        }
    }

    populateVector(method, parameters, method_parameters);

    cout << "Commencing Triangulation" << endl;
    CDT cdt;

    cdt = initial_triangulation(points, region_boundary, additional_constraints, delaunay);
    vector<Point_2> region_points = create_region_points(points, region_boundary);

    cout << "Total obtuse angles found after initial triangulation: " << analyze_obtuse_angles(cdt) << "\n";

    if (method == "local")
    {
        local_search(cdt, region_points, method_parameters);
    }
    else if (method == "sa")
    {
        simulated_annealing(cdt, region_points, method_parameters);
    }
    else if (method == "ant")
    {
        ant_colony_optimization(cdt, region_points, method_parameters);
    }
    else
    {
        cout << "Method not recognised..." << endl;
    }

    cout << "Went Well...." << endl;

    cout << "Analyzing triangles for obtuse angles...\n";
    cout << "Total obtuse angles found: " << analyze_obtuse_angles(cdt) << "\n";

    std::string filename = "../output.json"; // Specify your desired output filename
    create_json_output(cdt, filename, instance_uid, analyze_obtuse_angles(cdt), method, method_parameters);

    export_to_svg(cdt, "output.svg");

    return 0;
}