#include <iostream>
#include "./func.h"

int main()
{
    string file_path = "../test_instances/instance_test_22_2.json";
    string instance_uid;
    int num__constraints = 0;
    vector<Point_2> points;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;
    string method;
    ptree parameters;
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

    cout << "Commencing Triangulation" << endl;
    CDT cdt;

    cdt = triangulation(points, region_boundary, additional_constraints, parameters);

    cout << "Went Well...." << endl;

    analyze_obtuse_angles(cdt);

    std::string filename = "../output.json"; // Specify your desired output filename
    create_json_output(cdt, filename);

    export_to_svg(cdt, "output.svg");

    return 0;
}