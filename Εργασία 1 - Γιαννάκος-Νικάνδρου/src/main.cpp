#include <iostream>
#include "./func.h"

int main()
{
    // string file_path = "test_instances/instance_test_8.json";
    // string file_path = "instance_test_0.json";
    string file_path = "../input.json";
    string instance_uid;
    int num__constraints = 0;
    vector<Point_2> points;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;

    if (read_json_file(file_path, instance_uid, points, region_boundary, num__constraints, additional_constraints))
    {
        cout << "Instance UID: " << instance_uid << endl;
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
    }

    cout << "Commencing Triangulation" << endl;
    CDT cdt;
    cdt = triangulation(points, region_boundary);
    cout << "Went Well...." << endl;

    std::string filename = "../output.json"; // Specify your desired output filename
    create_json_output(cdt, filename);

    // // Visualize the triangulation
    // CGAL::draw(cdt);

    return 0;
}