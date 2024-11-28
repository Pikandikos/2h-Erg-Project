#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Simple_cartesian.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;

void export_to_svg(const CDT &cdt, const std::string &filename)
{
    std::ofstream ofs(filename);
    if (!ofs)
    {
        std::cerr << "Error: Cannot open file " << filename << " for writing." << std::endl;
        return;
    }

    // Calculate the bounding box of the triangulation
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        auto point = vit->point();
        min_x = std::min(min_x, point.x());
        max_x = std::max(max_x, point.x());
        min_y = std::min(min_y, point.y());
        max_y = std::max(max_y, point.y());
    }

    // Scale factor to fit the entire triangulation in a 500x500 SVG canvas
    double scale = std::min(500.0 / (max_x - min_x), 500.0 / (max_y - min_y));

    // Calculate canvas size
    double width = (max_x - min_x) * scale;
    double height = (max_y - min_y) * scale;

    // Write the SVG header
    ofs << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width << "\" height=\"" << height << "\">\n";

    // Draw edges
    for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit)
    {
        auto segment = cdt.segment(*eit);

        // Apply scaling and translation
        double x1 = (segment.source().x() - min_x) * scale;
        double y1 = (segment.source().y() - min_y) * scale;
        double x2 = (segment.target().x() - min_x) * scale;
        double y2 = (segment.target().y() - min_y) * scale;

        ofs << "<line x1=\"" << x1 << "\" y1=\"" << height - y1 // Invert y-axis for SVG coordinate system
            << "\" x2=\"" << x2 << "\" y2=\"" << height - y2
            << "\" style=\"stroke:black;stroke-width:1\" />\n";
    }

    // Draw vertices
    for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
    {
        auto point = vit->point();

        // Apply scaling and translation
        double x = (point.x() - min_x) * scale;
        double y = (point.y() - min_y) * scale;

        ofs << "<circle cx=\"" << x << "\" cy=\"" << height - y // Invert y-axis for SVG coordinate system
            << "\" r=\"3\" fill=\"red\" />\n";
    }

    ofs << "</svg>\n";
    ofs.close();
    std::cout << "Triangulation exported to " << filename << std::endl;
}
