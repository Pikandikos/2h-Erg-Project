#ifndef INSERT_NO_FLIP_H
#define INSERT_NO_FLIP_H

#include <CGAL/Constrained_triangulation_2.h>

template <typename CDT>
void insert_no_flip(CDT &cdt, const typename CDT::Point &point)
{
    // Insert without triggering edge flips
    cdt.insert(point); // Insert the point, without triggering flips
}

template <typename CDT>
bool remove_constraint_no_flip(CDT &cdt, const typename CDT::Edge &edge)
{
    // Extract the face and edge index from the input edge
    auto face = edge.first;  // The face containing the edge
    int index = edge.second; // The edge index within the face

    // Check if the edge is constrained before attempting removal
    bool edge_exists_before = cdt.is_constrained(edge);
    if (!edge_exists_before)
    {
        // Edge is not constrained, nothing to remove
        return false;
    }

    // Remove the constrained edge
    cdt.remove_constraint(face, index);

    // Verify the edge is no longer constrained
    bool edge_exists_after = cdt.is_constrained(edge);
    return !edge_exists_after; // Return true if removal was successful
}

#endif // INSERT_NO_FLIP_H
