#include <cmath>
#include <list>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_msg.h>

#include <geometry_msgs/Polygon.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/partition_2.h>

#include <meshproc_2dg/typedefs.h>
#include <meshproc_2dg/PolygonEntry.h>

namespace meshproc_2dg
{

bool PolygonEntry::loadFromMsg(geometry_msgs::Polygon const& msg)
{
    return PolygonEntry::loadFromMsg(msg, polygon);
}

bool PolygonEntry::writeToMsg(geometry_msgs::Polygon & msg) const
{
    return PolygonEntry::writeToMsg(msg, polygon);
}

bool PolygonEntry::loadFromMsg(const geometry_msgs::Polygon &msg, Polygon_2 &polygon)
{
    polygon.clear();
    int maxK = msg.points.size();
    for(int k = 0; k < maxK; k++)
        polygon.push_back(Point_2(msg.points[k].x, msg.points[k].y));
    return true;
}

bool PolygonEntry::writeToMsg(geometry_msgs::Polygon &msg, Polygon_2 const& polygon)
{
    msg.points.clear(); msg.points.reserve(polygon.size());
    for(Polygon_2::Vertex_const_iterator it = polygon.vertices_begin();
        it != polygon.vertices_end(); it++)
    {
        geometry_msgs::Point32 aux;
        aux.x = ::CGAL::to_double(it->x());
        aux.y = ::CGAL::to_double(it->y());
        msg.points.push_back(aux);
    }
    return true;
}

bool PolygonEntry::convexDecomposition(std::vector<Polygon_2> & results) const
{

    return PolygonEntry::convexDecomposition(results, polygon);
}

bool PolygonEntry::convexDecomposition(std::vector<Polygon_2> & results, Polygon_2 const& polygon)
{
    Traits partitionTraits;
    Validity_traits validityTraits;
    results.clear();
    results.reserve(polygon.size());
    if(polygon.is_clockwise_oriented())
    {
        Polygon_2 polygonRev = polygon;
        polygonRev.reverse_orientation();
        CGAL::optimal_convex_partition_2(polygonRev.vertices_begin(), polygonRev.vertices_end(),
                                         std::back_inserter(results), partitionTraits);
    }
    else
        CGAL::optimal_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(),
                                         std::back_inserter(results), partitionTraits);

    return CGAL::partition_is_valid_2(polygon.vertices_begin(), polygon.vertices_end(),
                                      results.begin(), results.end(), validityTraits);
}

}
