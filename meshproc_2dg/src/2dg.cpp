#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>

#include <ros/ros.h>

#include <meshproc_msgs/PolygonConvexDecomposition.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_msg.h>

#include <geometry_msgs/Polygon.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>

#include <meshproc_2dg/typedefs.h>
#include <meshproc_2dg/PolygonEntry.h>

bool do_PolygonConvexDecomposition(meshproc_msgs::PolygonConvexDecomposition::Request &req,
                                   meshproc_msgs::PolygonConvexDecomposition::Response &res)
{
    meshproc_2dg::Polygon_2 polygon;
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon, polygon);
    std::vector<meshproc_2dg::Polygon_2> results;
    res.operation_performed = meshproc_2dg::PolygonEntry::convexDecomposition(results, polygon);
    int maxK = results.size();
    res.convex_parts.clear();
    res.convex_parts.resize(maxK);
    for(int k = 0; k < maxK; k++)
    {

        meshproc_2dg::PolygonEntry::writeToMsg(res.convex_parts[k], results[k]);
    }
    return true;
}

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "meshproc_2dg");
  ros::NodeHandle n;

  ROS_INFO("Advertising services ...");
  ros::ServiceServer PolygonConvexDecomposition_service = n.advertiseService("meshproc_2dg/PolygonConvexDecomposition", do_PolygonConvexDecomposition);
  ROS_INFO(" ... all done.");

  ros::spin();
  return 0;
}

