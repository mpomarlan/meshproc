#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>

#include <ros/ros.h>

#include <meshproc_msgs/PolygonConvexDecomposition.h>
#include <meshproc_msgs/PolygonCSGRequest.h>
#include <meshproc_msgs/PolygonVisibility.h>
#include <meshproc_msgs/GetPolygonSkeleton.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_msg.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>

#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>

#include <meshproc_2dg/typedefs.h>
#include <meshproc_2dg/PolygonEntry.h>

bool do_PolygonConvexDecomposition(meshproc_msgs::PolygonConvexDecomposition::Request &req,
                                   meshproc_msgs::PolygonConvexDecomposition::Response &res)
{
    meshproc_2dg::Polygon_with_holes_2 polygon;
    //meshproc_2dg::Polygon_2 polygon;
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon, polygon);
    std::vector<meshproc_2dg::Polygon_2> results;
    res.operation_performed = meshproc_2dg::PolygonEntry::convexDecomposition(results, polygon, req.ignore_holes, req.triangulate);
    int maxK = results.size();
    res.convex_parts.clear();
    res.convex_parts.resize(maxK);
    for(int k = 0; k < maxK; k++)
    {

        meshproc_2dg::PolygonEntry::writeToMsg(res.convex_parts[k], results[k]);
    }
    return true;
}

bool do_PolygonCSGRequest(meshproc_msgs::PolygonCSGRequest::Request &req,
                          meshproc_msgs::PolygonCSGRequest::Response &res)
{
    meshproc_2dg::Polygon_with_holes_2 polygon_A, polygon_B;
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon_A, polygon_A);
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon_B, polygon_B);
    std::list<meshproc_2dg::Polygon_with_holes_2> result;
    res.operation_performed = meshproc_2dg::PolygonEntry::csgRequest(polygon_A, polygon_B, result, req.operation);
    int maxK = result.size();
    std::list<meshproc_2dg::Polygon_with_holes_2>::const_iterator it = result.begin();
    res.polygons_R.resize(maxK);
    for(int k = 0; k < maxK; k++, it++)
        meshproc_2dg::PolygonEntry::writeToMsg(res.polygons_R[k], *it);
    return true;
}

bool do_PolygonVisibility(meshproc_msgs::PolygonVisibility::Request &req,
                          meshproc_msgs::PolygonVisibility::Response &res)
{
    meshproc_2dg::Polygon_with_holes_2 polygon;
    std::vector<meshproc_2dg::Point_2> points;
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon, polygon);
    int maxK = req.points.size();
    points.resize(maxK);
    for(int k = 0; k < maxK; k++)
    {
        points[k] = meshproc_2dg::Point_2(req.points[k].x, req.points[k].y);
    }
    std::vector<meshproc_2dg::Polygon_2> results;
    res.operation_performed = meshproc_2dg::PolygonEntry::visibility(polygon, points, results);
    maxK = results.size();
    res.results.resize(maxK);
    for(int k = 0; k < maxK; k++)
    {
        meshproc_2dg::PolygonEntry::writeToMsg(res.results[k], results[k]);
    }
    return true;
}

bool do_GetPolygonSkeleton(meshproc_msgs::GetPolygonSkeleton::Request &req,
                           meshproc_msgs::GetPolygonSkeleton::Response &res)
{
    meshproc_2dg::Polygon_with_holes_2 polygon;
    meshproc_2dg::PolygonEntry::loadFromMsg(req.polygon, polygon);
    meshproc_2dg::SkeletonPtr skeleton = CGAL::create_interior_straight_skeleton_2(polygon);
    meshproc_2dg::PolygonEntry::writeToMsg(res.skeleton, *skeleton);
    return true;
}

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "meshproc_2dg");
  ros::NodeHandle n;

  ROS_INFO("Advertising services ...");
  ros::ServiceServer PolygonConvexDecomposition_service = n.advertiseService("meshproc_2dg/PolygonConvexDecomposition", do_PolygonConvexDecomposition);
  ros::ServiceServer PolygonCSGRequest_service = n.advertiseService("meshproc_2dg/PolygonCSGRequest", do_PolygonCSGRequest);
  ros::ServiceServer PolygonVisibility_service = n.advertiseService("meshproc_2dg/PolygonVisibility", do_PolygonVisibility);
  ros::ServiceServer PolygonSkeleton_service = n.advertiseService("meshproc_2dg/GetPolygonSkeleton", do_GetPolygonSkeleton);
  ROS_INFO(" ... all done.");

  ros::spin();
  return 0;
}

