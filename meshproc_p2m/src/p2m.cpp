#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>

#include <ros/ros.h>

#include <meshproc_msgs/Cloud2Mesh.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_msg.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Polyhedron_3.h>

#include <trimesh/TriMesh.h>
#include <trimesh/TriMesh_algo.h>

#include <meshproc_p2m/typedefs.h>
#include <meshproc_p2m/PointCloudEntry.h>

bool do_Cloud2Mesh(meshproc_msgs::Cloud2Mesh::Request &req,
                   meshproc_msgs::Cloud2Mesh::Response &res)
{
    meshproc_p2m::PointCloudEntry pce;
    res.input_found = pce.loadFromFile(req.input_filename);
    res.operation_done = false;
    if(res.input_found)
        res.operation_done = pce.cloud2Mesh(req.cell_size) && pce.writeMeshToFile(req.output_filename);
    return true;
}

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "meshproc_p2m");
  ros::NodeHandle n;

  ROS_INFO("Advertising services ...");
  ros::ServiceServer Cloud2Mesh_service = n.advertiseService("meshproc_p2m/Cloud2Mesh", do_Cloud2Mesh);
  ROS_INFO(" ... all done.");

  ros::spin();
  return 0;
}

