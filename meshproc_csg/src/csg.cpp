#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>

#include <trimesh/TriMesh.h>
#include <trimesh/TriMesh_algo.h>

#include <ros/ros.h>

#include <meshproc_msgs/CSGRequest.h>
#include <meshproc_msgs/GetLoadedMeshNames.h>
#include <meshproc_msgs/GetMeshProps.h>
#include <meshproc_msgs/LoadMesh.h>
#include <meshproc_msgs/UnloadMesh.h>
#include <meshproc_msgs/AffineTransformMesh.h>
#include <meshproc_msgs/GetMeshAABB.h>
#include <meshproc_msgs/GetNearMeshVertices.h>
#include <meshproc_msgs/GetMesh.h>
#include <meshproc_msgs/ConvexDecomposition.h>
#include <meshproc_msgs/GetMeshSkeleton.h>
#include <meshproc_msgs/ConvexHull.h>
#include <meshproc_msgs/ProjectMesh.h>
#include <meshproc_msgs/SolidifyMesh.h>

#include <meshproc_csg/kdtree++/kdtree.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <eigen_conversions/eigen_msg.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>

#include <meshproc_csg/typedefs.h>

#include <meshproc_csg/csg.h>

using namespace meshproc_csg;

/*Magic numbers here: default parameter values for the various stages. To be overwritten
  by rosparams/service request messages, if these provide alternative values.*/

MeshMap loadedMeshes;

bool do_LoadMesh(meshproc_msgs::LoadMesh::Request &req, meshproc_msgs::LoadMesh::Response &res)
{
    ROS_INFO("Loading mesh %s", req.mesh_name.c_str());
    res.loaded_mesh = res.io_error = res.mesh_already_loaded = false;
    MeshMap::iterator it = loadedMeshes.find(req.mesh_name);
    if(it != loadedMeshes.end())
        res.mesh_already_loaded = true;
    else
    {
        if((0 == req.mesh_filenames.size()) && (0 == req.mesh_msgs.size()))
            return true;
        it = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_name, new MeshEntry()));
        MeshEntry *R = it->second;
        int maxK = req.mesh_filenames.size();
        bool goOn = true;
        for(int k = 0; (k < maxK) && goOn; k++)
        {
            ROS_INFO("    loading part from file %s", req.mesh_filenames[k].c_str());
            goOn = R->loadFromFile(req.mesh_filenames[k], req.duplicate_dist);
        }
        if(!goOn)
        {
            ROS_INFO("    Encountered I/O error (file might not exist or is inaccessible), cancelling load.");
            res.io_error = true;
            delete R;
            loadedMeshes.erase(it);
            return true;
        }
        res.loaded_mesh = true;
        maxK = req.mesh_msgs.size();
        for(int k = 0; k < maxK; k++)
        {
            ROS_INFO("    loading part from message");
            R->loadFromMsg(req.mesh_msgs[k], req.duplicate_dist);
        }
    }
    ROS_INFO("    Loading done.");
    return true;
}

bool do_UnloadMesh(meshproc_msgs::UnloadMesh::Request &req, meshproc_msgs::UnloadMesh::Response &res)
{
    res.unloaded_mesh = 0;
    MeshMap::iterator it = loadedMeshes.find(req.mesh_name);
    if(loadedMeshes.end() != it)
    {
        delete it->second;
        res.unloaded_mesh = 1;
        loadedMeshes.erase(it);
    }
    return true;
}

bool do_GetLoadedMeshNames(meshproc_msgs::GetLoadedMeshNames::Request &req,
                           meshproc_msgs::GetLoadedMeshNames::Response &res)
{
    res.mesh_names.clear();
    MeshMap::iterator it = loadedMeshes.begin();
    for(it = loadedMeshes.begin(); it != loadedMeshes.end(); it++)
    {
        res.mesh_names.push_back(it->first);
    }
    return true;
}

bool do_GetMeshProps(meshproc_msgs::GetMeshProps::Request &req, meshproc_msgs::GetMeshProps::Response &res)
{
    MeshMap::iterator it = loadedMeshes.find(req.mesh_name);
    if(it != loadedMeshes.end())
    {
        res.is_loaded = true;
        res.is_closed = it->second->isClosed();
        res.is_manifold = it->second->isManifold();
        res.is_orientable = it->second->isOrientable();
        res.is_safe_for_csg = it->second->isCSGSafe();
        res.num_connected_components = it->second->getNrConnectedComponents();
        res.num_vertices = it->second->getNrVertices();
        res.num_edges = it->second->getNrEdges();
        res.num_faces = it->second->getNrFaces();
        res.Euler_characteristic = it->second->getEulerCharacteristic();
    }
    else
    {
        res.is_loaded = false;
    }
    return true;
}

MeshEntry* checkMeshAvailability(std::string const& meshName,
                                 meshproc_msgs::CSGRequest::Response::_mesh_A_loaded_type &isLoaded,
                                 meshproc_msgs::CSGRequest::Response::_mesh_A_csg_safe_type &isCSGSafe)
{
    isLoaded = isCSGSafe = false;
    MeshMap::iterator it = loadedMeshes.find(meshName);
    if(it != loadedMeshes.end())
    {
        isLoaded = true;
        isCSGSafe = it->second->isCSGSafe();
        return it->second;
    }
    return NULL;
}

bool canPerform(bool A_loaded, bool B_loaded, bool A_CSGSafe, bool B_CSGSafe, int operation)
{
    switch(operation)
    {
    case 0:
    case 1:
    case 2:
    case 3:
        return(A_CSGSafe && B_CSGSafe);
    case 4:
        return(A_loaded && B_loaded);
    case 5:
        return(A_CSGSafe && B_loaded);
    }
    return false;
}

bool do_CSGRequest(meshproc_msgs::CSGRequest::Request &req, meshproc_msgs::CSGRequest::Response &res)
{
    res.operation_performed = false;
    res.mesh_A_csg_safe = res.mesh_A_loaded = res.mesh_B_csg_safe = res.mesh_B_loaded = false;
    res.result = shape_msgs::Mesh();
    bool already_had_R = true;
    MeshEntry *A, *B, *R;
    A = checkMeshAvailability(req.mesh_A, res.mesh_A_loaded, res.mesh_A_csg_safe);
    B = checkMeshAvailability(req.mesh_B, res.mesh_B_loaded, res.mesh_B_csg_safe);
    MeshMap::iterator it = loadedMeshes.find(req.mesh_R);
    if(it == loadedMeshes.end())
    {
        already_had_R = false;
        it = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_R, new MeshEntry()));
    }
    R = it->second;
    if(!canPerform(res.mesh_A_loaded, res.mesh_B_loaded, res.mesh_A_csg_safe, res.mesh_B_csg_safe, req.operation))
    {
        //Operation can't be performed, so let's not leave a result mesh (if created just now) loaded as a side-effect.
        if(!already_had_R)
        {
            delete R;
            loadedMeshes.erase(req.mesh_R);
        }
        return true;
    }
    switch(req.operation)
    {
    case 0:
        res.operation_performed = R->setFromUnion(*A, *B);
        break;
    case 1:
        res.operation_performed = R->setFromIntersection(*A, *B);
        break;
    case 2:
        res.operation_performed = R->setFromDifference(*A, *B);
        break;
    case 3:
        res.operation_performed = R->setFromSymmetricDifference(*A, *B);
        break;
    case 4:
        res.operation_performed = R->setFromMinkowskiSum(*A, *B);
        break;
    case 5:
        res.operation_performed = R->setFromMinkowskiErosion(*A, *B);
        break;
    }

    if(!res.operation_performed)
    {
        //Operation failed, so let's not leave a result mesh (if created just now) loaded as a side-effect.
        if(!already_had_R)
        {
            delete R;
            loadedMeshes.erase(req.mesh_R);
        }
        return true;
    }

    if(req.return_result)
        R->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.result_to_file)
    {
        res.file_written = R->writeToFile(req.result_filename);
    }
    return true;
}

bool do_AffineTransformMesh(meshproc_msgs::AffineTransformMesh::Request &req,
                            meshproc_msgs::AffineTransformMesh::Response &res)
{
    MeshEntry *A, *R;
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_A);
    MeshMap::iterator itR = loadedMeshes.find(req.mesh_R);
    res.mesh_A_loaded = true;
    res.file_written = false;
    if(itA == loadedMeshes.end())
    {
        res.mesh_A_loaded = false;
        return true;
    }
    if(itR == loadedMeshes.end())
    {
        itR = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_R, new MeshEntry()));
    }
    A = itA->second;
    R = itR->second;
    Eigen::Affine3d M = Eigen::Translation3d(req.transform.translation.x,
                                             req.transform.translation.y,
                                             req.transform.translation.z)*
                        Eigen::Quaterniond(req.transform.rotation.w,
                                           req.transform.rotation.x,
                                           req.transform.rotation.y,
                                           req.transform.rotation.z);
    R->applyTransform(M);
    if(req.return_result)
        R->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.result_to_file)
    {
        res.file_written = R->writeToFile(req.result_filename);
    }
    return true;
}

bool do_GetMeshAABB(meshproc_msgs::GetMeshAABB::Request &req,
                    meshproc_msgs::GetMeshAABB::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_name);
    res.mesh_loaded = true;
    if(itA == loadedMeshes.end())
    {
        res.mesh_loaded = false;
        return true;
    }
    double maxX, minX, maxY, minY, maxZ, minZ;
    itA->second->getBoundingBox(maxX, minX, maxY, minY, maxZ, minZ);
    res.max_x = maxX;
    res.min_x = minX;
    res.max_y = maxY;
    res.min_y = minY;
    res.max_z = maxZ;
    res.min_z = minZ;
    return true;
}

bool do_GetMesh(meshproc_msgs::GetMesh::Request &req,
                meshproc_msgs::GetMesh::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_name);
    res.mesh_loaded = true;
    res.file_written = false;
    if(itA == loadedMeshes.end())
    {
        res.mesh_loaded = false;
        return true;
    }
    if(req.return_result)
        itA->second->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.result_to_file)
    {
        res.file_written = itA->second->writeToFile(req.result_filename);
    }
    return true;
}

bool do_GetNearMeshVertices(meshproc_msgs::GetNearMeshVertices::Request &req,
                            meshproc_msgs::GetNearMeshVertices::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_name);
    res.mesh_loaded = true;
    res.neighbors.clear();
    if(itA == loadedMeshes.end())
    {
        res.mesh_loaded = false;
        return true;
    }
    std::vector<MeshEntry::XYZTriplet> points; points.clear();
    itA->second->getNearVertices(req.point.x, req.point.y, req.point.z, req.distance, points);
    int maxK = points.size();
    for(int k = 0; k < maxK; k++)
    {
        geometry_msgs::Point aux;
        aux.x = points[k].x;
        aux.y = points[k].y;
        aux.z = points[k].z;
        res.neighbors.push_back(aux);
    }
    return true;
}

bool do_ConvexDecomposition(meshproc_msgs::ConvexDecomposition::Request &req,
                            meshproc_msgs::ConvexDecomposition::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_A);
    res.mesh_loaded = true;
    res.nr_components = 0;
    res.mesh_names.clear();
    res.file_written.clear();
    res.components.clear();
    if(itA == loadedMeshes.end())
    {
        res.mesh_loaded = false;
        return true;
    }
    std::vector<MeshEntry*> components; components.clear();
    std::cerr << "Getting components" << std::endl;
    itA->second->getConvexComponents(components);
    std::cerr << "Components obtained" << std::endl;
    int maxK = res.nr_components = components.size();
    for(int k = 0; k < maxK; k++)
    {
        char nr_str[20];
        std::sprintf(nr_str, "%d", k);
        std::string cName = req.mesh_R + nr_str;
        MeshMap::iterator itR = loadedMeshes.find(cName);
        if(itR == loadedMeshes.end())
        {
            itR = loadedMeshes.insert(loadedMeshes.begin(),
                                      std::pair<std::string, MeshEntry*>(cName, components[k]));
        }
        else
        {
            delete itR->second;
            itR->second = components[k];
        }
        MeshEntry* R = itR->second;
        res.mesh_names.push_back(cName);
        if(req.return_result)
        {
            res.components.push_back(shape_msgs::Mesh());
            R->writeToMsg(res.components[k]);
        }
        if(req.result_to_file)
        {
            res.file_written.push_back(R->writeToFile(req.result_filename + nr_str + req.result_filename_extension));
            res.file_names.push_back(req.result_filename + nr_str + req.result_filename_extension);
        }

    }
    return true;
}

bool do_GetMeshSkeleton(meshproc_msgs::GetMeshSkeleton::Request &req,
                        meshproc_msgs::GetMeshSkeleton::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_name);
    res.mesh_loaded = true;
    res.vertices.clear();
    res.edge_ends_L.clear();
    res.edge_ends_R.clear();
    if(itA == loadedMeshes.end())
    {
        res.mesh_loaded = false;
        return true;
    }
    std::vector<double> x, y, z;
    x.clear(); y.clear(); z.clear();
    itA->second->getMeshSkeleton(req.approximate, req.duplicate_distance, x, y, z, res.edge_ends_L, res.edge_ends_R);
    int maxK = x.size();
    res.vertices.reserve(maxK);
    for(int k = 0; k < maxK; k++)
    {
        geometry_msgs::Point aux;
        aux.x = x[k];
        aux.y = y[k];
        aux.z = z[k];
        res.vertices.push_back(aux);
    }
    return true;
}

bool do_ConvexHull(meshproc_msgs::ConvexHull::Request &req,
                   meshproc_msgs::ConvexHull::Response &res)
{
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_A);
    res.mesh_A_loaded = true;
    res.file_written = false;
    res.result = shape_msgs::Mesh();
    if(itA == loadedMeshes.end())
    {
        res.mesh_A_loaded = false;
        return true;
    }

    MeshMap::iterator itR = loadedMeshes.find(req.mesh_R);
    if(itR == loadedMeshes.end())
    {
        itR = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_R, new MeshEntry()));
    }
    MeshEntry *A, *R;
    A = itA->second;
    R = itR->second;
    R->setFromConvexHull(*A);
    if(req.return_result)
        R->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.result_to_file)
    {
        res.file_written = R->writeToFile(req.result_filename);
    }
    return true;
}

bool do_ProjectMesh(meshproc_msgs::ProjectMesh::Request &req,
                    meshproc_msgs::ProjectMesh::Response &res)
{
    bool already_had_R = true;
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_A);
    res.mesh_A_loaded = true;
    res.operation_performed = false;
    res.file_written = false;
    res.result = shape_msgs::Mesh();
    if(itA == loadedMeshes.end())
    {
        res.mesh_A_loaded = false;
        return true;
    }

    MeshMap::iterator itR = loadedMeshes.find(req.mesh_R);
    if(itR == loadedMeshes.end())
    {
        already_had_R = false;
        itR = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_R, new MeshEntry()));
    }
    MeshEntry *A, *R;
    A = itA->second;
    R = itR->second;
    res.operation_performed = R->setFromProjection(*A, req.normal.x, req.normal.y, req.normal.z, req.fill_holes);
    if(!res.operation_performed)
    {
        if(!already_had_R)
        {
            delete R;
            loadedMeshes.erase(req.mesh_R);
        }
        return true;
    }
    if(req.return_result_as_mesh)
        R->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.return_result_as_polygon)
    {
        std::vector<double> x, y, z;
        R->getBoundaryPolygon(x, y, z, res.edge_L, res.edge_R);
        int maxK = x.size();
        for(int k = 0; k < maxK; k++)
        {
            geometry_msgs::Point aux;
            aux.x = x[k];
            aux.y = y[k];
            aux.z = z[k];
            res.points.push_back(aux);
        }
    }
    if(req.result_to_file)
    {
        res.file_written = R->writeToFile(req.result_filename);
    }
    return true;
}

bool do_SolidifyMesh(meshproc_msgs::SolidifyMesh::Request &req,
                    meshproc_msgs::SolidifyMesh::Response &res)
{
    bool already_had_R = true;
    MeshMap::iterator itA = loadedMeshes.find(req.mesh_A);
    res.mesh_A_loaded = true;
    res.file_written = false;
    res.result = shape_msgs::Mesh();
    if(itA == loadedMeshes.end())
    {
        res.mesh_A_loaded = false;
        return true;
    }

    MeshMap::iterator itR = loadedMeshes.find(req.mesh_R);
    if(itR == loadedMeshes.end())
    {
        already_had_R = false;
        itR = loadedMeshes.insert(loadedMeshes.begin(),
                                 std::pair<std::string, MeshEntry*>(req.mesh_R, new MeshEntry()));
    }
    MeshEntry *A, *R;
    A = itA->second;
    R = itR->second;
    R->setFromSolidification(*A, req.thickness);
    if(req.return_result)
        R->writeToMsg(res.result);
    else
        res.result = shape_msgs::Mesh();
    if(req.result_to_file)
    {
        res.file_written = R->writeToFile(req.result_filename);
    }
    return true;
}

int main(int argc, char *argv[])
{
  ros::init(argc, argv, "meshproc_csg");
  ros::NodeHandle n;

  ROS_INFO("Advertising services ...");
  ros::ServiceServer loadMesh_service = n.advertiseService("meshproc_csg/LoadMesh", do_LoadMesh);
  ros::ServiceServer unloadMesh_service = n.advertiseService("meshproc_csg/UnloadMesh", do_UnloadMesh);
  ros::ServiceServer getLoadedMeshNams_service = n.advertiseService("meshproc_csg/GetLoadedMeshNames", do_GetLoadedMeshNames);
  ros::ServiceServer getMeshProps_service = n.advertiseService("meshproc_csg/GetMeshProps", do_GetMeshProps);
  ros::ServiceServer CSGRequest_service = n.advertiseService("meshproc_csg/CSGRequest", do_CSGRequest);
  ros::ServiceServer AffineTransformMesh_service = n.advertiseService("meshproc_csg/AffineTransformMesh", do_AffineTransformMesh);
  ros::ServiceServer GetMeshAABB_service = n.advertiseService("meshproc_csg/GetMeshAABB", do_GetMeshAABB);
  ros::ServiceServer GetNearMeshVertices_service = n.advertiseService("meshproc_csg/GetNearMeshVertices", do_GetNearMeshVertices);
  ros::ServiceServer GetMesh_service = n.advertiseService("meshproc_csg/GetMesh", do_GetMesh);
  ros::ServiceServer ConvexDecomposition_service = n.advertiseService("meshproc_csg/ConvexDecomposition", do_ConvexDecomposition);
  ros::ServiceServer GetMeshSkeleton_service = n.advertiseService("meshproc_csg/GetMeshSkeleton", do_GetMeshSkeleton);
  ros::ServiceServer ConvexHull_service = n.advertiseService("meshproc_csg/ConvexHull", do_ConvexHull);
  ros::ServiceServer ProjectMesh_service = n.advertiseService("meshproc_csg/ProjectMesh", do_ProjectMesh);
  ros::ServiceServer SolidifyMesh_service = n.advertiseService("meshproc_csg/SolidifyMesh", do_SolidifyMesh);
  ROS_INFO(" ... all done.");

  ros::spin();
  return 0;
}


