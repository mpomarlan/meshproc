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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/IO/read_off_points.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <CGAL/property_map.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/wlop_simplify_and_regularize_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/bilateral_smooth_point_set.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/tuple.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <trimesh/TriMesh.h>
#include <trimesh/TriMesh_algo.h>

#include <meshproc_p2m/typedefs.h>
#include <meshproc_p2m/PointCloudEntry.h>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>

namespace meshproc_p2m
{

struct Perimeter {
  double bound;
  Perimeter(double bound)
    : bound(bound)
  {}
  // The point type that will be injected here will be
  // CGAL::Exact_predicates_inexact_constructions_kernel::Point_3
  template <typename Point>
  bool operator()(const Point& p, const Point& q, const Point& r) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return false;
    }
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};

bool PointCloudEntry::loadFromFile(std::string const& filename)
{
    points.clear();

    std::string extension(filename.substr(filename.find_last_of(".") + 1));
    for(int k = 0; k < extension.size(); k++)
        extension[k] = std::tolower(extension[k]);

    std::cerr << "Opening file " << filename << std::endl;

    FILE *b = fopen(filename.c_str(), "rt");
    if(NULL == b)
    {
        std::cerr << "    File could not be opened for reading. Check if it exists and has the right permissions" << std::endl;
        return false;
    }

    if(0 == extension.compare("xyz"))
    {
        std::ifstream stream(filename.c_str());
        CGAL::read_xyz_points(stream, std::back_inserter(points));
    }
    else if(0 == extension.compare("off"))
    {
        std::ifstream stream(filename.c_str());
        CGAL::read_off_points(stream, std::back_inserter(points));
    }
    else
    {
        pcl::PointCloud<pcl::PointXYZ> cloud;
        pcl::io::loadPCDFile(filename, cloud);
        int maxK = cloud.size();
        for(int k = 0; k < maxK; k++)
        {
            points.push_back(Point(cloud[k].data[0], cloud[k].data[1], cloud[k].data[2]));
        }
    }
    return true;
}

bool PointCloudEntry::cloud2Mesh(double cellSize)
{
    mesh.clear();
    // Removes outliers using erase-remove idiom.
    // The Identity_property_map property map can be omitted here as it is the default value.
    const double removed_percentage = 5.0; // percentage of points to remove
    const int nb_neighbors = 30; // considers 24 nearest neighbor points
    points.erase(CGAL::remove_outliers(points.begin(), points.end(),
                                       CGAL::Identity_property_map<Point>(),
                                       nb_neighbors, removed_percentage),
                 points.end());

    std::cout << "\tgrid based simplification ..." << std::endl;
    // simplification by clustering using erase-remove idiom
    points.erase(CGAL::grid_simplify_point_set(points.begin(), points.end(), cellSize),
                 points.end());

    //mesh.insert(points.begin(), points.end());
    //mesh.increase_scale(1);
    //mesh.reconstruct_surface();
    facets.clear();

    Perimeter perimeter(0);
    CGAL::advancing_front_surface_reconstruction(points.begin(), points.end(), std::back_inserter(facets), perimeter);
    return true;
}

bool PointCloudEntry::writeMeshToFile(double gridsize, std::string const& filename) const
{
    trimesh::TriMesh M;
    PointCloudEntry::write_to_trimesh(gridsize, facets, points, &M);
    return M.write(filename.c_str());
}

double getDistance(Point_collection const& points, int A, int B)
{
    double xA, yA, zA, xB, yB, zB;
    xA = ::CGAL::to_double(points[A].x());
    yA = ::CGAL::to_double(points[A].y());
    zA = ::CGAL::to_double(points[A].z());
    xB = ::CGAL::to_double(points[B].x());
    yB = ::CGAL::to_double(points[B].y());
    zB = ::CGAL::to_double(points[B].z());
    double dx = xA - xB;
    double dy = yA - yB;
    double dz = zA - zB;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

bool PointCloudEntry::write_to_trimesh(double gridsize, std::vector<AFSFacet> const& R, Point_collection const& P, trimesh::TriMesh *M)
{
    M->vertices.clear();
    M->faces.clear();
    M->vertices.reserve(P.size());
    M->faces.reserve(R.size());
    for(Point_collection::const_iterator it = P.begin();
        it != P.end(); it++)
    {
        M->vertices.push_back(trimesh::point(::CGAL::to_double(it->x()),
                                             ::CGAL::to_double(it->y()),
                                             ::CGAL::to_double(it->z())));

    }

    double dt = 2.5*gridsize;
    for(std::vector<AFSFacet>::const_iterator it = R.begin(); it != R.end(); it++)
    {
        double dAB, dBC, dAC;
        dAB = getDistance(P, (*it)[0], (*it)[1]);
        dBC = getDistance(P, (*it)[1], (*it)[2]);
        dAC = getDistance(P, (*it)[0], (*it)[2]);
        if((dAB < dt) && (dBC < dt) && (dAC < dt))
            M->faces.push_back(trimesh::TriMesh::Face((*it)[0], (*it)[1], (*it)[2]));
    }
    return true;
}

}
