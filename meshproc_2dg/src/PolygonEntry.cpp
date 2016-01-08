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

#include <meshproc_msgs/PolygonWithHoles.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/Partition_is_valid_traits_2.h>
#include <CGAL/polygon_function_objects.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/partition_2.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_batched_point_location.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Polygon_convex_decomposition_2.h>
#include <CGAL/Polygon_vertical_decomposition_2.h>
#include <CGAL/Polygon_triangulation_decomposition_2.h>

#include <CGAL/Triangular_expansion_visibility_2.h>

#include <meshproc_2dg/typedefs.h>
#include <meshproc_2dg/PolygonEntry.h>

namespace meshproc_2dg
{

typedef CGAL::Segment_2<Kernel>                           Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                Arr_traits_2;
typedef CGAL::Arrangement_2<Arr_traits_2>                 Arrangement_2;
typedef Arrangement_2::Halfedge_const_handle              Halfedge_const_handle;
typedef Arrangement_2::Face_handle                        Face_handle;
typedef Arrangement_2::Face_const_handle                  Face_const_handle;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2> TEV;
typedef CGAL::Arr_point_location_result<Arrangement_2>    Point_location_result;
typedef std::pair<Point_2, Point_location_result::Type>   Query_result;


bool PolygonEntry::loadFromMsg(meshproc_msgs::PolygonWithHoles const& msg)
{
    return PolygonEntry::loadFromMsg(msg, polygon);
}

bool PolygonEntry::writeToMsg(meshproc_msgs::PolygonWithHoles & msg) const
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

bool PolygonEntry::loadFromMsg(const meshproc_msgs::PolygonWithHoles &msg, Polygon_with_holes_2 &polygon)
{
    polygon.clear();
    Polygon_2 boundary;
    PolygonEntry::loadFromMsg(msg.boundary, boundary);
    polygon = Polygon_with_holes_2(boundary);
    int maxK = msg.holes.size();
    for(int k = 0; k < maxK; k++)
    {
        Polygon_2 hole;
        PolygonEntry::loadFromMsg(msg.holes[k], hole);
        polygon.add_hole(hole);
    }
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

bool PolygonEntry::writeToMsg(meshproc_msgs::PolygonWithHoles &msg, Polygon_with_holes_2 const& polygon)
{
    PolygonEntry::writeToMsg(msg.boundary, polygon.outer_boundary());
    int k = 0;
    msg.holes.resize(polygon.number_of_holes());
    for(Polygon_with_holes_2::Hole_const_iterator it = polygon.holes_begin(); it != polygon.holes_end(); it++)
    {
        PolygonEntry::writeToMsg(msg.holes[k], *it);
        k++;
    }
    return true;
}

bool PolygonEntry::convexDecomposition(std::vector<Polygon_2> & results, bool ignoreHoles, bool triangulate) const
{
    return PolygonEntry::convexDecomposition(results, polygon, ignoreHoles, triangulate);
}

Polygon_2 PolygonEntry::Polygon_2_Conversion(Traits::Polygon_2 const& polygon)
{
    Polygon_2 retq;
    retq.clear();
    for(Traits::Polygon_2::Vertex_const_iterator it = polygon.vertices_begin(); it != polygon.vertices_end(); it++)
    {
        //Polygon_2::Point_2 aux;
        retq.push_back(*it);
    }
    return retq;
}

bool PolygonEntry::convexDecomposition(std::vector<Polygon_2> & results, Polygon_with_holes_2 const& polygon, bool ignoreHoles, bool triangulate)
{
    results.clear();
    if(ignoreHoles)
    {
        if(triangulate)
        {
            CGAL::Polygon_triangulation_decomposition_2<Kernel> triangulator;
            triangulator(polygon.outer_boundary(), std::back_inserter(results));
        }
        else
        {
            Polygon_2 boundary = polygon.outer_boundary();
            std::list<Traits::Polygon_2> decomp; decomp.clear();
            Traits traits;
            CGAL::optimal_convex_partition_2(boundary.vertices_begin(),
                                             boundary.vertices_end(),
                                             std::back_inserter(decomp), traits);
            for(std::list<Traits::Polygon_2>::const_iterator it = decomp.begin(); it!= decomp.end(); it++)
            {
                results.push_back(PolygonEntry::Polygon_2_Conversion(*it));
            }
        }
    }
    else
    {
        if(triangulate)
        {
            CGAL::Polygon_triangulation_decomposition_2<Kernel> triangulator;
            triangulator(polygon, std::back_inserter(results));
        }
        else
        {
            CGAL::Polygon_vertical_decomposition_2<Kernel> vertDecomp;
            vertDecomp(polygon, std::back_inserter(results));
        }
    }
    return true;
}

bool PolygonEntry::csgRequest(Polygon_with_holes_2 const& polygon_A, Polygon_with_holes_2 const& polygon_B, std::list<Polygon_with_holes_2> & result, int operation)
{

    CGAL::Polygon_vertical_decomposition_2<Kernel> decomp;
    switch(operation)
    {
        case 0: //Union
          result.push_back(Polygon_with_holes_2());
          CGAL::join(polygon_A, polygon_B, *(result.begin()));
          break;
        case 1: //Intersection
          CGAL::intersection(polygon_A, polygon_B, std::back_inserter(result));
          break;
        case 2: //Difference
          CGAL::difference(polygon_A, polygon_B, std::back_inserter(result));
          break;
        case 3: //Symmetric difference
            CGAL::symmetric_difference(polygon_A, polygon_B, std::back_inserter(result));
          break;
        case 4: //Minkowski sum
            result.push_back(CGAL::minkowski_sum_2(polygon_A, polygon_B, decomp));
          break;
        case 5: //Minkowski erosion
            return false; //Not supported yet.
          break;
        default:
          return false;
    }
    return true;
}

bool compareBoundary(Polygon_2 const& polygon_A, Arrangement_2::Ccb_halfedge_const_circulator const& ccb)
{
    Arrangement_2::Ccb_halfedge_const_circulator circ = ccb;
    Polygon_2 polygon_B;
    polygon_B.clear();
    do
    {
        polygon_B.push_back(circ->target()->point());
        circ++;
    }while(circ != ccb);

    int maxK_A = polygon_A.size();
    int maxK_B = polygon_B.size();
    if(maxK_A != maxK_B)
        return false;
    bool found_first = false;
    int index_0;
    for(int k = 0; (!found_first) && (k < maxK_A); k++)
    {
        found_first = (polygon_A.vertex(k) == polygon_B.vertex(0));
        if(found_first)
            index_0 = k;
    }
    if(!found_first)
        return false;
    bool retq = true;
    for(int k = 0; retq && (k < maxK_B); k++)
    {
        int l = index_0 + k;
        if(maxK_A <= l)
            l -= maxK_A;
        retq = (polygon_A.vertex(l) == polygon_B.vertex(k));
    }
    return retq;
}

bool PolygonEntry::visibility(Polygon_with_holes_2 const& polygon, std::vector<Point_2> const& points, std::vector<Polygon_2> & results)
{
    Arrangement_2 env;
    Arrangement_2 output;
    results.clear();
    //convert polygon to vector or list of segments
    std::list<Segment_2> segments;
    Polygon_2 outer_boundary = Polygon_2(polygon.outer_boundary());
    int maxK = outer_boundary.size();
    for(int k = 0; k < maxK; k++)
    {
        segments.push_back(outer_boundary.edge(k));
    }
    for(Polygon_with_holes_2::Hole_const_iterator it = polygon.holes_begin(); it != polygon.holes_end(); it++)
    {
        maxK = it->size();
        for(int k = 0; k < maxK; k++)
        {
            segments.push_back(it->edge(k));
        }
    }
    CGAL::insert_non_intersecting_curves(env, segments.begin(), segments.end());
    TEV tev(env);

    //find faceHandle: use the face whose outer boundary is the outer boundary of the polygon
    Face_const_handle faceHandle;
    Face_handle outputFaceHandle;
    bool found_face = false;

    for(Arrangement_2::Face_const_iterator it = env.faces_begin(); (!found_face) && (it != env.faces_end()); it++)
    {
        if(!it->is_unbounded())
            found_face = compareBoundary(outer_boundary, it->outer_ccb());
        if(found_face)
            faceHandle = it;
    }

    if(found_face)
    {
        maxK = points.size();
        for(int k = 0; k < maxK; k++)
        {
            outputFaceHandle = tev.compute_visibility(points[k], faceHandle, output);
            Arrangement_2::Ccb_halfedge_const_circulator boundary_start = outputFaceHandle->outer_ccb();
            Arrangement_2::Ccb_halfedge_const_circulator circ = boundary_start;
            Polygon_2 aux; aux.clear();
            do
            {
                aux.push_back(circ->target()->point());
                circ++;
            }while(circ != boundary_start);
            results.push_back(aux);
        }
    }

    return found_face;

#if 0
    std::list<Query_result> handles;
    CGAL::locate(env, points.begin(), points.end(), std::back_inserter(handles));

    maxK = points.size();
    std::list<Query_result>::const_iterator it = handles.begin();
    for(int k = 0; k < maxK; k++)
    {
        if(const Face_const_handle* faceHandle = boost::get<Face_const_handle>(&(it->second)))
        {
            //location in a face
            outputFaceHandle = tev.compute_visibility(points[k], faceHandle, output);
        }
        else if(const Halfedge_const_handle* edgeHandle = boost::get<Halfedge_const_handle>(&(it->second)))
        {
            //location on an edge
        }
        else
        {
            //Location in an isolated vertex: return just the point;
        }
        it++;
    }
#endif
}

}
