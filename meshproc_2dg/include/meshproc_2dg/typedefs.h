#ifndef __MESHPROC_2DG_TYPEDEFS__

#define __MESHPROC_2DG_TYPEDEFS__

namespace meshproc_2dg
{
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Partition_traits_2<Kernel>                  Traits;
typedef CGAL::Is_convex_2<Traits>                         Is_convex_2;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef Polygon_2::Point_2                                Point_2;
typedef Polygon_2::Vertex_const_iterator                  Vertex_const_iterator;
typedef std::list<Polygon_2>                              Polygon_list;
typedef std::vector<Polygon_2>                            Polygon_vector;
typedef CGAL::Partition_is_valid_traits_2<Traits, Is_convex_2>
                                                          Validity_traits;
typedef CGAL::Straight_skeleton_2<Kernel>                 Skeleton;
typedef boost::shared_ptr<Skeleton>                       SkeletonPtr;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

}

#endif
