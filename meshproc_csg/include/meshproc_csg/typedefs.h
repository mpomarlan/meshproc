#ifndef __MESHPROC_CSG_TYPEDEFS__

#define __MESHPROC_CSG_TYPEDEFS__

namespace meshproc_csg
{

#define __USE_EXACT_KERNEL__
#ifdef __USE_EXACT_KERNEL__
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#endif
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_2 Point_2;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;

#if 0
typedef CGAL::Bounded_kernel<Kernel> Bounded_Kernel;
typedef CGAL::Nef_polyhedron_2<Bounded_Kernel> Nef_polygon;
typedef Nef_polygon::Point Point2D;
typedef Nef_polygon::Explorer Explorer;
typedef Explorer::Face_const_iterator Polygon_Face_const_iterator;
typedef Explorer::Hole_const_iterator Polygon_Hole_const_iterator;
typedef Explorer::Halfedge_around_face_const_circulator Polygon_Halfedge_around_face_const_circulator;
typedef Explorer::Vertex_const_handle Polygon_Vertex_const_handle;
#endif

}

#endif
