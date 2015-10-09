#ifndef __MESHPROC_P2M_TYPEDEFS__

#define __MESHPROC_P2M_TYPEDEFS__

namespace meshproc_p2m
{
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef CGAL::Scale_space_surface_reconstruction_3<Kernel>      Reconstruction;
typedef Reconstruction::Point                                   Point;
typedef std::vector< Point >                                    Point_collection;
typedef Reconstruction::Triple_const_iterator                   Triple_iterator;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

}

#endif
