#ifndef __TRIANGULATE_H__

#define __TRIANGULATE_H__

/*A copy-paste of CGAL's triangulate_polyhedron.h, necessary to support triangulation
of the Exact_predicates_exact_constructions kernel: the call to compute_facet_normal
now calls a custom version, which uses ::CGAL::to_double.*/

namespace CGAL {

template <class Facet, class Kernel>
typename Kernel::Vector_3 compute_exact_facet_normal_internal(const Facet& f, double &len)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Facet::Halfedge_around_facet_const_circulator HF_circulator;
  Vector normal = CGAL::NULL_VECTOR;
  HF_circulator he = f.facet_begin();
  HF_circulator end = he;
  CGAL_For_all(he,end)
  {
    const Point& prev = he->prev()->vertex()->point();
    const Point& curr = he->vertex()->point();
    const Point& next = he->next()->vertex()->point();
    Vector n = CGAL::cross_product(next-curr,prev-curr);
    normal = normal + n;
  }

  double aux = std::sqrt(::CGAL::to_double(normal * normal));
  if(aux < 0.0001)
  {
      len  = 0.0;
      return CGAL::NULL_VECTOR;
  }
  len = 1.0;
  return normal / aux;
}


template <class Facet, class Kernel>
typename Kernel::Vector_3 compute_exact_facet_normal(const Facet& f)
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Facet::Halfedge_around_facet_const_circulator HF_circulator;
  double len;
  Vector normal = compute_exact_facet_normal_internal<Facet, Kernel>(f, len);
  if(len < 0.0001)
  {
      std::cerr << "WARNING: degenerate face detected, attempting to compute normal based on neighbours." << std::endl;
      HF_circulator he = f.facet_begin();
      HF_circulator end = he;
      CGAL_For_all(he, end)
      {
          std::cerr << "    neighbour:" << std::endl;
          double lenOpp = 0.0;
          Vector normalOpp = compute_exact_facet_normal_internal<Facet, Kernel>(*(he->opposite()->facet()), lenOpp);
          std::cerr << "    " << lenOpp << ": "
                    << ::CGAL::to_double(normalOpp.x())
                    << " "  << ::CGAL::to_double(normalOpp.y())
                    << " "  << ::CGAL::to_double(normalOpp.z())
                    << std::endl;
          normal = normal + normalOpp;
      }
      len = std::sqrt(::CGAL::to_double(normal * normal));
      if(len < 0.0001)
      {
          std::cerr << "    WARNING: spike detected." << std::endl;
          return Vector(1, 0, 0);
      }
      std::cerr << "    Returning a neighbour-based normal of len " << len
                << ": " << ::CGAL::to_double(normal.x())
                << " " << ::CGAL::to_double(normal.y())
                << " " << ::CGAL::to_double(normal.z())
                << std::endl;
      return normal/len;
  }
  return normal;
}

template <class Polyhedron>
class Triangulate_modifier_exact
  : public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
{
  typedef typename Polyhedron::HalfedgeDS HDS;
  typedef typename Polyhedron::Traits Traits;

  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet Facet;
  typedef typename Polyhedron::Facet_iterator Facet_iterator;
  typedef typename Polyhedron::Facet_handle Facet_handle;

  typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits>   P_traits;

  typedef CGAL::Triangulation_vertex_base_with_info_2<Halfedge_handle,
                                                      P_traits>        Vb;

  struct Face_info {
    typename Polyhedron::Halfedge_handle e[3];
    bool is_external;
  };

  typedef CGAL::Triangulation_face_base_with_info_2<Face_info,
                                                    P_traits>          Fb1;

  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>   Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                  TDS;
  typedef CGAL::No_intersection_tag                                    Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits,
                                                     TDS,
                                                     Itag>             CDTbase;
  typedef CGAL::Constrained_triangulation_plus_2<CDTbase>              CDT;

public:
  Triangulate_modifier_exact()
  {
  }

  bool is_external(typename CDT::Face_handle fh) const {
    return fh->info().is_external;
  }

  void operator()(HDS& hds) {
    CGAL::HalfedgeDS_decorator<HDS> decorator(hds);
    typedef typename HDS::Halfedge Halfedge;

    // One need to store facet handles into a vector, because the list of
    // facets of the polyhedron will be modified during the loop, and
    // that invalidates the range [facets_begin(), facets_end()[.
    std::vector<Facet_handle> facets;
    facets.reserve(hds.size_of_faces());
    for(Facet_iterator
          fit = hds.faces_begin(),
          end = hds.faces_end();
        fit != end; ++fit) {
      facets.push_back(fit);
    }

    // Iterates on the vector of facet handles
    for(typename std::vector<Facet_handle>::iterator
          fit_it = facets.begin(),
          end = facets.end();
        fit_it != end; ++fit_it)
    {
      Facet_handle fit = *fit_it;
      typename Traits::Vector_3 normal =
        compute_exact_facet_normal<Facet,Traits>(*fit);

      P_traits cdt_traits(normal);
      CDT cdt(cdt_traits);

      typename Facet::Halfedge_around_facet_circulator
        he_circ = fit->facet_begin(),
        he_circ_end(he_circ);
      typename CDT::Vertex_handle previous, first;
      do {
        typename CDT::Vertex_handle vh = cdt.insert(he_circ->vertex()->point());
        if(first == 0) {
          first = vh;
        }
        vh->info() = he_circ;
        if(previous != 0 && previous != vh) {
          cdt.insert_constraint(previous, vh);
        }
        previous = vh;
      } while( ++he_circ != he_circ_end );
      cdt.insert_constraint(previous, first);

      // sets mark is_external
      for(typename CDT::All_faces_iterator
            fit = cdt.all_faces_begin(),
            end = cdt.all_faces_end();
          fit != end; ++fit)
      {
        fit->info().is_external = false;
      }
      std::queue<typename CDT::Face_handle> face_queue;
      face_queue.push(cdt.infinite_vertex()->face());
      while(! face_queue.empty() ) {
        typename CDT::Face_handle fh = face_queue.front();
        face_queue.pop();
        if(fh->info().is_external) continue;
        fh->info().is_external = true;
        for(int i = 0; i <3; ++i) {
          if(!cdt.is_constrained(std::make_pair(fh, i)))
          {
            face_queue.push(fh->neighbor(i));
          }
        }
      }

      // then modify the polyhedron
      decorator.make_hole(fit->halfedge());
      for(typename CDT::Finite_edges_iterator
            eit = cdt.finite_edges_begin(),
            end = cdt.finite_edges_end();
          eit != end; ++eit)
      {
        typename CDT::Face_handle fh = eit->first;
        const int index = eit->second;
        typename CDT::Face_handle opposite_fh = fh->neighbor(eit->second);

        if((fh == 0) || (opposite_fh == 0))
            std::cerr << "ERROR: either fh is zero (" << (fh == 0) << ") or opposite_fh is zero (" << (opposite_fh == 0) << ")" << std::endl;

        const int opposite_index = opposite_fh->index(fh);
        const typename CDT::Vertex_handle va = fh->vertex(cdt. cw(index));
        const typename CDT::Vertex_handle vb = fh->vertex(cdt.ccw(index));

        if( ! (is_external(fh) && is_external(opposite_fh)) &&
            ! cdt.is_constrained(*eit) )
        {
          // strictly internal edge
          Halfedge_handle h = hds.edges_push_back(Halfedge(),
                                                  Halfedge());
          fh->info().e[index] = h;
          opposite_fh->info().e[opposite_index] = h->opposite();

          decorator.set_vertex(h, va->info()->vertex());
          decorator.set_vertex(h->opposite(), vb->info()->vertex());
        }
        if( cdt.is_constrained(*eit) )
        {
          if(!is_external(fh)) {
            fh->info().e[index] = va->info();
          }
          if(!is_external(opposite_fh)) {
            opposite_fh->info().e[opposite_index] = vb->info();
          }
        }
      }
      for(typename CDT::Finite_faces_iterator
            fit = cdt.finite_faces_begin(),
            end = cdt.finite_faces_end();
          fit != end; ++fit)
      {
        if(!is_external(fit))
        {
          Halfedge_handle h0 = fit->info().e[0];
          Halfedge_handle h1 = fit->info().e[1];
          Halfedge_handle h2 = fit->info().e[2];
          CGAL_assertion( h0 != Halfedge_handle() );
          CGAL_assertion( h1 != Halfedge_handle() );
          CGAL_assertion( h2 != Halfedge_handle() );

          typedef typename Halfedge::Base HBase;
          h0->HBase::set_next(h1);
          decorator.set_prev(h1, h0);
          h1->HBase::set_next(h2);
          decorator.set_prev(h2, h1);
          h2->HBase::set_next(h0);
          decorator.set_prev(h0, h2);

          decorator.fill_hole(h0);
        }
      }
    } // end loop on facets of the input polyhedron
  }
}; // end class Triangulate_modifier_exact
}

#endif

