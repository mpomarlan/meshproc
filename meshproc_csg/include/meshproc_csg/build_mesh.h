#ifndef __BUILD_MESH_H__

#define __BUILD_MESH_H__

namespace meshproc_csg
{

// A modifier to create a polyhedron from another, while switching between kernels
template <class HDSTarg, class PolyhedronSrc>
class Convert_mesh_kernel : public CGAL::Modifier_base<HDSTarg> {
public:
    Convert_mesh_kernel(PolyhedronSrc const& P):P(P) {}
    void operator()(HDSTarg& hds)
    {
        CGAL::Polyhedron_incremental_builder_3<HDSTarg> B(hds, true);
        B.begin_surface(P.size_of_vertices(), P.size_of_facets(), 0);
        typedef typename HDSTarg::Vertex Vertex;
        typedef typename Vertex::Point Point;
        for(typename PolyhedronSrc::Vertex_const_iterator it = P.vertices_begin(); it != P.vertices_end(); it++)
        {
            B.add_vertex( Point( ::CGAL::to_double(it->point().x()),
                                 ::CGAL::to_double(it->point().y()),
                                 ::CGAL::to_double(it->point().z())));
        }
        typedef CGAL::Inverse_index<typename PolyhedronSrc::Vertex_const_iterator> Index;
        Index index(P.vertices_begin(), P.vertices_end());
        for(typename PolyhedronSrc::Facet_const_iterator it = P.facets_begin();
            it != P.facets_end(); it++)
        {
            typename Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
            typename Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end;
            std::size_t n = circulator_size(hc_end);
            B.begin_facet();
            do
            {
                B.add_vertex_to_facet(index[hc_a->vertex()]);
                hc_a++;
            }while(hc_a != hc_end);
            B.end_facet();
        }
        B.end_surface();
    }

protected:
    PolyhedronSrc const& P;
};

// A modifier creating a polyhedron with the incremental builder from a triangulated mesh.
template <class HDS>
class Build_mesh : public CGAL::Modifier_base<HDS> {
public:
    Build_mesh(trimesh::TriMesh *M):M(M) {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

        B.begin_surface( M->vertices.size(), M->faces.size(), 0);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        int maxK = M->vertices.size();
        for(int k = 0; k < maxK; k++)
        {
            B.add_vertex( Point( M->vertices[k][0], M->vertices[k][1], M->vertices[k][2]));
        }
        maxK = M->faces.size();
        for(int k = 0; k < maxK; k++)
        {
            B.begin_facet();
            B.add_vertex_to_facet(M->faces[k].v[0]);
            B.add_vertex_to_facet(M->faces[k].v[1]);
            B.add_vertex_to_facet(M->faces[k].v[2]);
            B.end_facet();
        }
        B.end_surface();
    }
protected:
    trimesh::TriMesh *M;
};

void getNormal(Polyhedron::Face_const_handle const& f, double &x, double &y, double &z)
{
    std::cerr << "DBG getting normal " << std::endl;
    Polyhedron::Halfedge_around_facet_const_circulator hc_a = f->facet_begin();
    std::cerr << "DBG circs " << std::endl;
    Polyhedron::Halfedge_around_facet_const_circulator hc_b = hc_a; hc_b++;
    std::cerr << "DBG circs " << std::endl;
    Polyhedron::Halfedge_around_facet_const_circulator hc_c = hc_b; hc_c++;
    std::cerr << "DBG circs " << std::endl;

    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    x1 = ::CGAL::to_double(hc_a->vertex()->point().x());
    y1 = ::CGAL::to_double(hc_a->vertex()->point().y());
    z1 = ::CGAL::to_double(hc_a->vertex()->point().z());
    x2 = ::CGAL::to_double(hc_b->vertex()->point().x());
    y2 = ::CGAL::to_double(hc_b->vertex()->point().y());
    z2 = ::CGAL::to_double(hc_b->vertex()->point().z());
    x3 = ::CGAL::to_double(hc_c->vertex()->point().x());
    y3 = ::CGAL::to_double(hc_c->vertex()->point().y());
    z3 = ::CGAL::to_double(hc_c->vertex()->point().z());

    double vxX, vxY, vyX, vyY, vzX, vzY;

    vxX = x1-x2;
    vyX = y1-y2;
    vzX = z1-z2;

    vxY = x3-x2;
    vyY = y3-y2;
    vzY = z3-z2;

    x = -(vyX*vzY) + (vzX*vyY);
    y = -(vzX*vxY) + (vxX*vzY);
    z = -(vxX*vyY) + (vyX*vxY);

    double l = std::sqrt(x*x + y*y + z*z);
    std::cerr << "DBG normal " << x << " " << y << " " << z << " " << l << std::endl;
    x /= l;
    y /= l;
    z /= l;
}

template <class HDS>
class Build_solidification : public CGAL::Modifier_base<HDS> {
public:
    Build_solidification(Polyhedron const& mesh, double thickness): mesh(mesh), thickness(thickness) {}
    void operator()(HDS& hds)
    {
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        B.begin_surface(2*mesh.size_of_vertices(), 2*mesh.size_of_facets(), 0);

        typedef CGAL::Inverse_index<Polyhedron::Vertex_const_iterator> Index;
        Index index( mesh.vertices_begin(), mesh.vertices_end());

        std::vector<double> nxV, nyV, nzV; nxV.clear(); nyV.clear(); nzV.clear();
        nxV.reserve(mesh.size_of_vertices()); nyV.reserve(mesh.size_of_vertices()); nzV.reserve(mesh.size_of_vertices());

        for(Polyhedron::Vertex_const_iterator it = mesh.vertices_begin();
            it != mesh.vertices_end(); it++)
        {
            double x, y, z;
            x = ::CGAL::to_double(it->point().x());
            y = ::CGAL::to_double(it->point().y());
            z = ::CGAL::to_double(it->point().z());
            B.add_vertex(Point(x, y, z));
            double nx, ny, nz;
            nx = ny = nz = 0.0;
            Polyhedron::Vertex::Halfedge_around_vertex_const_circulator vit = it->vertex_begin();
            do
            {
                double xp, yp, zp, l;
                if(!vit->is_border())
                {
                    getNormal(vit->facet(), xp, yp, zp);
                    nx += (xp);
                    ny += (yp);
                    nz += (zp);
                }
                vit++;
            }while(vit != it->vertex_begin());
            double l = std::sqrt(nx*nx + ny*ny + nz*nz);
            std::cerr << "DBG nrsum " << nx << " " << ny << " " << nz << " " << l << std::endl;
            nx /= l;
            ny /= l;
            nz /= l;
            nxV.push_back(nx);
            nyV.push_back(ny);
            nzV.push_back(nz);
        }

        int k = 0;
        int maxK = nxV.size();
        for(Polyhedron::Vertex_const_iterator it = mesh.vertices_begin();
            it != mesh.vertices_end(); it++)
        {
            double x, y, z;
            x = ::CGAL::to_double(it->point().x());
            y = ::CGAL::to_double(it->point().y());
            z = ::CGAL::to_double(it->point().z());
            double sx, sy, sz;
            sx = nxV[k]*thickness;
            sy = nyV[k]*thickness;
            sz = nzV[k]*thickness;
            B.add_vertex(Point(x + sx, y + sy, z + sz));
            k++;
        }

        for(Polyhedron::Facet_const_iterator it = mesh.facets_begin();
            it != mesh.facets_end(); it++)
        {
            Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
            Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end;
            B.begin_facet();
            do
            {
                B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->vertex())]); hc_a++;
            }while(hc_a != hc_end);
            B.end_facet();
            B.begin_facet();
            do
            {
                B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->vertex())] + maxK); hc_a--;
            }while(hc_a != hc_end);
            B.end_facet();

            do
            {
                if(hc_a->opposite()->is_border())
                {
                    B.begin_facet();
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->vertex())]);
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->prev()->vertex())]);
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->vertex())] + maxK);
                    B.end_facet();

                    B.begin_facet();
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->prev()->vertex())]);
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->prev()->vertex())] + maxK);
                    B.add_vertex_to_facet(index[Polyhedron::Vertex_const_iterator(hc_a->vertex())] + maxK);
                    B.end_facet();
                }
                hc_a++;
            }while(hc_a != hc_end);
        }

        B.end_surface();
    }

protected:
    Polyhedron const& mesh;
    double thickness;
};

template <class HDS>
class Build_from_triangulation : public CGAL::Modifier_base<HDS> {
public:
    Build_from_triangulation(CDT & cdt, Eigen::Affine3d const& R, bool fillHoles): cdt(cdt), R(R), fillHoles(fillHoles) {}
    void operator()(HDS& hds)
    {
        CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

        int maxK = 0;
        for(CDT::Finite_faces_iterator it = cdt.finite_faces_begin();
            it != cdt.finite_faces_end(); it++)
            if(it->info().in_domain())
                maxK++;

        B.begin_surface(cdt.number_of_vertices(), maxK, 0);

        std::map<CDT::Vertex_handle, int> vhMap;

        int k = 0;
        for(CDT::Finite_vertices_iterator it = cdt.finite_vertices_begin();
            it != cdt.finite_vertices_end(); it++)
        {
            double x, y;
            double xR, yR, zR;
            x = ::CGAL::to_double(it->point().x());
            y = ::CGAL::to_double(it->point().y());
            xR = R(0,0)*x + R(0,1)*y;
            yR = R(1,0)*x + R(1,1)*y;
            zR = R(2,0)*x + R(2,1)*y;
            vhMap.insert(std::pair<CDT::Vertex_handle, int>(it, k));
            B.add_vertex(Point(xR, yR, zR));
            k++;
        }

        for(CDT::Finite_faces_iterator it = cdt.finite_faces_begin();
            it != cdt.finite_faces_end(); it++)
        {
            if(fillHoles || (it->info().in_domain()))
            {
                B.begin_facet();
                B.add_vertex_to_facet(vhMap.find(it->vertex(0))->second);
                B.add_vertex_to_facet(vhMap.find(it->vertex(1))->second);
                B.add_vertex_to_facet(vhMap.find(it->vertex(2))->second);
                B.end_facet();
            }
        }
        B.end_surface();
    }

protected:
    CDT & cdt;
    Eigen::Affine3d const& R;
    bool fillHoles;
};

template <class HDS>
class Build_polygon : public CGAL::Modifier_base<HDS> {
public:
    Build_polygon(std::vector<double> const& xP, std::vector<double> const& yP,
                  std::vector<double> const& zP, Eigen::Affine3d const& R): xV(xP), yV(yP), zV(zP), R(R) {}
    void operator()(HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

    int maxK = xV.size();
    B.begin_surface(maxK, 1, 0);
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    B.begin_facet();
    for(int k = 0; k < maxK; k++)
    {
        double x, y, z;
        double xP, yP;
        double xR, yR, zR;
        x = xV[k];
        y = yV[k];
        z = zV[k];
        xP = R(0, 0)*x + R(1, 0)*y + R(2, 0)*z;
        yP = R(0, 1)*x + R(1, 1)*y + R(2, 1)*z;
        xR = R(0, 0)*xP + R(0, 1)*yP;
        yR = R(1, 0)*xP + R(1, 1)*yP;
        zR = R(2, 0)*xP + R(2, 1)*yP;
        B.add_vertex(Point(xR, yR, zR));
        B.add_vertex_to_facet(k);
    }
    B.end_facet();
    B.end_surface();
    }
protected:
    std::vector<double> const& xV;
    std::vector<double> const& yV;
    std::vector<double> const& zV;
    Eigen::Affine3d const& R;
};

template <class HDS>
class Build_box : public CGAL::Modifier_base<HDS> {
public:
    Build_box(double maxX, double minX, double maxY, double minY, double maxZ, double minZ):
        maxX(maxX), minX(minX), maxY(maxY), minY(minY), maxZ(maxZ), minZ(minZ)
    {}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

        B.begin_surface( 8, 6, 0);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        B.add_vertex( Point( maxX, maxY, maxZ));
        B.add_vertex( Point( maxX, maxY, minZ));
        B.add_vertex( Point( maxX, minY, maxZ));
        B.add_vertex( Point( maxX, minY, minZ));
        B.add_vertex( Point( minX, maxY, maxZ));
        B.add_vertex( Point( minX, maxY, minZ));
        B.add_vertex( Point( minX, minY, maxZ));
        B.add_vertex( Point( minX, minY, minZ));
        B.begin_facet();
        B.add_vertex_to_facet(0);
        B.add_vertex_to_facet(1);
        B.add_vertex_to_facet(3);
        B.add_vertex_to_facet(2);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(4);
        B.add_vertex_to_facet(6);
        B.add_vertex_to_facet(7);
        B.add_vertex_to_facet(5);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(1);
        B.add_vertex_to_facet(5);
        B.add_vertex_to_facet(7);
        B.add_vertex_to_facet(3);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(2);
        B.add_vertex_to_facet(3);
        B.add_vertex_to_facet(7);
        B.add_vertex_to_facet(6);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(4);
        B.add_vertex_to_facet(0);
        B.add_vertex_to_facet(2);
        B.add_vertex_to_facet(6);
        B.end_facet();
        B.begin_facet();
        B.add_vertex_to_facet(5);
        B.add_vertex_to_facet(1);
        B.add_vertex_to_facet(0);
        B.add_vertex_to_facet(4);
        B.end_facet();
        B.end_surface();
    }
protected:
    double maxX, minX, maxY, minY, maxZ, minZ;
};

}

#endif

