#ifndef __BUILD_MESH_H__

#define __BUILD_MESH_H__

namespace meshproc_csg
{
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

