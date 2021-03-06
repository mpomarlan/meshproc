#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/bounding_box.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Bounded_kernel.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>


#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <boost/foreach.hpp>

#include <trimesh/TriMesh.h>
#include <trimesh/TriMesh_algo.h>

#include <meshproc_csg/kdtree++/kdtree.hpp>

#include <shape_msgs/Mesh.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <meshproc_csg/typedefs.h>

namespace meshproc_csg{
struct FaceInfo2
{
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){
    return (nesting_level%2 == 1) || (!nesting_level);
  }
};

typedef CGAL::Triangulation_vertex_base_2<Kernel>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Kernel>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    TDS;
typedef CGAL::Exact_predicates_tag                                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>  CDT;
}

#include <meshproc_csg/build_mesh.h>
#include <meshproc_csg/triangulate.h>

#include <meshproc_csg/csg.h>

namespace meshproc_csg
{

typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes;
typedef CGAL::Polygon_2<Kernel> Polygon;

typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron> Traits;
// default property maps
typedef boost::property_map<Polyhedron,
                            boost::vertex_external_index_t>::type  Vertex_index_map;
typedef boost::property_map<Polyhedron,
                            CGAL::halfedge_external_index_t>::type Halfedge_index_map;
typedef boost::property_map<Polyhedron,
                            CGAL::face_external_index_t>::type     Face_index_map;
typedef CGAL::Surface_mesh_shortest_path<Traits,
                                         Vertex_index_map,
                                         Halfedge_index_map,
                                         Face_index_map>  Surface_mesh_shortest_path;
typedef boost::graph_traits<Polyhedron> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::halfedge_iterator halfedge_iterator;
typedef Graph_traits::face_iterator face_iterator;


/******************************************************************************************/
/*** The code here is taken from the CGAL example Triangulation_2/polygon_triangulation ***/
/******************************************************************************************/

void mark_domains(CDT& ct,
             CDT::Face_handle start,
             int index,
             std::list<CDT::Edge>& border )
{
  if(start->info().nesting_level != -1){
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty()){
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1){
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++){
        CDT::Edge e(fh,i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1){
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(CDT& cdt)
{
  for(CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it){
    it->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while(! border.empty()){
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1){
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}
/******************************************************************************************/
/** The code above was taken from the CGAL example Triangulation_2/polygon_triangulation **/
/******************************************************************************************/

double duplicate_threshold = 0.00001; /*Square of a distance. If vertices are closer than this, they are
                                        treated as the same vertex by remove_duplicates.*/

/*Psst, this function doesn't exist :P

Not used anymore. Way too slow on usual-size meshes.*/
void remove_duplicates(trimesh::TriMesh *mesh)
{
  std::vector<int> replacements;
  replacements.clear();
  replacements.push_back(0);
  int maxK = mesh->vertices.size();
  for(int k = 1; k < maxK; k++)
  {
    int maxJ = k;
    bool replaced = false;
    for(int j = 0; j < maxJ; j++)
      if((j == replacements[j]) && (trimesh::dist2(mesh->vertices[k], mesh->vertices[j]) < duplicate_threshold))
      {
        replacements.push_back(j);
        replaced = true;
      }
    if(!replaced)
      replacements.push_back(k);
  }

  maxK = mesh->faces.size();
  for(int k = 0; k< maxK; k++)
  {
    mesh->faces[k].v[0] = replacements[mesh->faces[k].v[0]];
    mesh->faces[k].v[1] = replacements[mesh->faces[k].v[1]];
    mesh->faces[k].v[2] = replacements[mesh->faces[k].v[2]];
  }
  trimesh::remove_unused_vertices(mesh);
}

MeshEntry::MeshEntry()
{
    props_updated = false;
    triangulated = false;
}

MeshEntry::MeshEntry(Polyhedron const& polyhedron)
{
    props_updated = false;
    triangulated = false;
    mesh_data = polyhedron;
}

MeshEntry::~MeshEntry()
{
}

void MeshEntry::clear(void)
{
    props_updated = false;
    triangulated = false;
    mesh_data.clear();
    nef_polyhedron.clear();
}
bool MeshEntry::loadFromFile(std::string const& filename, double duplicate_dist)
{
    trimesh::TriMesh *M;
    bool fileAccessible = false;
    {
        std::ifstream f(filename.c_str());
        fileAccessible = f.good();
    }
    if(!fileAccessible)
        return false;
    M = trimesh::TriMesh::read(filename.c_str());
    return loadFromTrimesh(M, duplicate_dist);
}
bool MeshEntry::loadFromMsg(shape_msgs::Mesh const& message, double duplicate_dist)
{
    trimesh::TriMesh M;
    M.vertices.clear(); M.vertices.reserve(message.vertices.size());
    M.faces.clear(); M.faces.reserve(message.triangles.size());
    int maxK = message.vertices.size();
    for(int k = 0; k < maxK; k++)
        M.vertices.push_back(trimesh::point(message.vertices[k].x,
                                            message.vertices[k].y,
                                            message.vertices[k].z));
    maxK = message.triangles.size();
    for(int k = 0; k < maxK; k++)
        M.faces.push_back(trimesh::TriMesh::Face(message.triangles[k].vertex_indices[0],
                                                 message.triangles[k].vertex_indices[1],
                                                 message.triangles[k].vertex_indices[2]));
    return loadFromTrimesh(&M, duplicate_dist);
}
bool MeshEntry::loadFromTrimesh(trimesh::TriMesh *M, double duplicate_dist)
{
    clear();
    MeshEntry::remove_duplicates(M, duplicate_dist);
    Build_mesh<Polyhedron::HalfedgeDS> build_mesh(M);
    mesh_data.delegate(build_mesh);
    nef_polyhedron = Nef_polyhedron(mesh_data);
    return true;
}

bool MeshEntry::setFromMeshEntry(MeshEntry const& A)
{
    name = "";
    mesh_data = A.mesh_data;
    nef_polyhedron = A.nef_polyhedron;
    props_updated = A.props_updated;
    triangulated = A.triangulated;
    is_closed = A.is_closed;
    is_manifold = A.is_manifold;
    is_orientable = A.is_orientable;
    is_csg_safe = A.is_csg_safe;
    nr_vertices = A.nr_vertices;
    nr_edges = A.nr_edges;
    nr_faces = A.nr_faces;
    nr_connected_components = A.nr_connected_components;
    Euler_characteristic = A.Euler_characteristic;
    volume = A.volume;
    return true;
}

bool MeshEntry::setFromUnion(MeshEntry const& A, MeshEntry const& B)
{
    nef_polyhedron = A.nef_polyhedron + B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromIntersection(MeshEntry const& A, MeshEntry const& B)
{
    nef_polyhedron = A.nef_polyhedron * B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromDifference(MeshEntry const& A, MeshEntry const& B)
{
    nef_polyhedron = A.nef_polyhedron - B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromSymmetricDifference(MeshEntry const& A, MeshEntry const& B)
{
    nef_polyhedron = A.nef_polyhedron ^ B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromMinkowskiSum(MeshEntry const& A, MeshEntry const& B)
{
    Nef_polyhedron cA = A.nef_polyhedron;
    Nef_polyhedron cB = B.nef_polyhedron;
    std::cout << "Copied nefs ..." << std::endl;
    nef_polyhedron = CGAL::minkowski_sum_3(cA, cB);
    return update_mesh();
}
bool MeshEntry::setFromMinkowskiErosion(MeshEntry const& A, MeshEntry const& B)
{
    Polyhedron bbox = A.getBoundingBox(2);
    Nef_polyhedron nefRA(bbox);
    Nef_polyhedron cB = B.nef_polyhedron;
    nefRA = nefRA - A.nef_polyhedron;
    nefRA = CGAL::minkowski_sum_3(nefRA, cB);
    nef_polyhedron = A.nef_polyhedron - nefRA;
    return update_mesh();
}
bool MeshEntry::setFromSelectComponent(MeshEntry const& A, meshproc_csg::Point const& P)
{
    return true;
}
bool MeshEntry::setFromConvexHull(MeshEntry & A)
{
    typedef Polyhedron::Point_const_iterator Point_const_iterator;
    //typedef typename std::iterator_traits<Point_const_iterator>::value_type Point_3;
    typedef typename CGAL::internal::Convex_hull_3::Default_traits_for_Chull_3<Point>::type Traits;
    Traits traits;
    typename Traits::Collinear_3 collinear = traits.collinear_3_object();
    Polyhedron P;
    bool notAllCollinear = false;
    if(2 < A.mesh_data.size_of_vertices())
    {
        Point_const_iterator p1, p2, p3;
        p1 = A.mesh_data.points_begin();
        p2 = p1; p2++;
        p3 = p2; p3++;
        for(; (p3 != A.mesh_data.points_end()) && (!notAllCollinear); p3++)
        {
            notAllCollinear = !collinear(*p1, *p2, *p3);
        }
    }
    if((3 < A.mesh_data.size_of_vertices()) && (notAllCollinear))
    {
        CGAL::convex_hull_3(A.mesh_data.points_begin(), A.mesh_data.points_end(), P);
        nef_polyhedron = Nef_polyhedron(P);
    }
    else
        nef_polyhedron = Nef_polyhedron(A.mesh_data);
    nef_polyhedron = nef_polyhedron.regularization();
    //std::string refName(dummyObject.type().name); // "N4CGAL12Polyhedron_3INS_5EpeckENS_18Polyhedron_items_3ENS_18HalfedgeDS_defaultESaIiEEE"
    return update_mesh();
}

bool MeshEntry::setFromMesh2Prism(MeshEntry const& A, double height, double depth, bool filter, double zcomp)
{
    Nef_polyhedron prism;
    for(Polyhedron::Facet_const_iterator it = A.mesh_data.facets_begin();
        it!= A.mesh_data.facets_end(); it++)
    {
        Polyhedron facetPrism;
        Build_prism<Polyhedron::HalfedgeDS> prismBuilder(it, height, depth, filter, zcomp);
        facetPrism.delegate(prismBuilder);
        Nef_polyhedron aux = Nef_polyhedron(facetPrism);
        if(!aux.is_empty())
            prism = prism + aux;
    }
    nef_polyhedron = prism;
    return update_mesh();
}

bool MeshEntry::filterByNormal(std::string const& filename, double nx, double ny, double nz, double toleratedAngle)
{
    triangulate_mesh();
    double zcomp = std::cos(toleratedAngle);
    trimesh::TriMesh M;
    M.vertices.clear();
    M.faces.clear();
    int vIndex = 0;
    for(Polyhedron::Facet_const_iterator it = mesh_data.facets_begin(); it != mesh_data.facets_end(); it++)
    {
        double fnx, fny, fnz;
        std::vector<double> fx, fy, fz;

        Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
        Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end;
        Polyhedron::Halfedge_around_facet_const_circulator hc_b = hc_a; hc_b++;

        fnx = 0; fny = 0; fnz = 0;
        do
        {
            double ax = ::CGAL::to_double(hc_a->vertex()->point().x());
            double ay = ::CGAL::to_double(hc_a->vertex()->point().y());
            double az = ::CGAL::to_double(hc_a->vertex()->point().z());
            fx.push_back(ax);
            fy.push_back(ay);
            fz.push_back(az);
            double bx = ::CGAL::to_double(hc_b->vertex()->point().x());
            double by = ::CGAL::to_double(hc_b->vertex()->point().y());
            double bz = ::CGAL::to_double(hc_b->vertex()->point().z());
            fnx += (ay - by)*(az + bz);
            fny += (az - bz)*(ax + bx);
            fnz += (ax - bx)*(ay + by);
            hc_a++;
            hc_b++;
        }while(hc_a != hc_end);

        double l = std::sqrt(fnx*fnx + fny*fny + fnz*fnz);
        if(l < 0.000001)
        {
            fnx = 0;
            fny = 0;
            fnz = 0;
        }
        else
        {
            fnx /= l;
            fny /= l;
            fnz /= l;
        }

        int maxK = fx.size();

        if(zcomp <= (fnx*nx + fny*ny + fnz*nz))
        {
            trimesh::point aux;
            for(int k = 0; k < maxK; k++)
            {
                aux[0] = fx[k]; aux[1] = fy[k]; aux[2] = fz[k];
                M.vertices.push_back(aux);
            }
            trimesh::TriMesh::Face auxF;
            auxF.v[0] = vIndex;
            auxF.v[1] = vIndex + 1;
            auxF.v[2] = vIndex + 2;
            M.faces.push_back(auxF);
            vIndex += maxK;
        }
    }
    bool wrote_file = M.write(filename.c_str());

    return wrote_file;
}

bool MeshEntry::getBoundaryPolygon(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                        std::vector<int> &edge_L, std::vector<int> &edge_R)
{
    mesh_data.normalize_border();
    typedef Polyhedron::Edge_const_iterator Edge_const_iterator;
    typedef Polyhedron::Vertex_const_handle Vertex_const_handle;
    typedef std::map<Vertex_const_handle, int> Vertex_const_handle_map;
    Edge_const_iterator it = mesh_data.border_edges_begin();
    Vertex_const_handle_map knownVertices;
    int maxK = mesh_data.size_of_border_edges();
    knownVertices.clear();
    x.clear(); x.resize(maxK);
    y.clear(); y.resize(maxK);
    z.clear(); z.resize(maxK);
    edge_L.clear(); edge_L.reserve(maxK);
    edge_R.clear(); edge_R.reserve(maxK);
    for(; it != mesh_data.edges_end(); it++)
    {
        Vertex_const_handle v1, v2;
        std::pair<Vertex_const_handle_map::iterator, bool> insOp;
        double x1, y1, z1;
        double x2, y2, z2;
        v1 = it->vertex();
        v2 = it->opposite()->vertex();
        int k1, k2;
        insOp = knownVertices.insert(std::pair<Vertex_const_handle, int>(v1, knownVertices.size()));
        k1 = insOp.first->second;
        insOp = knownVertices.insert(std::pair<Vertex_const_handle, int>(v2, knownVertices.size()));
        k2 = insOp.first->second;
        x1 = ::CGAL::to_double(v1->point().x());
        x2 = ::CGAL::to_double(v2->point().x());
        y1 = ::CGAL::to_double(v1->point().y());
        y2 = ::CGAL::to_double(v2->point().y());
        z1 = ::CGAL::to_double(v1->point().z());
        z2 = ::CGAL::to_double(v2->point().z());
        x[k1] = x1; x[k2] = x2;
        y[k1] = y1; y[k2] = y2;
        z[k1] = z1; z[k2] = z2;
        edge_L.push_back(k1);
        edge_R.push_back(k2);
    }
    return true;
}

void getDSeg(double x, double y, double xn, double yn, double f, double &dx, double &dy)
{
    dx = xn - x;
    dy = yn - y;
    double l = std::sqrt(dx*dx + dy*dy);
    if(l < 0.01*f)
    {
        return;
    }
    dx /= l;
    dy /= l;
    dx *= (-f);
    dy *= (-f);
}

void inflatePolygon(Polygon &poly, double f)
{
    int maxK = poly.size();
    int prev = maxK - 1;
    int next = 1;
    std::vector<double> xv, yv;
    xv.resize(maxK);
    yv.resize(maxK);
    for(int k = 0; k < maxK; k++)
    {
        double x, y;
        double xp, yp;
        double xn, yn;
        x = ::CGAL::to_double(poly.vertex(k).x());
        y = ::CGAL::to_double(poly.vertex(k).y());
        xn = ::CGAL::to_double(poly.vertex(next).x());
        yn = ::CGAL::to_double(poly.vertex(next).y());
        xp = ::CGAL::to_double(poly.vertex(prev).x());
        yp = ::CGAL::to_double(poly.vertex(prev).y());
        double dx1, dy1, dz1;
        double dx2, dy2, dz2;
        getDSeg(x, y, xn, yn, f, dx1, dy1);
        getDSeg(x, y, xp, yp, f, dx2, dy2);
        x = x + dx1 + dx2;
        y = y + dy1 + dy2;
        xv[k] = x;
        yv[k] = y;
        prev++;
        next++;
        if(prev == maxK)
            prev = 0;
        if(next == maxK)
            next = 0;
    }
    poly.clear();
    for(int k = 0; k < maxK; k++)
        poly.push_back(Point_2(xv[k], yv[k]));
}

bool MeshEntry::setFromProjection(MeshEntry const& A, double nx, double ny, double nz, bool fillHoles)
{
    int dbgK = 0;
    double n = std::sqrt(nx*nx + ny*ny + nz*nz);
    if(0.000001 > n)
    {
        return false;
    }
    nx /= n;
    ny /= n;
    nz /= n;
    double ns = std::sqrt(nx*nx + ny*ny);
    Eigen::Affine3d afftran;
    afftran(0,0) = 1.0; afftran(0,1) = 0.0; afftran(0,2) = 0.0; afftran(0,3) = 0.0;
    afftran(1,0) = 0.0; afftran(1,1) = 1.0; afftran(1,2) = 0.0; afftran(1,3) = 0.0;
    afftran(2,0) = 0.0; afftran(2,1) = 0.0; afftran(2,2) = 1.0; afftran(2,3) = 0.0;
    if(0.000001 < ns)
    {
        double nxa = nx/ns;
        double nya = ny/ns;
        afftran(0,0) =  nya; afftran(0,1) =          nz*nxa; afftran(0,2) = nx; afftran(0,3) = 0.0;
        afftran(1,0) = -nxa; afftran(1,1) =          nz*nya; afftran(1,2) = ny; afftran(1,3) = 0.0;
        afftran(2,0) =  0.0; afftran(2,1) = -nx*nxa -ny*nya; afftran(2,2) = nz; afftran(2,3) = 0.0;
    }

    Polygon_with_holes P;
    bool inited = false;

    for(Polyhedron::Facet_const_iterator it = A.mesh_data.facets_begin(); it != A.mesh_data.facets_end(); it++)
    {
        Polygon Pit; Pit.clear();
        Polygon_with_holes Pith; Pith.clear();
        Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
        Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end;
        double n = 0;
        bool h1 = false, h2 = false, h3 = false;
        double x1, x2, x3, y1, y2, y3;
        do
        {
            double x, y, z, xP, yP;
            x = ::CGAL::to_double(hc_a->vertex()->point().x());
            y = ::CGAL::to_double(hc_a->vertex()->point().y());
            z = ::CGAL::to_double(hc_a->vertex()->point().z());
            xP = afftran(0,0)*x + afftran(1,0)*y + afftran(2,0)*z;
            yP = afftran(0,1)*x + afftran(1,1)*y + afftran(2,1)*z;
            Pit.push_back(Point_2(xP, yP));
            if(!h1)
            {
                x1 = xP;
                y1 = yP;
                h1 = true;
            }
            else if(!h2)
            {
                x2 = xP;
                y2 = yP;
                h2 = true;
            }
            else if(!h3)
            {
                x3 = xP;
                y3 = yP;
                h3 = true;
            }
            hc_a++;
        }while(hc_a != hc_end);
        n = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
        if(0.0 < n)
        {
            Pit.reverse_orientation();
        }
        inflatePolygon(Pit, 0.1);

        if(!inited)
        {
            inited = true;
            //std::cerr << "  <" << std::endl;
            CGAL::join(Pit, Pit, P);
            //std::cerr << "  >" << std::endl;
        }
        else
        {
            /*A bit of black magic, seems like Bool Ops are a bit more reliable on Polygs with
            holes. Therefore, init a Polygon with holes from the projection of a triangle and
            join it with the previously accumulated projection.*/
            //if(dbgK < 23)
            //{
            //std::cerr << "  <[" << std::endl;
            CGAL::join(Pit, Pit, Pith);
            //std::cerr << "  ||" << std::endl;
            CGAL::join(Pith, P, P);
            //std::cerr << "  ]>" << std::endl;
            //}
            dbgK++;
        }
    }

    Polygon oB = P.outer_boundary();
    CDT cdt;
    cdt.insert_constraint(oB.vertices_begin(), oB.vertices_end(), true);
    for(Polygon_with_holes::Hole_const_iterator it = P.holes_begin();
        it != P.holes_end(); it++)
        cdt.insert_constraint(it->vertices_begin(), it->vertices_end(), true);
    mark_domains(cdt);

    Build_from_triangulation<Polyhedron::HalfedgeDS> triangulationBuilder(cdt, afftran, fillHoles);
    FaceInfo2 fi;
    fi.in_domain();
    mesh_data.clear();
    mesh_data.delegate(triangulationBuilder);
    //Warning: this is a "thin" (non-solid) mesh, on conversion from Nef to Poly the result is empty.
    //Mesh should be stabilized by the user with a call to solidify.
    nef_polyhedron = Nef_polyhedron(mesh_data);
    props_updated = false;
    triangulated = false;

    return true;
}

bool MeshEntry::setFromSolidification(MeshEntry const& A, double thickness)
{
    if(&A != this)
    {
        mesh_data.clear();
        Build_solidification<Polyhedron::HDS> solidifier(A.mesh_data, thickness);
        mesh_data.delegate(solidifier);
    }
    else
    {
        Polyhedron copy(mesh_data);
        mesh_data.clear();
        Build_solidification<Polyhedron::HDS> solidifier(copy, thickness);
        mesh_data.delegate(solidifier);
    }
    props_updated = false;
    triangulated = false;
    nef_polyhedron = Nef_polyhedron(mesh_data);
    return true;
}

bool MeshEntry::getConvexComponents(std::vector<MeshEntry*> &components) const
{
    Nef_polyhedron cN = nef_polyhedron.regularization();
    CGAL::convex_decomposition_3(cN);
    cN = cN.interior();
    int k = 0;
    int bk = components.size();
    std::cerr << "CGAL did decomp into " << cN.number_of_volumes() - 1 << " components." <<std::endl;
    for(Nef_polyhedron::Volume_const_iterator it = ++cN.volumes_begin();
        it != cN.volumes_end(); it++)
        if(it->mark() && it->is_valid())
        {
            //CGAL::SNC_in_place_list_volume<CGAL::SNC_indexed_items::Volume<CGAL::SNC_structure<CGAL::Epeck, CGAL::SNC_indexed_items, bool> > > v;
            //v.is_valid();
            Polyhedron P, P2;
            P.clear(); P2.clear();
            cN.convert_inner_shell_to_polyhedron(it->shells_begin(), P);
            Nef_polyhedron N(P);
            N = N.interior();
            N.convert_to_polyhedron(P2);
            //trimesh::TriMesh M;
            //MeshEntry::write_to_trimesh(P, &M);
            //MeshEntry::remove_duplicates(&M, 0.00001);
            //Build_mesh<Polyhedron::HalfedgeDS> build_mesh(&M);
            //if(M.faces.size())
            //    P2.delegate(build_mesh);
            if((!P2.is_empty()) && P2.is_valid() && P2.is_closed() && (3 < P2.size_of_vertices()))
            {
                std::cerr << k << std::endl;
                components.push_back(new MeshEntry());
                components[bk+k]->mesh_data = P2;
                components[bk+k]->nef_polyhedron = Nef_polyhedron(P2);
                k++;
            }
        }
    return true;
}

bool MeshEntry::write_to_trimesh(Polyhedron const& P, trimesh::TriMesh *M)
{
    M->vertices.clear();
    M->faces.clear();
    M->vertices.reserve(P.size_of_vertices());
    M->faces.reserve(P.size_of_facets());
    for(Polyhedron::Vertex_const_iterator it = P.vertices_begin();
        it != P.vertices_end(); it++)
    {
        M->vertices.push_back(trimesh::point(::CGAL::to_double(it->point().x()),
                                             ::CGAL::to_double(it->point().y()),
                                             ::CGAL::to_double(it->point().z())));

    }

    typedef CGAL::Inverse_index<Polyhedron::Vertex_const_iterator> Index;
    Index index( P.vertices_begin(), P.vertices_end());

    for(Polyhedron::Facet_const_iterator it = P.facets_begin();
        it != P.facets_end(); it++)
    {
        Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
        Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end; hc_a++;
        Polyhedron::Halfedge_around_facet_const_circulator hc_b = hc_a; hc_b++;
        std::size_t n = circulator_size(hc_end);
        //n >= 3 should be guaranteed by the mesh ops from CGAL returning sane meshes, and we shouldn't load
        //insane meshes anyway.
        //note: actually, n should be equal to 3, guaranteed by a call to triangulate_polyhedron in
        //triangulate_mesh()
        do
        {
            M->faces.push_back(trimesh::TriMesh::Face(index[Polyhedron::Vertex_const_iterator(hc_end->vertex())],
                                                      index[Polyhedron::Vertex_const_iterator(hc_a->vertex())],
                                                      index[Polyhedron::Vertex_const_iterator(hc_b->vertex())]));
            hc_a++;
            hc_b++;
        }while(hc_b != hc_end);
    }
    return true;
}


bool MeshEntry::writeToFile(std::string const& filename)
{
    triangulate_mesh();
    trimesh::TriMesh M;
    MeshEntry::write_to_trimesh(mesh_data, &M);
    bool wrote_file = M.write(filename.c_str());
    return wrote_file;
}

bool MeshEntry::writeToMsg(shape_msgs::Mesh &message)
{
    triangulate_mesh();
    message.vertices.clear(); message.vertices.reserve(mesh_data.size_of_vertices());
    message.triangles.clear(); message.triangles.reserve(mesh_data.size_of_facets());
    for(Polyhedron::Vertex_const_iterator it = mesh_data.vertices_begin();
        it != mesh_data.vertices_end(); it++)
    {
        geometry_msgs::Point p;
        p.x = ::CGAL::to_double(it->point().x());
        p.y = ::CGAL::to_double(it->point().y());
        p.z = ::CGAL::to_double(it->point().z());
//        p.x = (it->point().x().approx().sup() + it->point().x().approx().inf())/2.0;
//        p.y = (it->point().y().approx().sup() + it->point().x().approx().inf())/2.0;
//        p.z = (it->point().z().approx().sup() + it->point().x().approx().inf())/2.0;
        message.vertices.push_back(p);
    }
    typedef CGAL::Inverse_index<Polyhedron::Vertex_const_iterator> Index;
    Index index( mesh_data.vertices_begin(), mesh_data.vertices_end());

    for(Polyhedron::Facet_const_iterator it = mesh_data.facets_begin();
        it != mesh_data.facets_end(); it++)
    {
        shape_msgs::MeshTriangle t;
        Polyhedron::Halfedge_around_facet_const_circulator hc_end = it->facet_begin();
        Polyhedron::Halfedge_around_facet_const_circulator hc_a = hc_end; hc_a++;
        Polyhedron::Halfedge_around_facet_const_circulator hc_b = hc_a; hc_b++;
        std::size_t n = circulator_size(hc_end);
        //n >= 3 should be guaranteed by the mesh ops from CGAL returning sane meshes, and we shouldn't load
        //insane meshes anyway.
        //note: actually, n should be equal to 3, guaranteed by a call to triangulate_polyhedron in
        //triangulate_mesh()
        do
        {
            t.vertex_indices[0] = index[Polyhedron::Vertex_const_iterator(hc_end->vertex())];
            t.vertex_indices[1] = index[Polyhedron::Vertex_const_iterator(hc_a->vertex())]; hc_a++;
            t.vertex_indices[2] = index[Polyhedron::Vertex_const_iterator(hc_b->vertex())]; hc_b++;
            message.triangles.push_back(t);
        }while(hc_b != hc_end);
    }
    return true;
}

bool MeshEntry::isClosed(void)
{
    update_properties();
    return is_closed;
}
bool MeshEntry::isManifold(void)
{
    update_properties();
    return is_manifold;
}
bool MeshEntry::isOrientable(void)
{
    update_properties();
    return is_orientable;
}
bool MeshEntry::isCSGSafe(void)
{
    update_properties();
    return is_csg_safe;
}
int MeshEntry::getNrVertices(void)
{
    update_properties();
    return nr_vertices;
}
int MeshEntry::getNrEdges(void)
{
    update_properties();
    return nr_edges;
}
int MeshEntry::getNrFaces(void)
{
    update_properties();
    return nr_faces;
}
int MeshEntry::getNrConnectedComponents(void)
{
    update_properties();
    return nr_connected_components;
}
int MeshEntry::getEulerCharacteristic(void)
{
    update_properties();
    return Euler_characteristic;
}

double MeshEntry::getVolume(void)
{
    update_properties();
    return volume;
}

bool MeshEntry::getVertices(std::vector<Point> &vertices) const
{
    vertices.resize(mesh_data.size_of_vertices());
    int k = 0;
    for(Polyhedron::Vertex_const_iterator it = mesh_data.vertices_begin();
        it != mesh_data.vertices_end(); it++)
        vertices[k++] = it->point();
    return true;
}

Polyhedron const& MeshEntry::getMesh(void) const
{
    return mesh_data;
}

void MeshEntry::remove_duplicates(trimesh::TriMesh *M, double duplicate_dist)
{
    std::vector<size_t> replacements; replacements.clear();
    MeshEntry::tree_type kdtree;
    kdtree.clear();
    double threshold = duplicate_dist;
    int maxK = M->vertices.size();
    replacements.reserve(maxK);
    for(int k = 0; k < maxK; k++)
    {
        MeshEntry::kdtree_node aux;
        aux.index = k;
        replacements.push_back(k);
        aux.xyz[0] = M->vertices[k][0];
        aux.xyz[1] = M->vertices[k][1];
        aux.xyz[2] = M->vertices[k][2];
        kdtree.insert(aux);
    }
    MeshEntry::tree_type::iterator it;
    it = kdtree.begin();
    for(int k = 0; k < maxK; k++, it++)
    {
        if(it->index == replacements[it->index])
        {
            std::vector<MeshEntry::kdtree_node> output; output.clear();
            kdtree.find_within_range(*it, threshold, std::back_inserter(output));
            int maxJ = output.size();
            for(int j = 0; j < maxJ; j++)
            {
                replacements[output[j].index] = it->index;
            }
        }
    }
    maxK = M->faces.size();
    for(int k = 0; k < maxK; )
    {
        M->faces[k].v[0] = replacements[M->faces[k].v[0]];
        M->faces[k].v[1] = replacements[M->faces[k].v[1]];
        M->faces[k].v[2] = replacements[M->faces[k].v[2]];
        if((M->faces[k].v[0] != M->faces[k].v[1]) &&
           (M->faces[k].v[1] != M->faces[k].v[2]) &&
           (M->faces[k].v[2] != M->faces[k].v[0]))
        {
            k++;
        }
        else
        {
            maxK--;
            M->faces[k] = M->faces[maxK];
        }
    }
    M->faces.resize(maxK);
    trimesh::remove_unused_vertices(M);
    replacements.clear();
    kdtree.clear();
}
double getSignedTriangleVolume(Polyhedron::Facet const& f)
{
    Polyhedron::Halfedge_around_facet_const_circulator hc_a = f.facet_begin();
    Polyhedron::Halfedge_around_facet_const_circulator hc_b = hc_a; hc_b++;
    Polyhedron::Halfedge_around_facet_const_circulator hc_c = hc_b; hc_c++;

    double xA, yA, zA;
    double xB, yB, zB;
    double xC, yC, zC;

    xA = ::CGAL::to_double(hc_a->vertex()->point().x());
    yA = ::CGAL::to_double(hc_a->vertex()->point().y());
    zA = ::CGAL::to_double(hc_a->vertex()->point().z());

    xB = ::CGAL::to_double(hc_b->vertex()->point().x());
    yB = ::CGAL::to_double(hc_b->vertex()->point().y());
    zB = ::CGAL::to_double(hc_b->vertex()->point().z());

    xC = ::CGAL::to_double(hc_c->vertex()->point().x());
    yC = ::CGAL::to_double(hc_c->vertex()->point().y());
    zC = ::CGAL::to_double(hc_c->vertex()->point().z());

    double v321 = xC*yB*zA;
    double v231 = xB*yC*zA;
    double v312 = xC*yA*zB;
    double v132 = xA*yC*zB;
    double v213 = xB*yA*zC;
    double v123 = xA*yB*zC;

    return (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

void MeshEntry::update_properties(void)
{
    if(!props_updated)
    {
        props_updated = true;
        is_closed = mesh_data.is_closed();
        is_manifold = nef_polyhedron.is_simple();
        is_orientable = true; //TODO: !!!! Need to write a proper test here
        is_csg_safe = is_closed && is_manifold && is_orientable;
        nr_vertices = mesh_data.size_of_vertices();
        nr_edges = (mesh_data.size_of_halfedges() >> 1);
        nr_faces = mesh_data.size_of_facets();
        Euler_characteristic = nr_vertices - nr_edges + nr_faces;
        nr_connected_components = nef_polyhedron.number_of_volumes() - 1;
        volume = 0.0;
        if(is_closed)
        {
            triangulate_mesh();
            for(Polyhedron::Facet_const_iterator it = mesh_data.facets_begin();
                it != mesh_data.facets_end(); it++)
            {
                volume += getSignedTriangleVolume(*it);
            }
            if(volume < 0.0)
                volume = -volume;
        }
    }
}

bool MeshEntry::triangulate_mesh(void)
{
    if(!triangulated)
    {
        triangulated = true;
        //TODO: test the quality of default CGAL triangulation here.
        //Some talk about replacing it with Delaunay triangulation for better results, from 2014:
        //https://github.com/vitalif/openscad/commit/f1cf6cb1137adaa8351949346b6d9ab3d7be88c7
        //https://github.com/openscad/openscad/issues/414
        CGAL::Triangulate_modifier_exact<Polyhedron> modifier;
        if(!mesh_data.is_pure_triangle())
        {
            //TODO: Replace the CGAL triangulation code, it's kinda unstable.
            mesh_data.delegate(modifier);
        }
        else
            ; // Already triangulated, so thankfully nothing to do.
    }
    return true;
}

bool MeshEntry::update_mesh(void)
{
    if(!nef_polyhedron.is_simple())
    {
        nef_polyhedron = nef_polyhedron.regularization();
    }
    if(nef_polyhedron.is_simple())
    {
        nef_polyhedron.convert_to_Polyhedron(mesh_data);
        props_updated = false;
        triangulated = false;
        return true;
    }
    else
    {
        nef_polyhedron = Nef_polyhedron(mesh_data);
        return false;
    }
}

bool MeshEntry::applyTransform(Eigen::Affine3d const& M)
{
    Nef_polyhedron::Aff_transformation_3 afftran(M(0,0), M(0,1), M(0,2), M(0,3),
                                                 M(1,0), M(1,1), M(1,2), M(1,3),
                                                 M(2,0), M(2,1), M(2,2), M(2,3));
    nef_polyhedron.transform(afftran);

    return update_mesh();
}

void MeshEntry::get_bounding_box_limit(double &limit, double coordinate, bool &limitSet, bool lt) const
{
    if((lt && (limit < coordinate)) ||
            ((!lt) && (limit > coordinate)) ||
            (!limitSet))
    {
        limitSet = true;
        limit = coordinate;
    }
}

Polyhedron MeshEntry::getBoundingBox(double scale) const
{
    double maxX = 0, maxY = 0, maxZ = 0, minX = 0, minY = 0, minZ = 0;
    getBoundingBox(maxX, minX, maxY, minY, maxZ, minZ, scale);
    Polyhedron P;
    Build_box<Polyhedron::HalfedgeDS> bbox(maxX, minX, maxY, minY, maxZ, minZ);
    P.delegate(bbox);
    return P;
}

void MeshEntry::getBoundingBox(double &maxX, double &minX, double &maxY, double &minY, double &maxZ, double &minZ, double scale) const
{
    maxX = 0; maxY = 0; maxZ = 0; minX = 0; minY = 0; minZ = 0;
    bool setMaxX = false, setMaxY = false, setMaxZ = false;
    bool setMinX = false, setMinY = false, setMinZ = false;
    for(Polyhedron::Vertex_const_iterator it = mesh_data.vertices_begin();
        it != mesh_data.vertices_end(); it++)
    {
        get_bounding_box_limit(maxX, ::CGAL::to_double(it->point().x()), setMaxX, true);
        get_bounding_box_limit(minX, ::CGAL::to_double(it->point().x()), setMinX, false);
        get_bounding_box_limit(maxY, ::CGAL::to_double(it->point().y()), setMaxY, true);
        get_bounding_box_limit(minY, ::CGAL::to_double(it->point().y()), setMinY, false);
        get_bounding_box_limit(maxZ, ::CGAL::to_double(it->point().z()), setMaxZ, true);
        get_bounding_box_limit(minZ, ::CGAL::to_double(it->point().z()), setMinZ, false);
    }
    maxX *= scale;
    minX *= scale;
    maxY *= scale;
    minY *= scale;
    maxZ *= scale;
    minZ *= scale;
}

bool MeshEntry::getNearVertices(double x, double y, double z, double distance, std::vector<MeshEntry::XYZTriplet> &points) const
{
    double distanceSq = distance*distance;
    for(Polyhedron::Vertex_const_iterator it = mesh_data.vertices_begin();
        it!= mesh_data.vertices_end(); it++)
    {
        double px, py, pz;
        px = ::CGAL::to_double(it->point().x());
        py = ::CGAL::to_double(it->point().y());
        pz = ::CGAL::to_double(it->point().z());
        double dx, dy, dz;
        dx = px - x; if(dx < 0) dx = -dx;
        dy = py - y; if(dy < 0) dy = -dy;
        dz = pz - z; if(dz < 0) dz = -dz;
        if((dx < distance) && (dy < distance) && (dz < distance) &&
                ((dx*dx + dy*dy + dz*dz) < distanceSq))
        {
            MeshEntry::XYZTriplet aux;
            aux.x = px;
            aux.y = py;
            aux.z = pz;
            points.push_back(aux);
        }
    }
    return true;
}

bool MeshEntry::getIntersectionPoint(Point const& A, Point const& B, Point & R, double &bA, double &bB, double &bC, int & fh)
{
    bool retq = false;
    double distance = 1000000.0;

    double d00;
    double d01;
    double d11;
    double d20;
    double d21;
    double denom;

    fh = 0;
    int k = 0;
    for(Polyhedron::Facet_iterator fit = mesh_data.facets_begin();
        fit != mesh_data.facets_end(); fit++, k++)
    {
        Polyhedron::Halfedge_around_facet_circulator ha = fit->facet_begin();
        Polyhedron::Halfedge_around_facet_circulator hb = ha; hb++;
        Polyhedron::Halfedge_around_facet_circulator hc = hb; hc++;
        Kernel::Triangle_3 triangle(ha->vertex()->point(), hb->vertex()->point(), hc->vertex()->point());
        Kernel::Line_3 line(A, B);
        CGAL::cpp11::result_of<Kernel::Intersect_3(Kernel::Line_3,  Kernel::Triangle_3)>::type result = CGAL::intersection(line, triangle);
        if(result)
        {
            Point candidate;
            if(const Kernel::Segment_3* s = boost::get<Kernel::Segment_3>(&*result))
            {
                double distS = ::CGAL::to_double(CGAL::squared_distance(s->source(), A));
                double distT = ::CGAL::to_double(CGAL::squared_distance(s->target(), A));
                if(distS < distT)
                    candidate = s->source();
                else
                    candidate = s->target();
            }
            else
                candidate = *(boost::get<Kernel::Point_3>(&*result));

            Kernel::Vector_3 v0 = hb->vertex()->point() - ha->vertex()->point();
            Kernel::Vector_3 v1 = hc->vertex()->point() - ha->vertex()->point();
            Kernel::Vector_3 v2 = candidate - ha->vertex()->point();

            double a00 = ::CGAL::to_double(v0*v0);
            double a01 = ::CGAL::to_double(v0*v1);
            double a11 = ::CGAL::to_double(v1*v1);
            double a20 = ::CGAL::to_double(v2*v0);
            double a21 = ::CGAL::to_double(v2*v1);
            double anom = a00*a11 - a01*a01;

            // A little hack: skip triangles in which it would be unstable to compute bar. coords.
            if(anom < 0.000001)
                continue;

            if(!retq)
            {
                R = candidate;
                fh = k;
                distance = ::CGAL::to_double(CGAL::squared_distance(candidate, A));
                d00 = a00;
                d01 = a01;
                d11 = a11;
                d20 = a20;
                d21 = a21;
                denom = anom;
            }
            else
            {
                double auxD = ::CGAL::to_double(CGAL::squared_distance(candidate, A));
                if(auxD < distance)
                {
                    distance = auxD;
                    R = candidate;
                    fh = k;
                    d00 = a00;
                    d01 = a01;
                    d11 = a11;
                    d20 = a20;
                    d21 = a21;
                    denom = anom;
                }
            }
            retq = true;
        }
        //else: no intersection
    }
    if(retq)
    {
        denom = 1.0/denom;
        bB = (d11*d20 - d01*d21)*denom;
        bC = (d00*d21 - d01*d20)*denom;
        bA = 1.0 - bB - bC;
    }
    return retq;
}

bool MeshEntry::getGeodesicPath(geometry_msgs::Point const& A, geometry_msgs::Point const& B, geometry_msgs::Point const& C, geometry_msgs::Point const& D, std::vector<geometry_msgs::Point> &path)
{

    Surface_mesh_shortest_path shortest_paths(mesh_data,
                                              get(boost::vertex_external_index, mesh_data),
                                              get(CGAL::halfedge_external_index, mesh_data),
                                              get(CGAL::face_external_index, mesh_data),
                                              get(CGAL::vertex_point, mesh_data));

    Point R;
    double bA, bB, bC;
    int faceSrc, faceDest;
    std::vector<Traits::Point_3> points;
    points.clear();
    path.clear();

    if(!getIntersectionPoint(Point(A.x, A.y, A.z), Point(B.x, B.y, B.z), R, bA, bB, bC, faceDest))
        return false;
    std::cerr << "Found intersection point ("
              << ::CGAL::to_double(R.x()) << " "
              << ::CGAL::to_double(R.y()) << " "
              << ::CGAL::to_double(R.z()) << ") baryc.("
              << bA << " " << bB << " " << bC << ")" << std::endl;
    Traits::Barycentric_coordinate faceLocationDest = {{bC, bA, bB}};
    if(!getIntersectionPoint(Point(C.x, C.y, C.z), Point(D.x, D.y, D.z), R, bA, bB, bC, faceSrc))
        return false;
    std::cerr << "Found intersection point ("
              << ::CGAL::to_double(R.x()) << " "
              << ::CGAL::to_double(R.y()) << " "
              << ::CGAL::to_double(R.z()) << ") baryc.("
              << bA << " " << bB << " " << bC << ")" << std::endl;
    Traits::Barycentric_coordinate faceLocationSrc = {{bC, bA, bB}};

    Polyhedron::Facet_iterator pfit = mesh_data.facets_begin();
    std::advance(pfit, faceSrc);
    face_iterator faceItSrc(pfit);
    pfit = mesh_data.facets_begin();
    std::advance(pfit, faceDest);
    face_iterator faceItDest(pfit);

    shortest_paths.add_source_point(*faceItSrc, faceLocationSrc);
    shortest_paths.shortest_path_points_to_source_points(*faceItDest, faceLocationDest, std::back_inserter(points));

    int maxK = points.size();
    for(int k = 0; k < maxK; k++)
    {
        geometry_msgs::Point aux;
        aux.x = ::CGAL::to_double(points[k].x());
        aux.y = ::CGAL::to_double(points[k].y());
        aux.z = ::CGAL::to_double(points[k].z());
        path.push_back(aux);
    }

    return true;
}

double MeshEntry::manhattan_distance(double xA, double yA, double zA,
                                     double xB, double yB, double zB)
{
    double dx = xA - xB; if(dx < 0.0) dx = -dx;
    double dy = yA - yB; if(dy < 0.0) dy = -dy;
    double dz = zA - zB; if(dz < 0.0) dz = -dz;
    double retq = dx;
    if(retq < dy)
        retq = dy;
    if(retq < dz)
        retq = dz;
    return retq;
}

double MeshEntry::euclidean_distance(double xA, double yA, double zA,
                                     double xB, double yB, double zB)
{
    double dx = xA - xB;
    double dy = yA - yB;
    double dz = zA - zB;
    double retq = std::sqrt(dx*dx + dy*dy + dz*dz);
    return retq;
}

bool MeshEntry::getMeshSkeleton(bool approximate, double duplicate_distance, std::vector<double> &x, std::vector<double> &y,
                                std::vector<double> &z, std::vector<int> &e_L, std::vector<int> &e_R) const
{
    x.clear();
    y.clear();
    z.clear();
    e_L.clear();
    e_R.clear();
    if(approximate)
    {
        trimesh::TriMesh M;
        M.vertices.clear();
        M.faces.clear();
        MeshEntry::write_to_trimesh(mesh_data, &M);
        int maxK = M.faces.size();
        std::vector<std::pair<int, int> > vertex_replacements;
        MeshEntry::init_sets(M.vertices.size(), vertex_replacements);
        for(int k = 0; k < maxK; k++)
        {
            double xA, xB, xC, yA, yB, yC, zA, zB, zC;
            int a, b, c;
            a = M.faces[k].v[0]; xA = M.vertices[a][0]; yA = M.vertices[a][1]; zA = M.vertices[a][2];
            b = M.faces[k].v[1]; xB = M.vertices[b][0]; yB = M.vertices[b][1]; zB = M.vertices[b][2];
            c = M.faces[k].v[2]; xC = M.vertices[c][0]; yC = M.vertices[c][1]; zC = M.vertices[c][2];
            if(MeshEntry::manhattan_distance(xA, yA, zA, xB, yB, zB) < duplicate_distance)
                MeshEntry::merge_sets(a, b, vertex_replacements);
            if(MeshEntry::manhattan_distance(xC, yC, zC, xB, yB, zB) < duplicate_distance)
                MeshEntry::merge_sets(c, b, vertex_replacements);
            if(MeshEntry::manhattan_distance(xA, yA, zA, xC, yC, zC) < duplicate_distance)
                MeshEntry::merge_sets(a, c, vertex_replacements);
        }

        int remaining_vertices = 0;
        maxK = M.vertices.size();
        for(int k = 0; k < maxK; k++)
            if(k == MeshEntry::get_set_index(k, vertex_replacements))
                remaining_vertices++;
        x.reserve(remaining_vertices);
        y.reserve(remaining_vertices);
        z.reserve(remaining_vertices);
        e_L.reserve(remaining_vertices);
        e_R.reserve(remaining_vertices);
        std::vector<int> new_indices;
        new_indices.clear(); new_indices.resize(maxK);
        for(int k = 0, ci = 0; k < maxK; k++)
        {
            if(k == MeshEntry::get_set_index(k, vertex_replacements))
            {
                new_indices[k] = ci; ci++;
                x.push_back(M.vertices[k][0]);
                y.push_back(M.vertices[k][1]);
                z.push_back(M.vertices[k][2]);
            }
        }
        std::vector<std::vector<bool> > marked_edge;
        marked_edge.clear(); marked_edge.resize(remaining_vertices);
        for(int k = 0; k < remaining_vertices; k++)
        {
            marked_edge[k].clear(); marked_edge[k].resize(remaining_vertices, false);
        }
        maxK = M.faces.size();
        for(int k = 0; k < maxK; k++)
        {
            int a, b, c;
            a = new_indices[MeshEntry::get_set_index(M.faces[k].v[0], vertex_replacements)];
            b = new_indices[MeshEntry::get_set_index(M.faces[k].v[1], vertex_replacements)];
            c = new_indices[MeshEntry::get_set_index(M.faces[k].v[2], vertex_replacements)];
            if((a != b) && (!marked_edge[a][b]))
            {
                marked_edge[a][b] = marked_edge[b][a] = true;
                e_L.push_back(a);
                e_R.push_back(b);
            }
            if((a != c) && (!marked_edge[a][c]))
            {
                marked_edge[a][c] = marked_edge[c][a] = true;
                e_L.push_back(a);
                e_R.push_back(c);
            }
            if((c != b) && (!marked_edge[c][b]))
            {
                marked_edge[c][b] = marked_edge[b][c] = true;
                e_L.push_back(c);
                e_R.push_back(b);
            }
        }
    }
    else
    {
        typedef CGAL::Simple_cartesian<double> KernelSimple;
        typedef KernelSimple::Point_3 PointSimple;
        typedef CGAL::Polyhedron_3<KernelSimple> PolyhedronSimple;
        typedef CGAL::Mean_curvature_flow_skeletonization<PolyhedronSimple>::Skeleton Skeleton;
        typedef Skeleton::vertex_descriptor SkeletonVertex;
        typedef Skeleton::edge_descriptor SkeletonEdge;

        PolyhedronSimple meshSimp;

        meshSimp.clear();
        Convert_mesh_kernel<PolyhedronSimple::HalfedgeDS, Polyhedron> build_mesh(mesh_data);
        meshSimp.delegate(build_mesh);

        Skeleton skeleton;
        std::vector<std::vector<bool> > wroteEdge; wroteEdge.clear();
        std::vector<bool> aux; aux.clear();
        std::map<SkeletonVertex, int> index; index.clear();
        CGAL::extract_mean_curvature_flow_skeleton(meshSimp, skeleton);
        int maxK = boost::num_vertices(skeleton);
        aux.resize(maxK, false);
        wroteEdge.resize(maxK, aux);
        x.reserve(maxK);
        y.reserve(maxK);
        z.reserve(maxK);
        e_L.reserve(boost::num_edges(skeleton));
        e_R.reserve(boost::num_edges(skeleton));
        int k = 0;
        BOOST_FOREACH(SkeletonVertex v, vertices(skeleton))
        {
            index.insert(std::pair<SkeletonVertex, int>(v, k));
            k++;
            x.push_back(::CGAL::to_double(skeleton[v].point.x()));
            y.push_back(::CGAL::to_double(skeleton[v].point.y()));
            z.push_back(::CGAL::to_double(skeleton[v].point.z()));
        }
        BOOST_FOREACH(SkeletonEdge e, edges(skeleton))
        {
            int indexA = index.find(source(e, skeleton))->second;
            int indexB = index.find(target(e, skeleton))->second;
            if(!wroteEdge[indexA][indexB])
            {
                wroteEdge[indexA][indexB] = true;
                wroteEdge[indexB][indexA] = true;
                e_L.push_back(indexA);
                e_R.push_back(indexB);
            }
        }
    }

    return true;
}

int MeshEntry::init_sets(int maxK, std::vector<std::pair<int, int> > &elements)
{
    elements.clear();
    elements.reserve(maxK);
    for(int k = 0; k < maxK; k++)
        elements.push_back(std::pair<int, int>(k, 0));
    return maxK;
}

int MeshEntry::get_set_index(int index, std::vector<std::pair<int, int> > &elements)
{
    if(index != elements[index].first)
    {
        elements[index].first = get_set_index(elements[index].first, elements);
    }
    return elements[index].first;
}

int MeshEntry::merge_sets(int index_A, int index_B, std::vector<std::pair<int, int> > &elements)
{
    int p_A = get_set_index(index_A, elements);
    int p_B = get_set_index(index_B, elements);

    if(p_A != p_B)
    {
        if(elements[p_A].second < elements[p_B].second)
        {
            elements[p_A].first = p_B;
            p_A = p_B;
        }
        else if(elements[p_B].second < elements[p_A].second)
        {
            elements[p_B].first = p_A;
        }
        else
        {
            elements[p_B].first = p_A;
            elements[p_A].second = elements[p_A].second + 1;
        }
        elements[p_B].first = p_A;
        elements[index_B].first = p_A;
    }
    return p_A;
}

}
