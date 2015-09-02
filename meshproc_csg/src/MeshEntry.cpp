#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/minkowski_sum_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/triangulate_polyhedron.h>
#include <CGAL/bounding_box.h>

#include <trimesh/TriMesh.h>
#include <trimesh/TriMesh_algo.h>

#include <meshproc_csg/kdtree++/kdtree.hpp>

#include <shape_msgs/Mesh.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <meshproc_csg/csg.h>

#include <meshproc_csg/build_mesh.h>
#include <meshproc_csg/triangulate.h>

namespace meshproc_csg
{

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
    transform_empty = true;
    base_mesh = this;
    transform = Eigen::Affine3d();
}

MeshEntry::~MeshEntry()
{
}

void MeshEntry::clear(void)
{
    props_updated = false;
    triangulated = false;
    mesh_data.clear();
    transform_empty = true;
    base_mesh = this;
    nef_polyhedron.clear();
}
bool MeshEntry::loadFromFile(std::string const& filename, double duplicate_sq_dist)
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
    loadFromTrimesh(M, duplicate_sq_dist);
    return true;
}
bool MeshEntry::loadFromMsg(shape_msgs::Mesh const& message, double duplicate_sq_dist)
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
    loadFromTrimesh(&M, duplicate_sq_dist);
}
bool MeshEntry::loadFromTrimesh(trimesh::TriMesh *M, double duplicate_sq_dist)
{
    remove_duplicates(M, duplicate_sq_dist);
    Build_mesh<Polyhedron::HalfedgeDS> build_mesh(M);
    mesh_data.delegate(build_mesh);
    nef_polyhedron = Nef_polyhedron(mesh_data);
}

bool MeshEntry::transformDependents(void)
{
    int maxK = dependents.size();
    for(int k = 0; k < maxK; k++)
        dependents[k]->transformMesh();
    dependents.clear();
    return true;
}

bool MeshEntry::process_transforms(void)
{
    transformDependents();
    base_mesh = this;
    transform_empty = true;
    return true;
}

bool MeshEntry::setFromUnion(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
    nef_polyhedron = A.nef_polyhedron + B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromIntersection(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
    nef_polyhedron = A.nef_polyhedron * B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromDifference(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
    nef_polyhedron = A.nef_polyhedron - B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromSymmetricDifference(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
    nef_polyhedron = A.nef_polyhedron ^ B.nef_polyhedron;
    return update_mesh();
}
bool MeshEntry::setFromMinkowskiSum(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
    Nef_polyhedron cA = A.nef_polyhedron;
    Nef_polyhedron cB = B.nef_polyhedron;
    std::cout << "Copied nefs ..." << std::endl;
    nef_polyhedron = CGAL::minkowski_sum_3(cA, cB);
    return update_mesh();
}
bool MeshEntry::setFromMinkowskiErosion(MeshEntry const& A, MeshEntry const& B)
{
    process_transforms();
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
    process_transforms();
    return true;
}

bool MeshEntry::writeToFile(std::string const& filename)
{
    transformMesh();
    triangulate_mesh();
    trimesh::TriMesh M;
    M.vertices.clear();
    M.faces.clear();
    M.vertices.reserve(mesh_data.size_of_vertices());
    M.faces.reserve(mesh_data.size_of_facets());
    for(Polyhedron::Vertex_const_iterator it = mesh_data.vertices_begin();
        it != mesh_data.vertices_end(); it++)
    {
        M.vertices.push_back(trimesh::point(::CGAL::to_double(it->point().x()),
                                            ::CGAL::to_double(it->point().y()),
                                            ::CGAL::to_double(it->point().z())));

    }

    typedef CGAL::Inverse_index<Polyhedron::Vertex_const_iterator> Index;
    Index index( mesh_data.vertices_begin(), mesh_data.vertices_end());

    for(Polyhedron::Facet_const_iterator it = mesh_data.facets_begin();
        it != mesh_data.facets_end(); it++)
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
            M.faces.push_back(trimesh::TriMesh::Face(index[Polyhedron::Vertex_const_iterator(hc_end->vertex())],
                                                     index[Polyhedron::Vertex_const_iterator(hc_a->vertex())],
                                                     index[Polyhedron::Vertex_const_iterator(hc_b->vertex())]));
            hc_a++;
            hc_b++;
        }while(hc_b != hc_end);
    }

    bool wrote_file = M.write(filename.c_str());
    return wrote_file;
}

bool MeshEntry::writeToMsg(shape_msgs::Mesh &message)
{
    transformMesh();
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

void MeshEntry::remove_duplicates(trimesh::TriMesh *M, double duplicate_sq_dist)
{
    std::vector<size_t> replacements; replacements.clear();
    MeshEntry::tree_type kdtree;
    kdtree.clear();
    double threshold = std::sqrt(duplicate_sq_dist);
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
    for(int k = 0; k < maxK; k++)
    {
        M->faces[k].v[0] = replacements[M->faces[k].v[0]];
        M->faces[k].v[1] = replacements[M->faces[k].v[1]];
        M->faces[k].v[2] = replacements[M->faces[k].v[2]];
    }
    trimesh::remove_unused_vertices(M);
    replacements.clear();
    kdtree.clear();
}
void MeshEntry::update_properties(void)
{
    transformMesh();
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
        nr_connected_components = 1; //TODO: !!!! Need to write a proper counter here
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
        mesh_data.delegate(modifier);
    }
}

bool MeshEntry::update_mesh(void)
{
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

bool MeshEntry::applyTransform(Eigen::Affine3d const& M, MeshEntry *baseMesh, bool incremental)
{
    if(baseMesh != this)
        baseMesh->transformMesh();
    if(transform_empty || (base_mesh != baseMesh) || (!incremental))
        transform = M;
    else
        transform = M*transform;
    base_mesh = baseMesh;
    transform_empty = false;
    return true;
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


bool MeshEntry::transformMesh(void)
{
    transformDependents();
    if(!transform_empty)
    {
        if(base_mesh != this)
            nef_polyhedron = base_mesh->nef_polyhedron;
        Nef_polyhedron::Aff_transformation_3 afftran(transform(0,0), transform(0,1), transform(0,2), transform(0,3),
                                                     transform(1,0), transform(1,1), transform(1,2), transform(1,3),
                                                     transform(2,0), transform(2,1), transform(2,2), transform(2,3));
        nef_polyhedron.transform(afftran);
        transform_empty = update_mesh();
        if(transform_empty)
            base_mesh = this;
    }
    return transform_empty;
}

bool MeshEntry::addDependent(MeshEntry * meshEntry)
{
    dependents.push_back(meshEntry);
    return true;
}

bool MeshEntry::getDependents(std::vector<MeshEntry *> & meshes) const
{
    int maxK = dependents.size();
    for(int k = 0; k < maxK; k++)
        meshes.push_back(dependents[k]);
    return true;
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

}
