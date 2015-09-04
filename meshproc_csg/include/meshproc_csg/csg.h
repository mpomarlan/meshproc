#ifndef __MESHPROC_CSG__

#define __MESHPROC_CSG__

namespace meshproc_csg{

#define __USE_EXACT_KERNEL__

#ifdef __USE_EXACT_KERNEL__
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
#endif
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;

class MeshEntry
{
public:

    typedef struct {double x, y, z;} XYZTriplet;

    MeshEntry();
    ~MeshEntry();

    void clear(void);
    bool loadFromFile(std::string const& filename, double duplicate_dist);
    bool loadFromMsg(shape_msgs::Mesh const& message, double duplicate_dist);

    bool setFromUnion(MeshEntry const& A, MeshEntry const& B);
    bool setFromIntersection(MeshEntry const& A, MeshEntry const& B);
    bool setFromDifference(MeshEntry const& A, MeshEntry const& B);
    bool setFromSymmetricDifference(MeshEntry const& A, MeshEntry const& B);
    bool setFromMinkowskiSum(MeshEntry const& A, MeshEntry const& B);
    bool setFromMinkowskiErosion(MeshEntry const& A, MeshEntry const& B);
    bool setFromSelectComponent(MeshEntry const& A, meshproc_csg::Point const& P);

    bool getConvexComponents(std::vector<MeshEntry*> &components) const;

    bool getMeshSkeleton(bool approximate, double duplicate_distance, std::vector<double> &x, std::vector<double> &y,
                         std::vector<double> &z, std::vector<int> &e_L, std::vector<int> &e_R) const;

    bool writeToFile(std::string const& filename);
    bool writeToMsg(shape_msgs::Mesh &message);

    bool isClosed(void);
    bool isManifold(void);
    bool isOrientable(void);
    bool isCSGSafe(void);
    int getNrVertices(void);
    int getNrEdges(void);
    int getNrFaces(void);
    int getNrConnectedComponents(void);
    int getEulerCharacteristic(void);

    Polyhedron getBoundingBox(double scale=1.0) const;
    void getBoundingBox(double &maxX, double &minX, double &maxY, double &minY, double &maxZ, double &minZ, double scale=1.0) const;

    /*The concept is: apply transformations lazily.
      applyTransform just changes the Eigen::Affine3d object transform.
      It's by a call to transformMesh that the transformation actually
      happens to the mesh and NEF polyhedron.*/
    bool applyTransform(Eigen::Affine3d const& M, MeshEntry *baseMesh, bool incremental);
    bool transformMesh(void);
    bool addDependent(MeshEntry * meshEntry);
    bool getDependents(std::vector<MeshEntry *> & meshes) const;
    bool transformDependents(void);

    bool getNearVertices(double x, double y, double z, double distance, std::vector<MeshEntry::XYZTriplet> &points) const;

private:
    bool loadFromTrimesh(trimesh::TriMesh *M, double duplicate_dist);
    static double manhattan_distance(double xA, double yA, double zA,
                                     double xB, double yB, double zB);
    static bool write_to_trimesh(Polyhedron const& P, trimesh::TriMesh *M);
    static void remove_duplicates(trimesh::TriMesh *M, double duplicate_dist);
    void update_properties(void);
    bool update_mesh(void);
    bool triangulate_mesh(void);
    bool process_transforms(void);

    void get_bounding_box_limit(double &limit, double coordinate, bool &limitSet, bool lt) const;

    static int init_sets(int maxK, std::vector<std::pair<int, int> > &elements);
    static int get_set_index(int index, std::vector<std::pair<int, int> > &elements);
    static int merge_sets(int index_A, int index_B, std::vector<std::pair<int, int> > &elements);

    std::vector<MeshEntry *> dependents;
    MeshEntry const* base_mesh;
    Eigen::Affine3d transform;
    bool transform_empty;

    struct kdtree_node
    {
        typedef double value_type;
        value_type xyz[3];
        size_t index;
        value_type operator[](size_t n) const
        {
            return xyz[n];
        }
        value_type distance(const kdtree_node &node)
        {
            value_type dx = xyz[0] - node.xyz[0];
            value_type dy = xyz[1] - node.xyz[1];
            value_type dz = xyz[2] - node.xyz[2];
            return std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    };
    typedef KDTree::KDTree<3, kdtree_node> tree_type;
    std::string name;
    meshproc_csg::Polyhedron mesh_data;
    meshproc_csg::Nef_polyhedron nef_polyhedron;
    bool props_updated;
    bool triangulated;
    bool is_closed;
    bool is_manifold;
    bool is_orientable;
    bool is_csg_safe;
    int nr_vertices;
    int nr_edges;
    int nr_faces;
    int nr_connected_components;
    int Euler_characteristic;
};

typedef std::map<std::string, MeshEntry*> MeshMap;

}

#endif
