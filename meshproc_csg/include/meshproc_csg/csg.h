#ifndef __MESHPROC_CSG__

#define __MESHPROC_CSG__

namespace meshproc_csg{

class MeshEntry
{
public:

    typedef struct {double x, y, z;} XYZTriplet;

    MeshEntry();
    MeshEntry(MeshEntry const& orig):
        name(""),
        mesh_data(orig.mesh_data),
        nef_polyhedron(orig.nef_polyhedron),
        props_updated(orig.props_updated),
        triangulated(orig.triangulated),
        is_closed(orig.is_closed),
        is_manifold(orig.is_manifold),
        is_orientable(orig.is_orientable),
        is_csg_safe(orig.is_csg_safe),
        nr_vertices(orig.nr_vertices),
        nr_edges(orig.nr_edges),
        nr_faces(orig.nr_faces),
        nr_connected_components(orig.nr_connected_components),
        Euler_characteristic(orig.Euler_characteristic),
        volume(orig.volume)
        {}
    MeshEntry(Polyhedron const& polyhedron);
    ~MeshEntry();

    void clear(void);
    bool loadFromFile(std::string const& filename, double duplicate_dist);
    bool loadFromMsg(shape_msgs::Mesh const& message, double duplicate_dist);

    bool setFromMeshEntry(MeshEntry const& A);

    bool setFromUnion(MeshEntry const& A, MeshEntry const& B);
    bool setFromIntersection(MeshEntry const& A, MeshEntry const& B);
    bool setFromDifference(MeshEntry const& A, MeshEntry const& B);
    bool setFromSymmetricDifference(MeshEntry const& A, MeshEntry const& B);
    bool setFromMinkowskiSum(MeshEntry const& A, MeshEntry const& B);
    bool setFromMinkowskiErosion(MeshEntry const& A, MeshEntry const& B);
    bool setFromSelectComponent(MeshEntry const& A, meshproc_csg::Point const& P);
    bool setFromConvexHull(MeshEntry & A);
    bool setFromProjection(MeshEntry const& A, double nx, double ny, double nz, bool fillHoles);
    bool setFromSolidification(MeshEntry const& A, double thickness);
    bool setFromMesh2Prism(MeshEntry const& A, double height, double depth, bool filter = false, double zcomp = 0);
    bool filterByNormal(std::string const& filename, double nx, double ny, double nz, double toleratedAngle);

    bool getConvexComponents(std::vector<MeshEntry*> &components) const;

    bool getBoundaryPolygon(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                            std::vector<int> &edge_L, std::vector<int> &edge_R);

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
    double getVolume(void);

    bool getVertices(std::vector<Point> &vertices) const;
    Polyhedron const& getMesh(void) const;

    Polyhedron getBoundingBox(double scale=1.0) const;
    void getBoundingBox(double &maxX, double &minX, double &maxY, double &minY, double &maxZ, double &minZ, double scale=1.0) const;

    bool applyTransform(Eigen::Affine3d const& M);

    bool getNearVertices(double x, double y, double z, double distance, std::vector<MeshEntry::XYZTriplet> &points) const;

protected:
    bool loadFromTrimesh(trimesh::TriMesh *M, double duplicate_dist);
    static double manhattan_distance(double xA, double yA, double zA,
                                     double xB, double yB, double zB);
    static double euclidean_distance(double xA, double yA, double zA,
                                     double xB, double yB, double zB);
    static bool write_to_trimesh(Polyhedron const& P, trimesh::TriMesh *M);
    static void remove_duplicates(trimesh::TriMesh *M, double duplicate_dist);
    void update_properties(void);
    bool update_mesh(void);
    bool triangulate_mesh(void);

    void get_bounding_box_limit(double &limit, double coordinate, bool &limitSet, bool lt) const;

    static int init_sets(int maxK, std::vector<std::pair<int, int> > &elements);
    static int get_set_index(int index, std::vector<std::pair<int, int> > &elements);
    static int merge_sets(int index_A, int index_B, std::vector<std::pair<int, int> > &elements);

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
    double volume;
};

typedef std::map<std::string, MeshEntry*> MeshMap;

}

#endif
