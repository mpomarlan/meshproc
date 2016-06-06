#ifndef __POINT_CLOUD_ENTRY__

#define __POINT_CLOUD_ENTRY__

namespace meshproc_p2m
{
class PointCloudEntry
{
public:
    PointCloudEntry():mesh(6, 200){}
    ~PointCloudEntry(){}

    bool loadFromFile(std::string const& filename);
    bool cloud2Mesh(double cellSize);
    bool writeMeshToFile(double gridsize, std::string const& filename) const;

private:
    static bool write_to_trimesh(double gridsize, std::vector<AFSFacet> const& R, Point_collection const& P, trimesh::TriMesh *M);
    Reconstruction mesh;
    Point_collection points;
    std::vector<AFSFacet> facets;
};

}

#endif
