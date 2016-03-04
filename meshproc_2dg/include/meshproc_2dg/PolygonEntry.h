#ifndef __POLYGON_ENTRY__

#define __POLYGON_ENTRY__

namespace meshproc_2dg
{
class PolygonEntry
{
public:
    PolygonEntry(){}
    PolygonEntry(Polygon_with_holes_2 const& data):polygon(data){}
    ~PolygonEntry(){}

    bool loadFromMsg(meshproc_msgs::PolygonWithHoles const& msg);
    bool writeToMsg(meshproc_msgs::PolygonWithHoles & msg) const;
    bool convexDecomposition(std::vector<Polygon_2> & results, bool ignoreHoles, bool triangulate) const;

    static bool loadFromMsg(const geometry_msgs::Polygon &msg, Polygon_2 & polygon);
    static bool loadFromMsg(const meshproc_msgs::PolygonWithHoles &msg, Polygon_with_holes_2 & polygon);
    static bool writeToMsg(geometry_msgs::Polygon &msg, Polygon_2 const& polygon);
    static bool writeToMsg(meshproc_msgs::PolygonWithHoles &msg, Polygon_with_holes_2 const& polygon);
    static bool writeToMsg(meshproc_msgs::Skeleton &msg, Skeleton const& skeleton);
    static bool convexDecomposition(std::vector<Polygon_2> & results, Polygon_with_holes_2 const& polygon, bool ignoreHoles, bool triangulate);
    static bool csgRequest(Polygon_with_holes_2 const& polygon_A, Polygon_with_holes_2 const& polygon_B, std::list<Polygon_with_holes_2> & result, int operation);
    static bool visibility(Polygon_with_holes_2 const& polygon, std::vector<Point_2> const& points, std::vector<Polygon_2> & results);
    static bool extrudePolygon(Polyhedron &mesh, std::vector<meshproc_2dg::Polygon_2> &triangles, double extrudeHeight, double extrudeDepth);
private:
    static Polygon_2 Polygon_2_Conversion(Traits::Polygon_2 const& polygon);

    Polygon_with_holes_2 polygon;
};

}

#endif
