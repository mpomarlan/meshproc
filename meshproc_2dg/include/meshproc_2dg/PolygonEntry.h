#ifndef __POLYGON_ENTRY__

#define __POLYGON_ENTRY__

namespace meshproc_2dg
{
class PolygonEntry
{
public:
    PolygonEntry(){}
    PolygonEntry(Polygon_2 const& data):polygon(data){}
    ~PolygonEntry(){}

    bool loadFromMsg(geometry_msgs::Polygon const& msg);
    bool writeToMsg(geometry_msgs::Polygon & msg) const;
    bool convexDecomposition(std::vector<Polygon_2> & results) const;

    static bool loadFromMsg(const geometry_msgs::Polygon &msg, Polygon_2 & polygon);
    static bool writeToMsg(geometry_msgs::Polygon &msg, Polygon_2 const& polygon);
    static bool convexDecomposition(std::vector<Polygon_2> & results, Polygon_2 const& polygon);
private:
    Polygon_2 polygon;
};

}

#endif
