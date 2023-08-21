#ifndef THREEDIMENSIONALSHAPE_H
#define THREEDIMENSIONALSHAPE_H

#include "Mesh.h"
#include "NonManifoldMesh/nonmanifoldmesh.h"
#include "SlabMesh.h"
class ThreeDimensionalShape
{
public:
    ThreeDimensionalShape();
    void ComputeInputNMM();

    // load the user simplified ma
    void LoadInputNMM(std::string fname);

    long PrepareSimplifySlabMesh();

    // initial the matrix of each face and vertex for slab mesh
    void InitialSlabMesh();

    double NearestPoint(Vector3d point, unsigned vid);

    // compute the Hausdorff distance
    void ComputeHausdorffDistance();

    void PruningSlabMesh();
    double distCenterToBoundary(Point_t center);
    double distCenterToBoundary(const Cell_handle_t & fci,Point center);
    int facetCount(const Finite_cells_iterator_t & fci);
private:
    FT nearestPointOfLine(Point p,Point a,Point b);
    int checkPointInTriangle(Point p,Point a,Point b,Point c);
    void createSharpSlabEdgeAndFace(int tag,std::pair<int, int> edge, std::map<int, int>& shapePointIndexMap, std::map<std::pair<int, int>, int>& shapeEdgeIndexMap, std::map<std::tuple<int, int, int>, int>& shapeFaceIndexMap);
    void createSharpSlabTetrahedron(int tag,std::tuple<int, int, int> face, std::map<int, int>& shapePointIndexMap, std::map<std::pair<int, int>, int>& shapeEdgeIndexMap, std::map<std::tuple<int, int, int>, int>& shapeFaceIndexMap);
    void createSlabFace(std::tuple<int,int,int> face,std::map<int,int>& shapePointIndexMap,std::map<std::pair<int,int>,int>& shapeEdgeIndexMap,std::map<std::tuple<int,int,int>,int>& shapeFaceIndexMap);
public:
    Mesh input;		// the mesh of the input shape 用来四面体化的数据
    Surface_mesh surface_mesh;
    unsigned num_vor_v, num_vor_e, num_vor_f;

    NonManifoldMesh input_nmm;//medial mesh data for compute and save

    SlabMesh slab_mesh;//medial mesh data for render and simplify

    bool slab_initial;
};

#endif // THREEDIMENSIONALSHAPE_H
