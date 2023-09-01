#ifndef _SLABMESH_H
#define _SLABMESH_H
#undef slots
#include <torch/torch.h>
#define slots Q_SLOTS
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "PrimMesh.h"
#include "geometry/Search/Triangle.h"
#include "geometry/Search/BVH.h"
#include <QVector3D>
class SlabPrim
{
public:
    Wm4::Matrix4d slab_A;//是论文矩阵的两倍，在计算cost时又乘以0.5
    Wm4::Vector4d slab_b;
    double slab_c;

    Wm4::Matrix4d add_A;//边界保护项或者中轴锥项
    Wm4::Vector4d add_b;//边界保护项或者中轴锥项
    double add_c;//边界保护项或者中轴锥项

    // hyperbolic weight
    double hyperbolic_weight;//边的稳定性系数

    SlabPrim() : slab_c(0.0), add_c(0.0), hyperbolic_weight(0.0){}
};

class SlabVertex : public PrimVertex, public SlabPrim
{
public:
    bool is_pole;
    bool is_non_manifold;
    bool is_disk;
    bool is_boundary;
};

class SlabEdge : public PrimEdge, public SlabPrim
{
};
class SlabEnvelop
{
public:
    std::vector<Vector3d> vertices;
    std::vector<std::array<unsigned int,3>> indexTris;
};

class SlabFace : public PrimFace, public SlabPrim
{
public:
    SlabEnvelop envelop;
};

typedef std::pair<bool, SlabVertex*> Bool_SlabVertexPointer;
typedef std::pair<bool, SlabEdge*> Bool_SlabEdgePointer;
typedef std::pair<bool, SlabFace*> Bool_SlabFacePointer;

class SlabMesh : public PrimMesh
{
public:
    Tree* tree;
    //MyCGAL::Primitives::BVHAccel<double>* tree;
    MyCGAL::Primitives::BVHAccel<double>* bvh=nullptr;
    MyCGAL::Primitives::BVHAccel<double>* facetsBvh=nullptr;
    Surface_mesh surface_mesh;
    double diag;
    std::map<unsigned,std::pair<Point,double>> surface2MatMap;
    //std::map<Point,Point> mat2SurfaceMap;
    std::map<MyCGAL::Primitives::Vector3d,std::vector<MyCGAL::Primitives::Vector3d>> mat2SurfaceMap;
    std::map<face_descriptor,Vector> fnormals;
    std::map<vertex_descriptor,Vector> vnormals;
public:
     std::map<unsigned,std::pair<Point,double>> optimSurface2MatMap;
public:
    SlabMesh():tree(nullptr),bvh(nullptr),facetsBvh(nullptr),compute_hausdorff(false),hasComputedHausdorff(false){

    }
    ~SlabMesh(){
        if(tree){
            delete tree;
        }
    }
    void constructTree(Surface_mesh &surface_mesh){
        this->surface_mesh=surface_mesh;
        int counter=0;
        for ( Surface_mesh::Vertex_iterator it = this->surface_mesh.vertices_begin (); it != this->surface_mesh.vertices_end (); ++it ) {
            it->id () = counter;
            counter++;
        }
        counter = 0;
        for (Surface_mesh::Face_iterator it = this->surface_mesh.facets_begin(); it != this->surface_mesh.facets_end(); ++it) {
            it->id() = counter;
            counter++;
        }

        CGAL::Polygon_mesh_processing::compute_normals(this->surface_mesh,
                                                       boost::make_assoc_property_map(vnormals),
                                                       boost::make_assoc_property_map(fnormals));
        tree=new Tree(CGAL::faces(this->surface_mesh).first,CGAL::faces(this->surface_mesh).second, this->surface_mesh);
        tree->accelerate_distance_queries();
        CGAL::Iso_cuboid_3<simple_kernel> bbox = CGAL::bounding_box(this->surface_mesh.points_begin(), this->surface_mesh.points_end());
        // calculate radius of bounding box

        diag = 0.5 * sqrt((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin())
                               + (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin())
                               + (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()));
        //auto center = (bbox.first + bbox.second) / 2.0;
//        std::vector<MyCGAL::Primitives::Object<double>*>objects;
//        for(face_descriptor fd: CGAL::faces(surface_mesh)){
//            Point p0=fd->halfedge()->vertex()->point();
//            Point p1=fd->halfedge()->next()->vertex()->point();
//            Point p2=fd->halfedge()->next()->next()->vertex()->point();
//            MyCGAL::Primitives::Vector3d v0(p0.x(),p0.y(),p0.z());
//            MyCGAL::Primitives::Vector3d v1(p1.x(),p1.y(),p1.z());
//            MyCGAL::Primitives::Vector3d v2(p2.x(),p2.y(),p2.z());
//            MyCGAL::Primitives::Triangle<double> *tri=new MyCGAL::Primitives::Triangle<double>(v0,v1,v2);
//            objects.push_back(tri);
//        }
//        tree=new MyCGAL::Primitives::BVHAccel<double>(objects);
    }
    void deleteWrongEdges();
    MyCGAL::Primitives::BVHAccel<double>* constructBVH();
    MyCGAL::Primitives::BVHAccel<double>* constructMeshFacetsBVH(Surface_mesh &surface_mesh)
    {
        int counter=0;
        std::vector<MyCGAL::Primitives::Object<double>*>objects;
        for(face_descriptor fd: CGAL::faces(surface_mesh)){
            Point p0=fd->halfedge()->vertex()->point();
            Point p1=fd->halfedge()->next()->vertex()->point();
            Point p2=fd->halfedge()->next()->next()->vertex()->point();
            MyCGAL::Primitives::Vector3d v0(p0.x(),p0.y(),p0.z());
            MyCGAL::Primitives::Vector3d v1(p1.x(),p1.y(),p1.z());
            MyCGAL::Primitives::Vector3d v2(p2.x(),p2.y(),p2.z());
            MyCGAL::Primitives::Triangle<double> *tri=new MyCGAL::Primitives::Triangle<double>(v0,v1,v2);
            tri->id=counter++;
            objects.push_back(tri);
        }
        MyCGAL::Primitives::BVHAccel<double> *tree=new MyCGAL::Primitives::BVHAccel<double>(objects);
        return tree;
    }
public:
    void project();
    void inverseProject();
    torch::Tensor laplacian_packed();
    torch::Tensor laplacian_weights();
public:
    std::vector<Bool_SlabVertexPointer> vertices;
    std::vector<Bool_SlabEdgePointer> edges;
    std::vector<Bool_SlabFacePointer> faces;
    std::vector<QVector3D> matColors;
    std::vector<QVector3D> matEdgeColors;
    std::vector<QVector3D> matFaceColors;

    std::vector<QVector3D> colors;
    std::vector<QVector3D> edgeColors;
    std::vector<QVector3D> faceColors;

    double hausdorff1;
    double hausdorff2;
public:
    void update();
    void updateSize();
    void compute();
public:
    // 1. preserve method one
    // 2. preserve method two
    // 3. preserve method three
    int preserve_boundary_method;

    // 0. no hyperbolic weight
    // 1. add hyperbolic distance to the related edges
    // 2. add hyperbolic area to the related face
    // 3. add ratio of hyperbolic and Euclid to the related edges
    int hyperbolic_weight_type;

    // 1.compute the boundary vertices only
    // 2.compute the boundary vertices and its related vertices
    int boundary_compute_scale;

    // the influence factor to control the ratio between hyperbolic and Euclid distance
    double k;

    // whether clear the error before initialization
    bool clear_error;

    bool preserve_saved_vertices;

    bool compute_hausdorff=false;

    bool hasComputedHausdorff=false;

    bool prevent_inversion;

    bool initial_boundary_preserve;

    double m_min[3];
    double m_max[3];

    unsigned simplified_inside_edges;
    unsigned simplified_boundary_edges;

    double bound_weight;//边界边权值取0.1
    Mesh_domain * domain=nullptr;
public:
    void AdjustStorage();

public:
    bool ValidVertex(unsigned vid);
    bool Edge(unsigned vid0, unsigned vid1, unsigned & eid);
    bool Face(const std::set<unsigned> & vset, unsigned & fid);
    void UpdateCentroid(unsigned fid);
    void ComputeFacesCentroid();
    void UpdateNormal(unsigned fid);
    void ComputeFacesNormal();
    void UpdateVertexNormal(unsigned vid);
    void ComputeVerticesNormal();
    void GetNeighborVertices(unsigned vid, std::set<unsigned> & neighborvertices);
    void GetLinkedEdges(unsigned eid, std::set<unsigned> & neighboredges);
    void GetAdjacentFaces(unsigned fid, std::set<unsigned> & neighborfaces);
    bool Contractible(unsigned vid_src, unsigned vid_tgt);
    bool Contractible(unsigned vid_src1, unsigned vid_src2, Vector3d v_tgt);
    bool MergeVertices(unsigned vid_src1, unsigned vid_src2, unsigned &vid_tgt);

    unsigned VertexIncidentEdgeCount(unsigned vid);
    unsigned VertexIncidentFaceCount(unsigned vid);
    unsigned EdgeIncidentFaceCount(unsigned eid);


public:
    void DeleteFace(unsigned fid);
    void DeleteEdge(unsigned eid);
    void DeleteVertex(unsigned vid);

    void InsertVertex(SlabVertex *vertex, unsigned &vid);
    void InsertEdge(unsigned vid0, unsigned vid1, unsigned & eid);
    void InsertFace(std::set<unsigned> vset);

    void ComputeEdgeCone(unsigned eid);
    void ComputeEdgesCone();
    void ComputeSlabEnvelop(unsigned fid);
    void ComputeVertexProperty(unsigned vid);
    void ComputeVerticesProperty();
    void ComputeFaceSimpleTriangles(unsigned fid);
    void ComputeFacesSimpleTriangles();

public:
    void initBoundaryCollapseQueue();
    void initCollapseQueue();
    void Simplify(int threshold);
    void SimplifyBoudary(int threshold);
    bool MinCostBoundaryEdgeCollapse(unsigned & eid);
    bool MinCostEdgeCollapse(unsigned & eid);
    void EvaluateEdgeCollapseCost(unsigned eid);
    void EvaluateEdgeHausdorffCost(unsigned eid);
    void ReEvaluateEdgeHausdorffCost(unsigned eid);
    void checkSkeletonAndTube();
    double computeHausdorff1();
    double computeHausdorff2();
public:
    void DistinguishVertexType();
    unsigned GetSavedPointNumber();
    unsigned GetConnectPointNumber();
    void InsertSavedPoint(unsigned vid);
    double NearestPoint(Vector3d point, unsigned vid);

public:
    void PreservBoundaryMethodOne();
    void PreservBoundaryMethodTwo();
    void PreservBoundaryMethodThree();
    void PreservBoundaryMethodFour();
    void PreservBoundaryMethodFive();

    void GetEnvelopeSet(const Vector4d & lamder, const set<unsigned> & neighbor_v, const set< std::set<unsigned> > & adj_faces, vector<Sphere> & sph_vec, vector<Cone> & con_vec, vector<SimpleTriangle> & st_vec);
    double EvaluateVertexDistanceErrorEnvelope(Vector4d & lamder, set<unsigned> & neighbor_vertices, set< set<unsigned> > & neighbor_faces, set<unsigned> & bplist);

    double GetHyperbolicLength(unsigned eid);
    double GetRatioHyperbolicEuclid(unsigned eid);

    void ExportSimplifyResult();
    void Export(std::string fname,bool backup=false);

public:
    void clear();
    void RecomputerVertexType();
    void computebb();

    void CleanIsolatedVertices();
    void InitialTopologyProperty(unsigned vid);
    void InitialTopologyProperty();
};
class PointToMatLoss
{
public:
    PointToMatLoss(){}
    torch::Tensor forward(SlabMesh *mesh,torch::Tensor& v0,torch::Tensor &r0,Surface_mesh& surface_mesh);
    torch::Tensor projectOnSphere(torch::Tensor&c,torch::Tensor&r,torch::Tensor&tp);
    MyCGAL::Primitives::BVHAccel<double>* constructBVH(SlabMesh*mesh,at::Tensor vn,torch::Tensor rn);
};
class Sphere2Boundary
{
public:
    Sphere2Boundary(MyCGAL::Primitives::BVHAccel<double>* tree):tree(tree)
    {

    }
    torch::Tensor forward(SlabMesh *mesh,at::Tensor &v0,at::Tensor &r0);
private:
    MyCGAL::Primitives::BVHAccel<double>* tree;
};

class DistToBoundaryLoss
{
public:
    DistToBoundaryLoss(Tree* tree/*MyCGAL::Primitives::BVHAccel<double>* tree*/):tree(tree){
        k=5;
    }
    void searchMatToSurfaceKNearest(SlabMesh *mesh);
    torch::Tensor forward(SlabMesh *mesh,std::map<face_descriptor,Vector> fnormals,std::map<vertex_descriptor,Vector> vnormals,at::Tensor &v0,at::Tensor &r0);
private:
    Tree* tree;
    int k=3;
    //MyCGAL::Primitives::BVHAccel<double>* tree;
    std::map<int,std::vector<MyCGAL::Primitives::Triangle<double>*>> matToSurfaceKNearest;
};

#endif
