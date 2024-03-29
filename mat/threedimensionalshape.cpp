#include "threedimensionalshape.h"
#include <map>
ThreeDimensionalShape::ThreeDimensionalShape():slab_initial(false)
{

}

void ThreeDimensionalShape::ComputeInputNMM()
{
    input_nmm.numVertices = 0;
    input_nmm.numEdges = 0;
    input_nmm.numFaces = 0;

    //
    input_nmm.vertices.clear();
    input_nmm.edges.clear();
    input_nmm.faces.clear();
    Triangulation * pt = &(input.dt);

    num_vor_v = 0;
    num_vor_e = 0;
    num_vor_f = 0;
    input_nmm.maxRadius=0;
    input_nmm.avgRadius=0;
    double sumRadius;
    double len[4];
    len[0] = input.m_max[0] - input.m_min[0];
    len[1] = input.m_max[1] - input.m_min[1];
    len[2] = input.m_max[2] - input.m_min[2];
    len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
    input_nmm.diameter = len[3];

    for(Finite_vertices_iterator_t fvi = pt->finite_vertices_begin(); fvi != pt->finite_vertices_end(); fvi ++){
        SamplePoint sample(fvi->point()[0], fvi->point()[1], fvi->point()[2]);
        sample.tag=fvi->info().tag;
        input_nmm.BoundaryPoints.push_back(sample);//可能有一部分点重合，所有这里没有包括所有网格点
    }
    int mas_vertex_count(0);
    for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci ++)
    {
        if(fci->info().inside == false)
        {
            fci->info().tag = -1;
            continue;
        }
        int onSurfaceCount = facetCount(fci);//统计四面体有几个边界面

        fci->info().tag = mas_vertex_count ++;
        Bool_VertexPointer bvp;//中轴点（四面体外接球球心）
        bvp.first = true;
        bvp.second = new NonManifoldMesh_Vertex;
        (*bvp.second).isShow=false;
        Triangulation::Tetrahedron cell=pt->tetrahedron(fci);
        double circumRadius=pt->TetCircumRadius(cell);
        Point_t circumCenter=CGAL::circumcenter(cell);
        double scribedRadius;
        Point_t scribed=pt->inscribeSphereCenter(pt->tetrahedron(fci),&scribedRadius);
        CGAL::Bounded_side side=cell.bounded_side(circumCenter);


        fci->info().faceCount=onSurfaceCount;
        Point center(circumCenter[0],circumCenter[1],circumCenter[2]);
        switch(onSurfaceCount){
        case 0:{
            switch(side){
            case CGAL::ON_UNBOUNDED_SIDE:{
                (*bvp.second).sphere.center = to_wm4(circumCenter);
                (*bvp.second).sphere.radius =std::min(circumRadius,distCenterToBoundary(circumCenter));
                //                (*bvp.second).sphere.center = to_wm4(scribed);
                //                (*bvp.second).sphere.radius =std::min(scribedRadius,distCenterToBoundary(scribed));
                (*bvp.second).isShow=true;
                //bvp.first=false;
                break;
            }
            case CGAL::ON_BOUNDARY:
            case CGAL::ON_BOUNDED_SIDE:
                (*bvp.second).sphere.center = to_wm4(circumCenter);
                (*bvp.second).sphere.radius =std::min(circumRadius,std::min(distCenterToBoundary(circumCenter),distCenterToBoundary(fci,center)));/*distCenterToBoundary(fci,center);*/
                //std::cout<<"wawa"<<std::endl;
                break;

            }
            //std::cout<<"0 facet"<<std::endl;
            break;
        }
        case 1:
            (*bvp.second).sphere.center = to_wm4(circumCenter);
            (*bvp.second).sphere.radius =std::min(circumRadius,distCenterToBoundary(fci,center));
            break;
        case 2:{
            (*bvp.second).sphere.center = to_wm4(circumCenter);
            (*bvp.second).sphere.radius =std::min(circumRadius,distCenterToBoundary(fci,center));
            break;
        }
        case 3:{
            (*bvp.second).sphere.center = to_wm4(circumCenter);
            (*bvp.second).sphere.radius =std::min(circumRadius,distCenterToBoundary(fci,center));
            break;
        }
        case 4:
            (*bvp.second).sphere.center = to_wm4(circumCenter);
            (*bvp.second).sphere.radius =std::min(circumRadius,distCenterToBoundary(fci,center));
            break;
        default:
            std::cerr<<fci->info().tag<<" error"<<std::endl;
        }

        //double inscribedRadius;
        //Point_t inscribed =pt->inscribeSphereCenter(pt->tetrahedron(fci),&inscribedRadius);
        (*bvp.second).is_pole = fci->info().is_pole;
        for (unsigned k = 0; k < 4; k++) {
            (*bvp.second).bplist.insert(fci->vertex(k)->info().id);//关联样本点id
        }
        //(*bvp.second).sphere.radius = pt->TetCircumRadius(pt->tetrahedron(fci));

        if((*bvp.second).sphere.radius>input_nmm.maxRadius){
            input_nmm.maxRadius=(*bvp.second).sphere.radius;
        }
        sumRadius+=(*bvp.second).sphere.radius;
        input_nmm.vertices.push_back(bvp);
        input_nmm.numVertices++;
        num_vor_v++;
    }
    input_nmm.avgRadius=sumRadius/num_vor_v;
    for(Finite_facets_iterator_t ffi = pt->finite_facets_begin(); ffi != pt->finite_facets_end(); ffi ++)//两个对偶面入射四面体的球心构成一条中轴边
    {
        Triangulation::Object o = pt->dual(*ffi);
        if(const Triangulation::Segment *s = CGAL::object_cast<Triangulation::Segment>(&o))
        {
            if( (ffi->first->info().inside == false) || (pt->mirror_facet(*ffi).first->info().inside == false) )
                continue;
            if(!input_nmm.vertices[ffi->first->info().tag].first ||
                    !input_nmm.vertices[pt->mirror_facet(*ffi).first->info().tag].first){//add by me
                continue;
            }
            Bool_EdgePointer bep;//中轴边
            bep.first = true;
            bep.second = new NonManifoldMesh_Edge;
            (*bep.second).vertices_.first = ffi->first->info().tag;
            (*bep.second).vertices_.second = pt->mirror_facet(*ffi).first->info().tag;
            (*input_nmm.vertices[ffi->first->info().tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
            (*input_nmm.vertices[pt->mirror_facet(*ffi).first->info().tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
            input_nmm.edges.push_back(bep);
            input_nmm.numEdges ++;
            num_vor_e ++;
        }
    }
    for(Finite_edges_iterator_t fei = pt->finite_edges_begin(); fei != pt->finite_edges_end(); fei ++)
    {
        //        if(-1==fei->first->info().tag){//add by me
        //            continue;
        //        }
        bool all_finite_inside = true;
        std::vector<Cell_handle_t> vec_ch;
        Cell_circulator_t cc = pt->incident_cells(*fei);//边的所有入射四面体
        do
        {
            if(pt->is_infinite(cc))
                all_finite_inside = false;
            else if(cc->info().inside == false)
                all_finite_inside = false;
            //            else if (cc->info().tag == -1)//add by me
            //                all_finite_inside = false;
            vec_ch.push_back(cc++);
        }while(cc != pt->incident_cells(*fei));
        if(!all_finite_inside)
            continue;
        for(unsigned k = 2; k < vec_ch.size() - 1; k ++)//0--1,1--2,2--3,...,(k-2)--(k-1),(k-1)--0 这些边是两个对偶面入射四面体的球心，在面的循环中已经加入
        {
            if(!input_nmm.vertices[vec_ch[0]->info().tag].first ||
                    !input_nmm.vertices[vec_ch[k]->info().tag].first){//add by me
                continue;
            }
            Bool_EdgePointer bep;
            bep.first = true;
            bep.second = new NonManifoldMesh_Edge;
            (*bep.second).vertices_.first = vec_ch[0]->info().tag;
            (*bep.second).vertices_.second = vec_ch[k]->info().tag;

            (*input_nmm.vertices[vec_ch[0]->info().tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
            (*input_nmm.vertices[vec_ch[k]->info().tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
            input_nmm.edges.push_back(bep);
            input_nmm.numEdges ++;
        }
        for(unsigned k = 1; k < vec_ch.size() - 1; k ++)
        {

            unsigned vid[3];
            vid[0] = vec_ch[0]->info().tag;
            vid[1] = vec_ch[k]->info().tag;
            vid[2] = vec_ch[k+1]->info().tag;

            unsigned eid[3];
            if (!input_nmm.Edge(vid[0], vid[1], eid[0])) {//是否存在这条中轴边
                continue;
            }
            if (!input_nmm.Edge(vid[0], vid[2], eid[1])) {
                continue;
            }
            if (!input_nmm.Edge(vid[1], vid[2], eid[2])) {
                continue;
            }
            Bool_FacePointer bfp;//中轴三角面
            bfp.first = true;
            bfp.second = new NonManifoldMesh_Face;
            (*bfp.second).edges_.insert(eid[0]);//中轴三角面关联的中轴边
            (*bfp.second).edges_.insert(eid[1]);
            (*bfp.second).edges_.insert(eid[2]);

            (*bfp.second).vertices_.insert(vec_ch[0]->info().tag);
            (*bfp.second).vertices_.insert(vec_ch[k]->info().tag);
            (*bfp.second).vertices_.insert(vec_ch[k+1]->info().tag);//中轴三角面的中轴点序号

            input_nmm.vertices[vid[0]].second->faces_.insert(input_nmm.faces.size());
            input_nmm.vertices[vid[1]].second->faces_.insert(input_nmm.faces.size());
            input_nmm.vertices[vid[2]].second->faces_.insert(input_nmm.faces.size());//中轴点关联的中轴三角面序号
            input_nmm.edges[eid[0]].second->faces_.insert(input_nmm.faces.size());
            input_nmm.edges[eid[1]].second->faces_.insert(input_nmm.faces.size());
            input_nmm.edges[eid[2]].second->faces_.insert(input_nmm.faces.size());//中轴边关联的中轴三角面序号
            input_nmm.faces.push_back(bfp);
            input_nmm.numFaces ++;
            num_vor_f ++;
        }
    }
    //=========================sharp point/edge/face=========================
    bool addSharpFeature=false;
    if(addSharpFeature){
        std::map<int,int> shapePointIndexMap;
        for(Finite_vertices_iterator_t fvi = pt->finite_vertices_begin(); fvi != pt->finite_vertices_end(); fvi ++){
            if(input.sharpPointMap.count(fvi->info().id)>0){//添加中轴点，冗余的孤立点事后会被删除
                Point_t p=fvi->point();
                Bool_VertexPointer bvp;//中轴点
                bvp.first = true;
                bvp.second = new NonManifoldMesh_Vertex;
                (*bvp.second).sphere.center=Vector3d(p.x(),p.y(),p.z());
                (*bvp.second).sphere.radius=0;
                shapePointIndexMap[fvi->info().id]=input_nmm.numVertices;
                input_nmm.vertices.push_back(bvp);
                input_nmm.numVertices++;
            }
        }
        std::map<std::pair<int,int>,int> shapeEdgeIndexMap;
        std::map<std::tuple<int,int,int>,int> shapeFaceIndexMap;
        for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci ++)
        {
            if(fci->info().inside == false)
            {
                fci->info().tag = -1;
                continue;
            }
            int i0=fci->vertex(0)->info().id;
            int i1=fci->vertex(1)->info().id;
            int i2=fci->vertex(2)->info().id;
            int i3=fci->vertex(3)->info().id;
            std::vector<int> tmp={i0,i1,i2,i3};
            std::sort(tmp.begin(),tmp.end());
            i0=tmp[0];
            i1=tmp[1];
            i2=tmp[2];
            i3=tmp[3];
            if(input.sharpEdgeMap.count(std::make_pair(i0,i1))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i0,i1),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
            if(input.sharpEdgeMap.count(std::make_pair(i0,i2))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i0,i2),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
            if(input.sharpEdgeMap.count(std::make_pair(i0,i3))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i0,i3),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
            if(input.sharpEdgeMap.count(std::make_pair(i1,i2))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i1,i2),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
            if(input.sharpEdgeMap.count(std::make_pair(i1,i3))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i1,i3),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
            if(input.sharpEdgeMap.count(std::make_pair(i2,i3))>0){
                createSharpSlabEdgeAndFace(fci->info().tag,std::make_pair(i2,i3),shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }


            std::tuple<int,int,int> face012=std::make_tuple(i0,i1,i2);
            if(input.sharpFaceMap.count(face012)>0){
                createSharpSlabTetrahedron(fci->info().tag,face012,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }

            std::tuple<int,int,int> face013=std::make_tuple(i0,i1,i3);
            if(input.sharpFaceMap.count(face013)>0){
                createSharpSlabTetrahedron(fci->info().tag,face013,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }

            std::tuple<int,int,int> face023=std::make_tuple(i0,i2,i3);
            if(input.sharpFaceMap.count(face023)>0){
                createSharpSlabTetrahedron(fci->info().tag,face023,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }

            std::tuple<int,int,int> face123=std::make_tuple(i1,i2,i3);
            if(input.sharpFaceMap.count(face123)>0){
                createSharpSlabTetrahedron(fci->info().tag,face123,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
            }
        }
    }
    //else{
    //    for(Finite_cells_iterator_t fci = pt->finite_cells_begin(); fci != pt->finite_cells_end(); fci ++)
    //    {
    //        if(fci->info().inside == false)
    //        {
    //            fci->info().tag = -1;
    //            continue;
    //        }
    //        int i0=fci->vertex(0)->info().id;
    //        int i1=fci->vertex(1)->info().id;
    //        int i2=fci->vertex(2)->info().id;
    //        int i3=fci->vertex(3)->info().id;
    //        std::vector<int> tmp={i0,i1,i2,i3};
    //        std::sort(tmp.begin(),tmp.end());
    //        i0=tmp[0];
    //        i1=tmp[1];
    //        i2=tmp[2];
    //        i3=tmp[3];



    //        std::pair<int,int> edge01=std::make_pair(i0,i1);
    //        if(input.sharpEdgeMap.count(edge01)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }

    //        std::pair<int,int> edge02=std::make_pair(i0,i2);
    //        if(input.sharpEdgeMap.count(edge02)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }

    //        std::pair<int,int> edge03=std::make_pair(i0,i3);
    //        if(input.sharpEdgeMap.count(edge03)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }

    //        std::pair<int,int> edge12=std::make_pair(i1,i2);
    //        if(input.sharpEdgeMap.count(edge12)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }

    //        std::pair<int,int> edge13=std::make_pair(i1,i3);
    //        if(input.sharpEdgeMap.count(edge13)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }

    //        std::pair<int,int> edge23=std::make_pair(i2,i3);
    //        if(input.sharpEdgeMap.count(edge23)>0){
    //            //input_nmm.vertices[fci->info().tag].second->saved_vertex=true;
    //            continue;
    //        }
    //    }
    //}




    input_nmm.Export(input_nmm.meshname);

    input_nmm.numVertices = 0;
    input_nmm.numEdges = 0;
    input_nmm.numFaces = 0;
}
void ThreeDimensionalShape::LoadInputNMM(std::string fname){
    std::ifstream mastream(fname.c_str());
    //    NonManifoldMesh newinputnmm;
    //    newinputnmm.numVertices = 0;
    //    newinputnmm.numEdges = 0;
    //    newinputnmm.numFaces = 0;
    int nv, ne, nf;
    mastream >> nv >> ne >> nf;

    // slab mesh
    /*slab_mesh.numVertices = 0;
    slab_mesh.numEdges = 0;
    slab_mesh.numFaces = 0;*/

    double len[4];
    len[0] = input.m_max[0] - input.m_min[0];
    len[1] = input.m_max[1] - input.m_min[1];
    len[2] = input.m_max[2] - input.m_min[2];
    len[3] = sqrt(len[0]*len[0]+len[1]*len[1]+len[2]*len[2]);
    //newinputnmm.diameter = len[3];
    slab_mesh.bound_weight = 0.1;

    //    for(unsigned i = 0; i < input.pVertexList.size(); i ++)
    //        newinputnmm.BoundaryPoints.push_back(SamplePoint(
    //                                                 input.pVertexList[i]->point()[0],
    //                                             input.pVertexList[i]->point()[1],
    //                input.pVertexList[i]->point()[2]
    //                ));

    for(unsigned i = 0; i < nv; i ++)
    {
        char ch;
        double x,y,z,r;
        mastream >> ch >> x >> y >> z >> r;

        // handle the slab mesh
        Bool_SlabVertexPointer bsvp2;//用于渲染的中轴点
        bsvp2.first = true;
        bsvp2.second = new SlabVertex;
        (*bsvp2.second).sphere.center[0] = x /*/ input.bb_diagonal_length*/;
        (*bsvp2.second).sphere.center[1] = y /*/ input.bb_diagonal_length*/;
        (*bsvp2.second).sphere.center[2] = z /*/ input.bb_diagonal_length*/;
        (*bsvp2.second).sphere.radius = r /*/ input.bb_diagonal_length*/;
        (*bsvp2.second).index = slab_mesh.vertices.size();
        slab_mesh.vertices.push_back(bsvp2);
        //slab_mesh.numVertices ++;
    }

    for(unsigned i = 0; i < ne; i ++)
    {
        char ch;
        unsigned ver[2];
        mastream >> ch;
        mastream >> ver[0];
        mastream >> ver[1];

        // handle the slab mesh
        Bool_SlabEdgePointer bsep2;
        bsep2.first = true;
        bsep2.second = new SlabEdge;
        (*bsep2.second).vertices_.first = ver[0];
        (*bsep2.second).vertices_.second = ver[1];
        (*slab_mesh.vertices[(*bsep2.second).vertices_.first].second).edges_.insert(slab_mesh.edges.size());
        (*slab_mesh.vertices[(*bsep2.second).vertices_.second].second).edges_.insert(slab_mesh.edges.size());
        (*bsep2.second).index = slab_mesh.edges.size();
        slab_mesh.edges.push_back(bsep2);
        //slab_mesh.numEdges ++;
    }

    for(unsigned i = 0; i < nf; i ++)
    {
        char ch;
        unsigned vid[3];
        unsigned eid[3];
        mastream >> ch >> vid[0] >> vid[1] >> vid[2];

        // handle the slab mesh
        Bool_SlabFacePointer bsfp2;
        bsfp2.first = true;
        bsfp2.second = new SlabFace;
        (*bsfp2.second).vertices_.insert(vid[0]);
        (*bsfp2.second).vertices_.insert(vid[1]);
        (*bsfp2.second).vertices_.insert(vid[2]);
        if (slab_mesh.Edge(vid[0], vid[1], eid[0])) {
            (*bsfp2.second).edges_.insert(eid[0]);
        }
        else {
            std::cerr << "no found" << std::endl;
        }
        if (slab_mesh.Edge(vid[0], vid[2], eid[1])) {
            (*bsfp2.second).edges_.insert(eid[1]);
        }
        else {
            std::cerr << "no found" << std::endl;
        }
        if (slab_mesh.Edge(vid[1], vid[2], eid[2])) {
            (*bsfp2.second).edges_.insert(eid[2]);
        }
        else {
            std::cerr << "no found" << std::endl;
        }
        (*bsfp2.second).index = slab_mesh.faces.size();
        slab_mesh.vertices[vid[0]].second->faces_.insert(slab_mesh.faces.size());
        //slab_mesh.vertices[vid[0]].second->related_face += 2;
        slab_mesh.vertices[vid[1]].second->faces_.insert(slab_mesh.faces.size());
        //slab_mesh.vertices[vid[1]].second->related_face += 2;
        slab_mesh.vertices[vid[2]].second->faces_.insert(slab_mesh.faces.size());
        //slab_mesh.vertices[vid[2]].second->related_face += 2;
        slab_mesh.edges[eid[0]].second->faces_.insert(slab_mesh.faces.size());
        slab_mesh.edges[eid[1]].second->faces_.insert(slab_mesh.faces.size());
        slab_mesh.edges[eid[2]].second->faces_.insert(slab_mesh.faces.size());
        slab_mesh.faces.push_back(bsfp2);
        //slab_mesh.numFaces++;
    }
    std::string saved_vertex;
    mastream>>saved_vertex;
    if(saved_vertex.rfind("#saved_vertex", 0)==0){
        char ch;
        for (unsigned i = 0; i < slab_mesh.vertices.size(); i++){
            mastream >> ch >>slab_mesh.vertices[i].second->saved_vertex;
        }
    }
    std::string bplist;
    mastream>>bplist;
    if(bplist.rfind("#bplist", 0)==0){
        for(unsigned i = 0; i < nv; i ++)
        {
            unsigned size;
            mastream>>size;
            unsigned tag;
            for(unsigned k=0;k<size;k++){
                mastream>>tag;
                slab_mesh.vertices[i].second->bplist.emplace(tag);
            }
        }
    }
    //newinputnmm.ComputeFacesNormal();
    //newinputnmm.ComputeFacesCentroid();
    //newinputnmm.ComputeFacesSimpleTriangles();
    //newinputnmm.ComputeEdgesCone();
    //input_nmm = newinputnmm;

    /*slab_mesh.iniNumVertices = slab_mesh.numVertices;
    slab_mesh.iniNumEdges = slab_mesh.numEdges;
    slab_mesh.iniNumFaces = slab_mesh.numFaces;*/

    slab_mesh.update();
}

long ThreeDimensionalShape::PrepareSimplifySlabMesh()
{
    slab_mesh.clear();//清空数据
    long startt = clock();
    InitialSlabMesh();//计算cost需要的矩阵
    slab_mesh.initCollapseQueue();//计算所有边的cost
    slab_mesh.checkSkeletonAndTube();//保护骨骼
    long endt = clock();
    return endt - startt;
}
void ThreeDimensionalShape::InitialSlabMesh()
{
    // handle each face
    for(unsigned i = 0; i < slab_mesh.vertices.size(); i++)
    {
        if(!slab_mesh.vertices[i].first)
            continue;

        SlabVertex sv = *slab_mesh.vertices[i].second;
        std::set<unsigned> fset = sv.faces_;
        Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(), sv.sphere.center.Z(), sv.sphere.radius);

        for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
        {
            SlabFace sf = *slab_mesh.faces[*si].second;

            if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) ||
                    sf.st[1].normal == Vector3d(0., 0., 0.))
                continue;

            Vector4d normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(), sf.st[0].normal.Z(), 1.0);//n1
            Vector4d normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(), sf.st[1].normal.Z(), 1.0);//n2

            // compute the matrix of A
            Matrix4d temp_A1, temp_A2;
            temp_A1.MakeTensorProduct(normal1, normal1);//temp_A1=(n_1,1) * (n_1,1)^T
            temp_A2.MakeTensorProduct(normal2, normal2);//temp_A2=(n_2,1) * (n_2,1)^T
            temp_A1 *= 2.0;//temp_A1=2*(n_1,1) * (n_1,1)^T
            temp_A2 *= 2.0;//temp_A2=2*(n_2,1) * (n_2,1)^T

            // compute the matrix of b
            double normal_mul_point1 = normal1.Dot(C1);//(n_1,1)^T*m
            double normal_mul_point2 = normal2.Dot(C1);//(n_2,1)^T*m
            Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;//2*(n_1,1)*(n_1,1)^T*m
            Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;//2*(n_2,1)*(n_2,1)^T*m

            //compute c
            double temp_c1 = normal_mul_point1 * normal_mul_point1;//m^T*(n_1,1)*(n_1,1)^T*m
            double temp_c2 = normal_mul_point2 * normal_mul_point2;//m^T*(n_2,1)*(n_2,1)^T*m

            slab_mesh.vertices[i].second->slab_A += temp_A1;
            slab_mesh.vertices[i].second->slab_A += temp_A2;
            slab_mesh.vertices[i].second->slab_b += temp_b1;
            slab_mesh.vertices[i].second->slab_b += temp_b2;
            slab_mesh.vertices[i].second->slab_c += temp_c1;
            slab_mesh.vertices[i].second->slab_c += temp_c2;

            slab_mesh.vertices[i].second->related_face += 2;
        }
    }

    switch(slab_mesh.preserve_boundary_method)//default 0
    {
    case 1 :
        slab_mesh.PreservBoundaryMethodOne();
        break;
    case 2 :
        //slab_mesh.PreservBoundaryMethodTwo();
        break;
    case 3 :
        slab_mesh.PreservBoundaryMethodThree();
        break;
    default:
        slab_mesh.PreservBoundaryMethodFour();
        break;
    }

}
double ThreeDimensionalShape::NearestPoint(Vector3d point, unsigned vid)
{
    set<unsigned> faces = slab_mesh.vertices[vid].second->faces_;
    set<unsigned> edges = slab_mesh.vertices[vid].second->edges_;
    double mind = DBL_MAX;

    // calculation of related faces
    for (set<unsigned>::iterator si = faces.begin(); si != faces.end(); si++)
    {
        if (!slab_mesh.faces[*si].first)
            continue;
        SlabFace sf = *slab_mesh.faces[*si].second;
        if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) ||
                sf.st[1].normal == Vector3d(0., 0., 0.))
            continue;

        Vector3d v[3];
        Vector3d tfp;
        double td;
        for (int i = 0; i < 2; i++)
        {
            v[0] = sf.st[i].v[0];
            v[1] = sf.st[i].v[1];
            v[2] = sf.st[i].v[2];
            ProjectOntoTriangle(point, v[0], v[1], v[2], tfp, td);

            if(td < mind)	mind = td;
        }
    }

    // calculation of related edges
    for (set<unsigned>::iterator si = edges.begin(); si != edges.end(); si++)
    {
        if (!slab_mesh.edges[*si].first)
            continue;
        SlabEdge se = *slab_mesh.edges[*si].second;
        if (se.valid_cone == false)
            continue;

        Vector3d tfp;
        double td;
        se.cone.ProjectOntoCone(point, tfp, td);
        if (td < 0)
            td = -td;

        if(td < mind)	mind = td;
    }

    return mind;
}

void ThreeDimensionalShape::ComputeHausdorffDistance()
{
#if 0
    long start_time = clock();

    ma_qem_mesh.maxhausdorff_distance = 0.;
    for(unsigned i = 0; i < ma_qem_mesh.faces.size(); i ++)
    {
        Vector3d ve[8];
        if (ma_qem_mesh.faces[i].second->valid_st == false ||
                ma_qem_mesh.faces[i].second->st[0].normal == Vector3d(0., 0., 0.) ||
                ma_qem_mesh.faces[i].second->st[1].normal == Vector3d(0., 0., 0.))
            continue;

        ve[0] = ma_qem_mesh.faces[i].second->st[0].v[0];
        ve[1] = ma_qem_mesh.faces[i].second->st[0].v[1];
        ve[2] = ma_qem_mesh.faces[i].second->st[0].v[2];
        ve[3] = ma_qem_mesh.faces[i].second->st[1].v[0];
        ve[4] = ma_qem_mesh.faces[i].second->st[1].v[1];
        ve[5] = ma_qem_mesh.faces[i].second->st[1].v[2];
        ve[6] = (ve[0] + ve[1] +ve[2]) / 3.0;
        ve[7] = (ve[3] + ve[4] +ve[5]) / 3.0;
        double face_haus = 0.;
        for (int j = 0; j < 8; j++)
        {
            Vector3d fp = input.NearestVertex(ve[j]);
            double len = (ve[j] - fp).Length();
            face_haus = max(len, face_haus);
        }
        ma_qem_mesh.faces[i].second->hausdorff_dist = face_haus;

        ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance,face_haus);
    }

    long end_time = clock();
    long result = end_time - start_time;
#endif

    //ma_qem_mesh.maxhausdorff_distance = 0.;
    slab_mesh.maxhausdorff_distance = 0;
    double sumhausdorff_distance = 0;
    for (unsigned i = 0; i < input.pVertexList.size(); i++)
    {
        double min_dis = DBL_MAX;
        unsigned min_index = -1;
        Vector3d bou_ver(input.pVertexList[i]->point()[0], input.pVertexList[i]->point()[1], input.pVertexList[i]->point()[2]);
        bou_ver /=  input.bb_diagonal_length;

        for (unsigned j = 0; j < slab_mesh.vertices.size(); j++)
        {
            Sphere ma_ver = slab_mesh.vertices[j].second->sphere;
            double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
            //if (temp_length >= 0 && temp_length < min_dis)
            if (temp_length < min_dis)
            {
                min_dis = temp_length;
                min_index = j;
            }

            //double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
            //if (temp_near_dis < min_dis)
            //{
            //	min_dis = temp_near_dis;
            //	min_index = j;
            //}

        }

        //// 为何这里得出的结果比前面得出的结果还要小？
        //double nearest_dis = NearestPoint(bou_ver, min_index);
        //min_dis = min(nearest_dis, min_dis);

        //sumhausdorff_distance += min_dis;


        if (min_index != -1)
        {
            double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
            min_dis = min(temp_near_dis, min_dis);

            sumhausdorff_distance += min_dis;

            //ma_qem_mesh.vertices[min_index].second->bplist.push_back(i);
            //ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance, min_dis);

            slab_mesh.vertices[min_index].second->bplist.insert(i);
            slab_mesh.maxhausdorff_distance = max(slab_mesh.maxhausdorff_distance, min_dis);

            //input.pVertexList[i]->vqem_hausdorff_dist = min_dis / input.bb_diagonal_length;
            input.pVertexList[i]->vqem_hausdorff_dist = min_dis;
            input.pVertexList[i]->vqem_hansdorff_index = min_index;

            //input.pVertexList[i]->slab_hausdorff_dist = min_dis / input.bb_diagonal_length;
            input.pVertexList[i]->slab_hausdorff_dist = min_dis;
            input.pVertexList[i]->slab_hansdorff_index = min_index;
        }
    }

    //ma_qem_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
    slab_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
    //ma_qem_mesh.initialhausdorff_distance = ma_qem_mesh.maxhausdorff_distance;
    slab_mesh.initialhausdorff_distance = slab_mesh.maxhausdorff_distance;
}

void ThreeDimensionalShape::PruningSlabMesh()
{
    slab_mesh.ComputeVerticesProperty();

    bool has_boundary_non_pole;
    do
    {
        has_boundary_non_pole = false;
        unsigned vid;
        for(unsigned i = 0; i < slab_mesh.vertices.size(); i ++)
            if(slab_mesh.vertices[i].first)
                if((slab_mesh.vertices[i].second->is_boundary) &&
                        (!slab_mesh.vertices[i].second->is_non_manifold) &&
                        (!slab_mesh.vertices[i].second->is_pole) &&
                        (slab_mesh.vertices[i].second->edges_.size() == 2))
                {
                    vid = i;
                    has_boundary_non_pole = true;
                    break;
                }
        if(has_boundary_non_pole)
            slab_mesh.DeleteVertex(vid);
    }while(has_boundary_non_pole);


    slab_mesh.CleanIsolatedVertices();
    slab_mesh.ComputeEdgesCone();
    slab_mesh.ComputeFacesSimpleTriangles();
    slab_mesh.DistinguishVertexType();
    slab_mesh.computebb();
}

double ThreeDimensionalShape::distCenterToBoundary(Point_t dtPoint)
{
    Locate_type lt;
    int li,lj;
    Cell_handle_t cell=input.dt.locate(dtPoint,lt,li,lj);
    Point center(dtPoint[0],dtPoint[1],dtPoint[2]);
    switch(lt){
    case Triangulation::VERTEX:
        std::cout<<"VERTEX"<<std::endl;
        return 0;
    case Triangulation::EDGE:{
        int vi=cell->vertex(li)->info().id;
        int vj=cell->vertex(lj)->info().id;
        std::pair<int,int> edgeij=std::make_pair(vi<vj?vi:vj,vi<vj?vj:vi);
        if(input.edgeMap.count(edgeij)>0){
            std::cout<<"EDGE"<<std::endl;
            return 0;
        }
        break;
    }
    case Triangulation::FACET:
        if(!cell->neighbor(li)->info().inside){
            std::cout<<"FACET"<<std::endl;
            return 0;
        }
        break;
    case Triangulation::CELL:
        //std::cout<<"CELL"<<std::endl;
        
        break;
    default:
        return 0;
    }
    return distCenterToBoundary(cell, center);
}
FT ThreeDimensionalShape::nearestPointOfLine(Point p,Point a,Point b){
    Vector_3 vap=p-a;
    Vector_3 vab=b-a;
    FT t=vap*vab/(vab*vab);
    return t;
}
int ThreeDimensionalShape::checkPointInTriangle(Point p,Point a,Point b,Point c)
{
    Vector_3 v01=b-a;
    Vector_3 v12=c-b;
    Vector_3 v20=a-c;
    Vector_3 n=CGAL::cross_product(v01,v12);
    FT f0=CGAL::cross_product(p-a,v01)*n;
    FT f1=CGAL::cross_product(p-b,v12)*n;
    FT f2=CGAL::cross_product(p-c,v20)*n;
    if(f0==0||f1==0||f2==0){
        return 0;
    }else if(f0>0&&f1>0&&f2>0){
        return 1;
    }else if(f0<0&&f1<0&&f2<0){
        return 1;
    }else{
        return -1;
    }
}

void ThreeDimensionalShape::createSharpSlabEdgeAndFace(int tag, std::pair<int, int> edge, std::map<int, int> &shapePointIndexMap, std::map<std::pair<int, int>, int> &shapeEdgeIndexMap, std::map<std::tuple<int, int, int>, int> &shapeFaceIndexMap)
{
    int i0=edge.first;
    int i1=edge.second;
    std::vector<int> tmp = { shapePointIndexMap[i0],shapePointIndexMap[i1]};
    std::sort(tmp.begin(), tmp.end());
    i0=tmp[0];
    i1=tmp[1];
    std::pair<int,int> edge01=std::make_pair(i0,i1);
    if(shapeEdgeIndexMap.count(edge01)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = i0;
        (*bep.second).vertices_.second = i1;
        (*input_nmm.vertices[i0].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i1].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edge01]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }
    std::pair<int,int> edget0=std::make_pair(tag,i0);
    if(shapeEdgeIndexMap.count(edget0)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = tag;
        (*bep.second).vertices_.second = i0;
        (*input_nmm.vertices[tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i0].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edget0]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }
    std::pair<int,int> edget1=std::make_pair(tag,i1);
    if(shapeEdgeIndexMap.count(edget1)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = tag;
        (*bep.second).vertices_.second = i1;
        (*input_nmm.vertices[tag].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i1].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edget1]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }
    std::tuple<int,int,int> face=std::make_tuple(tag,i0,i1);
    if(shapeFaceIndexMap.count(face)==0){
        Bool_FacePointer bfp;//中轴三角面
        bfp.first = true;
        bfp.second = new NonManifoldMesh_Face;

        (*bfp.second).vertices_.insert(tag);
        (*bfp.second).vertices_.insert(i0);
        (*bfp.second).vertices_.insert(i1);//中轴三角面的中轴点序号
        unsigned eid[3];
        if (input_nmm.Edge(i0, i1, eid[2])) {
            (*bfp.second).edges_.insert(eid[2]);
        }
        else {
            std::cerr << "wwww" << std::endl;
        }
        if (input_nmm.Edge(tag, i0, eid[0])) {//是否存在这条中轴边
            (*bfp.second).edges_.insert(eid[0]);//中轴三角面关联的中轴边
        }else {
            std::cerr << "wwww" << std::endl;
        }
        if (input_nmm.Edge(tag, i1, eid[1])) {
            (*bfp.second).edges_.insert(eid[1]);
        }else {
            std::cerr << "wwww" << std::endl;
        }
        
        input_nmm.vertices[tag].second->faces_.insert(input_nmm.faces.size());
        input_nmm.vertices[i0].second->faces_.insert(input_nmm.faces.size());
        input_nmm.vertices[i1].second->faces_.insert(input_nmm.faces.size());//中轴点关联的中轴三角面序号
        input_nmm.edges[eid[0]].second->faces_.insert(input_nmm.faces.size());
        input_nmm.edges[eid[1]].second->faces_.insert(input_nmm.faces.size());
        input_nmm.edges[eid[2]].second->faces_.insert(input_nmm.faces.size());//中轴边关联的中轴三角面序号
        input_nmm.faces.push_back(bfp);
        input_nmm.numFaces ++;
    }

}

void ThreeDimensionalShape::createSharpSlabTetrahedron(int tag,std::tuple<int, int, int> face, std::map<int, int> &shapePointIndexMap, std::map<std::pair<int, int>, int> &shapeEdgeIndexMap, std::map<std::tuple<int, int, int>, int> &shapeFaceIndexMap)
{
    int i0 = std::get<0>(face);
    int i1 = std::get<1>(face);
    int i2 = std::get<2>(face);
    std::vector<int> tmp = { shapePointIndexMap[i0],shapePointIndexMap[i1],shapePointIndexMap[i2] };
    std::sort(tmp.begin(), tmp.end());
    i0 = tmp[0];
    i1 = tmp[1];
    i2 = tmp[2];
    face=std::make_tuple(i0,i1,i2);
    if(shapeFaceIndexMap.count(face)==0){
        createSlabFace(face,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
    }

    std::tuple<int,int,int> facet01=std::make_tuple(tag,i0,i1);
    if(shapeFaceIndexMap.count(facet01)==0){
        createSlabFace(facet01,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
    }

    std::tuple<int,int,int> facet02=std::make_tuple(tag,i0,i2);
    if(shapeFaceIndexMap.count(facet02)==0){
        createSlabFace(facet02,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
    }
    std::tuple<int,int,int> facet12=std::make_tuple(tag,i1,i2);
    if(shapeFaceIndexMap.count(facet12)==0){
        createSlabFace(facet12,shapePointIndexMap,shapeEdgeIndexMap,shapeFaceIndexMap);
    }


}



void ThreeDimensionalShape::createSlabFace(std::tuple<int, int, int> face, std::map<int, int> &shapePointIndexMap, std::map<std::pair<int, int>, int> &shapeEdgeIndexMap, std::map<std::tuple<int, int, int>, int> &shapeFaceIndexMap)
{
    int i0=std::get<0>(face);
    int i1=std::get<1>(face);
    int i2=std::get<2>(face);


    std::pair<int,int> edge01=std::make_pair(i0,i1);
    if(shapeEdgeIndexMap.count(edge01)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = i0;
        (*bep.second).vertices_.second = i1;
        (*input_nmm.vertices[i0].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i1].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edge01]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }

    std::pair<int,int> edge02=std::make_pair(i0,i2);
    if(shapeEdgeIndexMap.count(edge02)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = i0;
        (*bep.second).vertices_.second = i2;
        (*input_nmm.vertices[i0].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i2].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edge02]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }

    std::pair<int,int> edge12=std::make_pair(i1,i2);
    if(shapeEdgeIndexMap.count(edge12)==0){
        Bool_EdgePointer bep;
        bep.first = true;
        bep.second = new NonManifoldMesh_Edge;
        (*bep.second).vertices_.first = i1;
        (*bep.second).vertices_.second = i2;
        (*input_nmm.vertices[i1].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        (*input_nmm.vertices[i2].second).edges_.insert(input_nmm.edges.size());//每个中轴点都存储关联的中轴边的序号
        shapeEdgeIndexMap[edge12]=input_nmm.edges.size();
        input_nmm.edges.push_back(bep);
        input_nmm.numEdges ++;
    }
    std::tuple<int, int, int> slabFace = std::make_tuple(i0, i1, i2);
    if(shapeFaceIndexMap.count(slabFace)==0){
        Bool_FacePointer bfp;//中轴三角面
        bfp.first = true;
        bfp.second = new NonManifoldMesh_Face;

        (*bfp.second).vertices_.insert(i0);
        (*bfp.second).vertices_.insert(i1);
        (*bfp.second).vertices_.insert(i2);//中轴三角面的中轴点序号
        unsigned eid[3];
        if (input_nmm.Edge(i0, i1, eid[0])) {//是否存在这条中轴边
            (*bfp.second).edges_.insert(eid[0]);//中轴三角面关联的中轴边
        }else {
            std::cerr << "wwww" << std::endl;
        }
        if (input_nmm.Edge(i0, i2, eid[1])) {
            (*bfp.second).edges_.insert(eid[1]);
        }else {
            std::cerr << "wwww" << std::endl;
        }
        if (input_nmm.Edge(i1, i2, eid[2])) {
            (*bfp.second).edges_.insert(eid[2]);
        }else {
            std::cerr << "wwww" << std::endl;
        }
        input_nmm.vertices[i0].second->faces_.insert(input_nmm.faces.size());
        input_nmm.vertices[i1].second->faces_.insert(input_nmm.faces.size());
        input_nmm.vertices[i2].second->faces_.insert(input_nmm.faces.size());//中轴点关联的中轴三角面序号
        input_nmm.edges[eid[0]].second->faces_.insert(input_nmm.faces.size());
        input_nmm.edges[eid[1]].second->faces_.insert(input_nmm.faces.size());
        input_nmm.edges[eid[2]].second->faces_.insert(input_nmm.faces.size());//中轴边关联的中轴三角面序号
        input_nmm.faces.push_back(bfp);
        input_nmm.numFaces ++;
    }
}

double ThreeDimensionalShape::distCenterToBoundary(const Cell_handle_t &fci,Point center)
{
    Point_and_primitive_id pp=input.tree->closest_point_and_primitive(center);
    Point closest_point = pp.first;
    Surface_mesh::Face_handle fh = pp.second; // closest primitive id
    Vector v=(closest_point-center);
    return std::sqrt(v.x()*v.x()+v.y()*v.y()+v.z()*v.z());
    //    FT minDis=std::numeric_limits<FT>::max();
    //    for(int k=0;k<4;++k){
    //        int vk=fci->vertex(k)->info().id;
    //        Point pk=input.pVertexList[vk]->point();
    //        Halfedge_around_vertex_circulator hav =input.pVertexList[vk]->vertex_begin();
    //        do{
    //            if (hav->facet() == nullptr) {
    //                //std::cerr << "it is a outlier" << std::endl;
    //            }else{
    //                Halfedge_handle h=hav->facet()->halfedge();
    //                Point p0=h->vertex()->point();
    //                Point p1=h->next()->vertex()->point();
    //                Point p2=h->next()->next()->vertex()->point();
    //                Plane plane(p0 ,
    //                            p1,
    //                            p2);
    //                //std::cout<<plane.a()<<"," << plane.b() << "," << plane.c() << ","<< plane.d()<<std::endl;
    //                Point prj= plane.projection(center);

    //                int result=checkPointInTriangle(prj,p0,p1,p2);
    //                if(result==-1){
    //                    FT minTriDis=std::numeric_limits<FT>::max();
    //                    FT t=nearestPointOfLine(center,p0,p1);
    //                    if(t>=0&&t<=1){
    //                        Point target=p0+t*(p1-p0);
    //                        FT tmp=std::sqrt(CGAL::squared_distance(target,center));
    //                        if(tmp<minTriDis){
    //                            minTriDis=tmp;
    //                        }
    //                    }
    //                    t=nearestPointOfLine(center,p1,p2);
    //                    if(t>=0&&t<=1){
    //                        Point target=p1+t*(p2-p1);
    //                        FT tmp=std::sqrt(CGAL::squared_distance(target,center));
    //                        if(tmp<minTriDis){
    //                            minTriDis=tmp;
    //                        }
    //                    }
    //                    t=nearestPointOfLine(center,p2,p0);
    //                    if(t>=0&&t<=1){
    //                        Point target=p2+t*(p0-p2);
    //                        FT tmp=std::sqrt(CGAL::squared_distance(target,center));
    //                        if(tmp<minTriDis){
    //                            minTriDis=tmp;
    //                        }
    //                    }
    //                    FT tmp=std::sqrt(CGAL::squared_distance(p0,center));
    //                    if(tmp<minTriDis){
    //                        minTriDis=tmp;
    //                    }
    //                    tmp=std::sqrt(CGAL::squared_distance(p1,center));
    //                    if(tmp<minTriDis){
    //                        minTriDis=tmp;
    //                    }
    //                    tmp=std::sqrt(CGAL::squared_distance(p2,center));
    //                    if(tmp<minTriDis){
    //                        minTriDis=tmp;
    //                    }
    //                    if(minTriDis<minDis){
    //                        minDis=minTriDis;
    //                    }
    //                }else{
    //                    FT dis=std::sqrt(CGAL::squared_distance(prj,center));
    //                    if(dis<minDis){
    //                        minDis=dis;
    //                    }
    //                }
    //            }

    //            hav ++;
    //        }while(hav !=input.pVertexList[vk]->vertex_begin());
    //    }
    //    //std::cout<<minDis<<std::endl;
    //    return minDis;

}
bool checkEqual(std::vector<int> &f1,std::vector<int> &f2){
    std::sort(f1.begin(),f1.end());
    std::sort(f2.begin(),f2.end());
    for(int k=0;k<f1.size();++k){
        if(f1[k]!=f2[k])return false;
    }
    return true;
}
int ThreeDimensionalShape::facetCount(const Finite_cells_iterator_t &fci)
{
    std::map<Facet_handle,bool> visitMap;
    int count=0;

    for(int k=0;k<4;++k){
        if(!fci->neighbor(k)->info().inside){
            count++;
        }
    }
    //    for(int k=0;k<4;++k){
    //        int vk0=fci->vertex(k)->info().id;
    //        int vk1=fci->vertex((k+1)%4)->info().id;
    //        int vk2=fci->vertex((k+2)%4)->info().id;
    //        std::vector<int> fk={vk0,vk1,vk2};
    //        Halfedge_around_vertex_circulator hav0 =input.pVertexList[vk0]->vertex_begin();
    //        Halfedge_around_vertex_circulator hav1 =input.pVertexList[vk1]->vertex_begin();
    //        Halfedge_around_vertex_circulator hav2 =input.pVertexList[vk2]->vertex_begin();
    //        bool skip=false;
    //        do{
    //            if (hav0->facet() != nullptr) {
    //                Facet_handle f=hav0->facet();
    //                if(visitMap.count(f)==0){
    //                    Halfedge_handle h=hav0->facet()->halfedge();
    //                    int v0=h->vertex()->id;
    //                    int v1=h->next()->vertex()->id;
    //                    int v2=h->next()->next()->vertex()->id;
    //                    std::vector<int> fh={v0,v1,v2};
    //                    if(checkEqual(fh,fk)){
    //                        count++;
    //                        skip=true;
    //                        break;
    //                    }
    //                }
    //                visitMap[f]=true;
    //            }
    //            hav0++;
    //        }while(hav0 !=input.pVertexList[vk0]->vertex_begin());
    //        if(skip)continue;
    //        do{
    //            if (hav1->facet() != nullptr) {
    //                Facet_handle f=hav1->facet();
    //                if(visitMap.count(f)==0){
    //                    Halfedge_handle h=hav1->facet()->halfedge();
    //                    int v0=h->vertex()->id;
    //                    int v1=h->next()->vertex()->id;
    //                    int v2=h->next()->next()->vertex()->id;
    //                    std::vector<int> fh={v0,v1,v2};
    //                    if(checkEqual(fh,fk)){
    //                        count++;
    //                        skip=true;
    //                        break;
    //                    }
    //                }
    //                visitMap[f]=true;
    //            }
    //            hav1++;
    //        }while(hav1 !=input.pVertexList[vk1]->vertex_begin());
    //        if(skip)continue;
    //        do{
    //            if (hav2->facet() != nullptr) {
    //                Facet_handle f=hav2->facet();
    //                if(visitMap.count(f)==0){
    //                    Halfedge_handle h=hav2->facet()->halfedge();
    //                    int v0=h->vertex()->id;
    //                    int v1=h->next()->vertex()->id;
    //                    int v2=h->next()->next()->vertex()->id;
    //                    std::vector<int> fh={v0,v1,v2};
    //                    if(checkEqual(fh,fk)){
    //                        count++;
    //                        break;
    //                    }
    //                }
    //                visitMap[f]=true;
    //            }
    //            hav2++;
    //        }while(hav2 !=input.pVertexList[vk2]->vertex_begin());
    //    }
    return count;
}
