#include "DistToBoundaryLoss.h"

at::Tensor DistToBoundaryLoss::forward(at::Tensor &m0)
{
    int64_t num=m0.sizes()[0];
    torch::Tensor sum=torch::zeros({1});
    //#pragma omp parallel for
    for(int64_t i=0;i<num;i++){
        Point query(m0[i][0].item<double>(),m0[i][1].item<double>(),m0[i][2].item<double>());

        Point_and_primitive_id pp=tree->closest_point_and_primitive(query);
        Point closest_point = pp.first;
        //Vector v=(closest_point-query);
        torch::Tensor c=torch::tensor({closest_point.x(),closest_point.y(),closest_point.z()});
        torch::Tensor q=m0[i].index_select(0,torch::tensor({0,1,2},torch::kLong));
        torch::Tensor v=c-q;
        //#pragma omp critical
        {
        sum+=(v.norm()-m0[i][3]).square();
        }
    }
    return sum/num;
}
