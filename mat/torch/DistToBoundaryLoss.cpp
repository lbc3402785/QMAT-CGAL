#include "DistToBoundaryLoss.h"

at::Tensor DistToBoundaryLoss::forward(at::Tensor &v0,at::Tensor &r0)
{
    int64_t num=v0.sizes()[0];
    torch::Tensor sum=torch::zeros({1});
    #pragma omp parallel for
    for(int64_t i=0;i<num;i++){
        Point query(v0[i][0].item<double>(),v0[i][1].item<double>(),v0[i][2].item<double>());//item不会求导数

        Point_and_primitive_id pp=tree->closest_point_and_primitive(query);
        Point closest_point = pp.first;
        //Vector v=(closest_point-query);
        torch::Tensor q=torch::tensor({closest_point.x(),closest_point.y(),closest_point.z()});
        torch::Tensor v=q-v0[i];
        #pragma omp critical
        {
        sum+=(v.norm()-r0[i]).square();
        }
    }
    return sum/num;
}
