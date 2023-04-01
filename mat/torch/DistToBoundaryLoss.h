#ifndef DISTTOBOUNDARYLOSS_H
#define DISTTOBOUNDARYLOSS_H
#undef slots
#include <torch/torch.h>
#define slots Q_SLOTS
#include "../PrimMesh.h"
class DistToBoundaryLoss
{
public:
    DistToBoundaryLoss(Tree* tree):tree(tree){

    }
    torch::Tensor forward(at::Tensor &v0,at::Tensor &r0);
private:
    Tree* tree;
};

#endif // DISTTOBOUNDARYLOSS_H
