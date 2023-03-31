INCLUDEPATH+=$$PWD
HEADERS += \
    $$PWD/GeometryObjects/GeometryObjects.h \
    $$PWD/LinearAlgebra/Wm4Math.h \
    $$PWD/LinearAlgebra/Wm4Matrix.h \
    $$PWD/LinearAlgebra/Wm4Vector.h \
    $$PWD/mat/Mesh.h \
    $$PWD/mat/NonManifoldMesh/nonmanifoldmesh.h \
    $$PWD/mat/PrimMesh.h \
    $$PWD/mat/SlabMesh.h \
    $$PWD/mat/threedimensionalshape.h \
    $$PWD/mat/torch/DistToBoundaryLoss.h

SOURCES += \
    $$PWD/GeometryObjects/GeometryObjects.cpp \
    $$PWD/LinearAlgebra/Wm4Math.cpp \
    $$PWD/LinearAlgebra/Wm4Matrix.cpp \
    $$PWD/LinearAlgebra/Wm4Vector.cpp \
    $$PWD/mat/Mesh.cpp \
    $$PWD/mat/NonManifoldMesh/nonmanifoldmesh.cpp \
    $$PWD/mat/PrimMesh.cpp \
    $$PWD/mat/SlabMesh.cpp \
    $$PWD/mat/threedimensionalshape.cpp \
    $$PWD/mat/torch/DistToBoundaryLoss.cpp
