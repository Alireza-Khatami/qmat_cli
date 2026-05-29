#include "QmatVisualizer.h"

#ifdef QMAT_WITH_POLYSCOPE

#include "SlabMesh.h"

void QmatVisualizer::Setup(SlabMesh& sm)
{
    SetupSimplificationViewer(sm, vs_);
}

#endif  // QMAT_WITH_POLYSCOPE
