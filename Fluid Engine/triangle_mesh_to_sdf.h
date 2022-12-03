#ifndef INCLUDE_JET_TRIANGLE_MESH_TO_SDF_H_
#define INCLUDE_JET_TRIANGLE_MESH_TO_SDF_H_

#include "scalar_grid3.h"
#include "triangle_mesh3.h"

#include "pch.h"

#include "array3.h"
#include "array_utils.h"

#include <algorithm>
#include <vector>

namespace jet {

    //!
    //! \brief Generates signed-distance field out of given triangle mesh.
    //!
    //! This function generates signed-distance field from a triangle mesh. The sign
    //! is determined by TriangleMesh3::isInside (negative means inside).
    //!
    //! \warning Parameter \p exactBand is no longer used and will be deprecated in
    //! next release (v2.x).
    //!
    //! \param[in]      mesh      The mesh.
    //! \param[in,out]  sdf       The output signed-distance field.
    //! \param[in]      exactBand This parameter is no longer used.
    //!
    void triangleMeshToSdf(const TriangleMesh3& mesh, ScalarGrid3* sdf,
        const unsigned int exactBand = 1);

    void triangleMeshToSdf(const TriangleMesh3& mesh, ScalarGrid3* sdf,
        const unsigned int) {
        Size3 size = sdf->dataSize();
        if (size.x * size.y * size.z == 0) {
            return;
        }

        const auto pos = sdf->dataPosition();
        mesh.updateQueryEngine();
        sdf->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            const Vector3D p = pos(i, j, k);
            const double d = mesh.closestDistance(p);
            const double sd = mesh.isInside(p) ? -d : d;
            (*sdf)(i, j, k) = sd;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_TRIANGLE_MESH_TO_SDF_H_