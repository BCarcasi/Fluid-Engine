#ifndef INCLUDE_JET_GRID_BLOCKED_BOUNDARY_CONDITION_SOLVER3_H_
#define INCLUDE_JET_GRID_BLOCKED_BOUNDARY_CONDITION_SOLVER3_H_

#include "grid_fractional_boundary_condition_solver3.h"

#include <memory>
#include "pch.h"
#include "physics_helpers.h"
#include "array_utils.h"
#include "grid_blocked_boundary_condition_solver3.h"
#include "level_set_utils.h"
#include "surface_to_implicit3.h"
#include <algorithm>

namespace jet {

    //!
    //! \brief Blocked 3-D boundary condition solver for grids.
    //!
    //! This class constrains the velocity field by projecting the flow to the
    //! blocked representation of the collider. A collider is rasterized into voxels
    //! and each face of the collider voxels projects the velocity field onto its
    //! face. This implementation should pair up with GridSinglePhasePressureSolver3
    //! since the pressure solver assumes blocked boundary representation as well.
    //!
    class GridBlockedBoundaryConditionSolver3 final
        : public GridFractionalBoundaryConditionSolver3 {
    public:
        //! Default constructor.
        GridBlockedBoundaryConditionSolver3();

        //!
        //! Constrains the velocity field to conform the collider boundary.
        //!
        //! \param velocity Input and output velocity grid.
        //! \param extrapolationDepth Number of inner-collider grid cells that
        //!     velocity will get extrapolated.
        //!
        void constrainVelocity(
            FaceCenteredGrid3* velocity,
            unsigned int extrapolationDepth = 5) override;

        //! Returns the marker which is 1 if occupied by the collider.
        const Array3<char>& marker() const;

    protected:
        //! Invoked when a new collider is set.
        void onColliderUpdated(
            const Size3& gridSize,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin) override;

    private:
        Array3<char> _marker;
    };

    //! Shared pointer type for the GridBlockedBoundaryConditionSolver3.
    typedef std::shared_ptr<GridBlockedBoundaryConditionSolver3>
        GridBlockedBoundaryConditionSolver3Ptr;

    static const char blocked_boundary_kfluid = 1;
    static const char kCollider = 0;

    GridBlockedBoundaryConditionSolver3::GridBlockedBoundaryConditionSolver3() {
    }

    void GridBlockedBoundaryConditionSolver3::constrainVelocity(
        FaceCenteredGrid3* velocity,
        unsigned int extrapolationDepth) {
        GridFractionalBoundaryConditionSolver3::constrainVelocity(
            velocity, extrapolationDepth);

        // No-flux: project the velocity at the marker interface
        Size3 size = velocity->resolution();
        auto u = velocity->uAccessor();
        auto v = velocity->vAccessor();
        auto w = velocity->wAccessor();
        auto uPos = velocity->uPosition();
        auto vPos = velocity->vPosition();
        auto wPos = velocity->wPosition();

        _marker.forEachIndex([&](size_t i, size_t j, size_t k) {
            if (_marker(i, j, k) == kCollider) {
                if (i > 0 && _marker(i - 1, j, k) == blocked_boundary_kfluid) {
                    Vector3D colliderVel = collider()->velocityAt(uPos(i, j, k));
                    u(i, j, k) = colliderVel.x;
                }
                if (i < size.x - 1 && _marker(i + 1, j, k) == blocked_boundary_kfluid) {
                    Vector3D colliderVel
                        = collider()->velocityAt(uPos(i + 1, j, k));
                    u(i + 1, j, k) = colliderVel.x;
                }
                if (j > 0 && _marker(i, j - 1, k) == blocked_boundary_kfluid) {
                    Vector3D colliderVel = collider()->velocityAt(vPos(i, j, k));
                    v(i, j, k) = colliderVel.y;
                }
                if (j < size.y - 1 && _marker(i, j + 1, k) == blocked_boundary_kfluid) {
                    Vector3D colliderVel
                        = collider()->velocityAt(vPos(i, j + 1, k));
                    v(i, j + 1, k) = colliderVel.y;
                }
                if (k > 0 && _marker(i, j, k - 1) == blocked_boundary_kfluid) {
                    Vector3D colliderVel = collider()->velocityAt(wPos(i, j, k));
                    w(i, j, k) = colliderVel.z;
                }
                if (k < size.z - 1 && _marker(i, j, k + 1) == blocked_boundary_kfluid) {
                    Vector3D colliderVel
                        = collider()->velocityAt(wPos(i, j, k + 1));
                    w(i, j, k + 1) = colliderVel.z;
                }
            }
            });
    }

    const Array3<char>& GridBlockedBoundaryConditionSolver3::marker() const {
        return _marker;
    }

    void GridBlockedBoundaryConditionSolver3::onColliderUpdated(
        const Size3& gridSize,
        const Vector3D& gridSpacing,
        const Vector3D& gridOrigin) {
        GridFractionalBoundaryConditionSolver3::onColliderUpdated(
            gridSize, gridSpacing, gridOrigin);

        const auto sdf
            = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(colliderSdf());

        _marker.resize(gridSize);
        _marker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf((*sdf)(i, j, k))) {
                _marker(i, j, k) = kCollider;
            }
            else {
                _marker(i, j, k) = blocked_boundary_kfluid;
            }
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_BLOCKED_BOUNDARY_CONDITION_SOLVER3_H_