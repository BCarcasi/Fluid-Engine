#ifndef INCLUDE_JET_GRID_BACKWARD_EULER_DIFFUSION_SOLVER3_H_
#define INCLUDE_JET_GRID_BACKWARD_EULER_DIFFUSION_SOLVER3_H_

#include "constant_scalar_field3.h"
#include "fdm_linear_system_solver3.h"
#include "grid_diffusion_solver3.h"
#include <limits>
#include <memory>

#include "pch.h"
#include "constants.h"
#include "fdm_iccg_solver3.h"
#include "fdm_utils.h"
#include "level_set_utils.h"

namespace jet {

    //!
    //! \brief 3-D grid-based backward Euler diffusion solver.
    //!
    //! This class implements 3-D grid-based forward Euler diffusion solver using
    //! second-order central differencing spatially. Since the method is following
    //! the implicit time-integration (i.e. backward Euler), larger time interval or
    //! diffusion coefficient can be used without breaking the result. Note, higher
    //! values for those parameters will still impact the accuracy of the result.
    //! To solve the backward Euler method, a linear system solver is used and
    //! incomplete Cholesky conjugate gradient method is used by default.
    //!
    class GridBackwardEulerDiffusionSolver3 final : public GridDiffusionSolver3 {
    public:
        enum BoundaryType {
            Dirichlet,
            Neumann
        };

        //! Constructs the solver with given boundary type.
        explicit GridBackwardEulerDiffusionSolver3(
            BoundaryType boundaryType = Neumann);

        //!
        //! Solves diffusion equation for a scalar field.
        //!
        //! \param source Input scalar field.
        //! \param diffusionCoefficient Amount of diffusion.
        //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
        //! \param dest Output scalar field.
        //! \param boundarySdf Shape of the solid boundary that is empty by default.
        //! \param boundarySdf Shape of the fluid boundary that is full by default.
        //!
        void solve(
            const ScalarGrid3& source,
            double diffusionCoefficient,
            double timeIntervalInSeconds,
            ScalarGrid3* dest,
            const ScalarField3& boundarySdf
            = ConstantScalarField3(kMaxD),
            const ScalarField3& fluidSdf
            = ConstantScalarField3(-kMaxD)) override;

        //!
        //! Solves diffusion equation for a collocated vector field.
        //!
        //! \param source Input collocated vector field.
        //! \param diffusionCoefficient Amount of diffusion.
        //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
        //! \param dest Output collocated vector field.
        //! \param boundarySdf Shape of the solid boundary that is empty by default.
        //! \param boundarySdf Shape of the fluid boundary that is full by default.
        //!
        void solve(
            const CollocatedVectorGrid3& source,
            double diffusionCoefficient,
            double timeIntervalInSeconds,
            CollocatedVectorGrid3* dest,
            const ScalarField3& boundarySdf
            = ConstantScalarField3(kMaxD),
            const ScalarField3& fluidSdf
            = ConstantScalarField3(-kMaxD)) override;

        //!
        //! Solves diffusion equation for a face-centered vector field.
        //!
        //! \param source Input face-centered vector field.
        //! \param diffusionCoefficient Amount of diffusion.
        //! \param timeIntervalInSeconds Small time-interval that diffusion occur.
        //! \param dest Output face-centered vector field.
        //! \param boundarySdf Shape of the solid boundary that is empty by default.
        //! \param boundarySdf Shape of the fluid boundary that is full by default.
        //!
        void solve(
            const FaceCenteredGrid3& source,
            double diffusionCoefficient,
            double timeIntervalInSeconds,
            FaceCenteredGrid3* dest,
            const ScalarField3& boundarySdf
            = ConstantScalarField3(kMaxD),
            const ScalarField3& fluidSdf
            = ConstantScalarField3(-kMaxD)) override;

        //! Sets the linear system solver for this diffusion solver.
        void setLinearSystemSolver(const FdmLinearSystemSolver3Ptr& solver);

    private:
        BoundaryType _boundaryType;
        FdmLinearSystem3 _system;
        FdmLinearSystemSolver3Ptr _systemSolver;
        Array3<char> _markers;

        void buildMarkers(
            const Size3& size,
            const std::function<Vector3D(size_t, size_t, size_t)>& pos,
            const ScalarField3& boundarySdf,
            const ScalarField3& fluidSdf);

        void buildMatrix(
            const Size3& size,
            const Vector3D& c);

        void buildVectors(
            const ConstArrayAccessor3<double>& f,
            const Vector3D& c);

        void buildVectors(
            const ConstArrayAccessor3<Vector3D>& f,
            const Vector3D& c,
            size_t component);
    };

    //! Shared pointer type for the GridBackwardEulerDiffusionSolver3.
    typedef std::shared_ptr<GridBackwardEulerDiffusionSolver3>
        GridBackwardEulerDiffusionSolver3Ptr;

    const char backward_euler_kfluid = 0;
    const char kAir = 1;
    const char kBoundary = 2;

    GridBackwardEulerDiffusionSolver3::GridBackwardEulerDiffusionSolver3(
        BoundaryType boundaryType) : _boundaryType(boundaryType) {
        _systemSolver = std::make_shared<FdmIccgSolver3>(100, kEpsilonD);
    }

    void GridBackwardEulerDiffusionSolver3::solve(
        const ScalarGrid3& source,
        double diffusionCoefficient,
        double timeIntervalInSeconds,
        ScalarGrid3* dest,
        const ScalarField3& boundarySdf,
        const ScalarField3& fluidSdf) {
        auto pos = source.dataPosition();
        Vector3D h = source.gridSpacing();
        Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);

        buildMarkers(source.dataSize(), pos, boundarySdf, fluidSdf);
        buildMatrix(source.dataSize(), c);
        buildVectors(source.constDataAccessor(), c);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachDataPointIndex(
                [&](size_t i, size_t j, size_t k) {
                    (*dest)(i, j, k) = _system.x(i, j, k);
                });
        }
    }

    void GridBackwardEulerDiffusionSolver3::solve(
        const CollocatedVectorGrid3& source,
        double diffusionCoefficient,
        double timeIntervalInSeconds,
        CollocatedVectorGrid3* dest,
        const ScalarField3& boundarySdf,
        const ScalarField3& fluidSdf) {
        auto pos = source.dataPosition();
        Vector3D h = source.gridSpacing();
        Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);

        buildMarkers(source.dataSize(), pos, boundarySdf, fluidSdf);
        buildMatrix(source.dataSize(), c);

        // u
        buildVectors(source.constDataAccessor(), c, 0);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachDataPointIndex(
                [&](size_t i, size_t j, size_t k) {
                    (*dest)(i, j, k).x = _system.x(i, j, k);
                });
        }

        // v
        buildVectors(source.constDataAccessor(), c, 1);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachDataPointIndex(
                [&](size_t i, size_t j, size_t k) {
                    (*dest)(i, j, k).y = _system.x(i, j, k);
                });
        }

        // w
        buildVectors(source.constDataAccessor(), c, 2);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachDataPointIndex(
                [&](size_t i, size_t j, size_t k) {
                    (*dest)(i, j, k).z = _system.x(i, j, k);
                });
        }
    }

    void GridBackwardEulerDiffusionSolver3::solve(
        const FaceCenteredGrid3& source,
        double diffusionCoefficient,
        double timeIntervalInSeconds,
        FaceCenteredGrid3* dest,
        const ScalarField3& boundarySdf,
        const ScalarField3& fluidSdf) {
        Vector3D h = source.gridSpacing();
        Vector3D c = timeIntervalInSeconds * diffusionCoefficient / (h * h);

        // u
        auto uPos = source.uPosition();
        buildMarkers(source.uSize(), uPos, boundarySdf, fluidSdf);
        buildMatrix(source.uSize(), c);
        buildVectors(source.uConstAccessor(), c);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachUIndex(
                [&](size_t i, size_t j, size_t k) {
                    dest->u(i, j, k) = _system.x(i, j, k);
                });
        }

        // v
        auto vPos = source.vPosition();
        buildMarkers(source.vSize(), vPos, boundarySdf, fluidSdf);
        buildMatrix(source.vSize(), c);
        buildVectors(source.vConstAccessor(), c);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachVIndex(
                [&](size_t i, size_t j, size_t k) {
                    dest->v(i, j, k) = _system.x(i, j, k);
                });
        }

        // w
        auto wPos = source.wPosition();
        buildMarkers(source.wSize(), wPos, boundarySdf, fluidSdf);
        buildMatrix(source.wSize(), c);
        buildVectors(source.wConstAccessor(), c);

        if (_systemSolver != nullptr) {
            // Solve the system
            _systemSolver->solve(&_system);

            // Assign the solution
            source.parallelForEachWIndex(
                [&](size_t i, size_t j, size_t k) {
                    dest->w(i, j, k) = _system.x(i, j, k);
                });
        }
    }

    void GridBackwardEulerDiffusionSolver3::setLinearSystemSolver(
        const FdmLinearSystemSolver3Ptr& solver) {
        _systemSolver = solver;
    }

    void GridBackwardEulerDiffusionSolver3::buildMarkers(
        const Size3& size,
        const std::function<Vector3D(size_t, size_t, size_t)>& pos,
        const ScalarField3& boundarySdf,
        const ScalarField3& fluidSdf) {
        _markers.resize(size);

        _markers.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(boundarySdf.sample(pos(i, j, k)))) {
                _markers(i, j, k) = kBoundary;
            }
            else if (isInsideSdf(fluidSdf.sample(pos(i, j, k)))) {
                _markers(i, j, k) = backward_euler_kfluid;
            }
            else {
                _markers(i, j, k) = kAir;
            }
            });
    }

    void GridBackwardEulerDiffusionSolver3::buildMatrix(
        const Size3& size,
        const Vector3D& c) {
        _system.A.resize(size);

        bool isDirichlet = (_boundaryType == Dirichlet);

        // Build linear system
        _system.A.parallelForEachIndex(
            [&](size_t i, size_t j, size_t k) {
                auto& row = _system.A(i, j, k);

                // Initialize
                row.center = 1.0;
                row.right = row.up = row.front = 0.0;

                if (_markers(i, j, k) == backward_euler_kfluid) {
                    if (i + 1 < size.x) {
                        if ((isDirichlet && _markers(i + 1, j, k) != kAir)
                            || _markers(i + 1, j, k) == backward_euler_kfluid) {
                            row.center += c.x;
                        }

                        if (_markers(i + 1, j, k) == backward_euler_kfluid) {
                            row.right -= c.x;
                        }
                    }

                    if (i > 0
                        && ((isDirichlet && _markers(i - 1, j, k) != kAir)
                            || _markers(i - 1, j, k) == backward_euler_kfluid)) {
                        row.center += c.x;
                    }

                    if (j + 1 < size.y) {
                        if ((isDirichlet && _markers(i, j + 1, k) != kAir)
                            || _markers(i, j + 1, k) == backward_euler_kfluid) {
                            row.center += c.y;
                        }

                        if (_markers(i, j + 1, k) == backward_euler_kfluid) {
                            row.up -= c.y;
                        }
                    }

                    if (j > 0
                        && ((isDirichlet && _markers(i, j - 1, k) != kAir)
                            || _markers(i, j - 1, k) == backward_euler_kfluid)) {
                        row.center += c.y;
                    }

                    if (k + 1 < size.z) {
                        if ((isDirichlet && _markers(i, j, k + 1) != kAir)
                            || _markers(i, j, k + 1) == backward_euler_kfluid) {
                            row.center += c.z;
                        }

                        if (_markers(i, j, k + 1) == backward_euler_kfluid) {
                            row.front -= c.z;
                        }
                    }

                    if (k > 0
                        && ((isDirichlet && _markers(i, j, k - 1) != kAir)
                            || _markers(i, j, k - 1) == backward_euler_kfluid)) {
                        row.center += c.z;
                    }
                }
            });
    }

    void GridBackwardEulerDiffusionSolver3::buildVectors(
        const ConstArrayAccessor3<double>& f,
        const Vector3D& c) {
        Size3 size = f.size();

        _system.x.resize(size, 0.0);
        _system.b.resize(size, 0.0);

        // Build linear system
        _system.x.parallelForEachIndex(
            [&](size_t i, size_t j, size_t k) {
                _system.b(i, j, k) = _system.x(i, j, k) = f(i, j, k);

                if (_boundaryType == Dirichlet && _markers(i, j, k) == backward_euler_kfluid) {
                    if (i + 1 < size.x && _markers(i + 1, j, k) == kBoundary) {
                        _system.b(i, j, k) += c.x * f(i + 1, j, k);
                    }

                    if (i > 0 && _markers(i - 1, j, k) == kBoundary) {
                        _system.b(i, j, k) += c.x * f(i - 1, j, k);
                    }

                    if (j + 1 < size.y && _markers(i, j + 1, k) == kBoundary) {
                        _system.b(i, j, k) += c.y * f(i, j + 1, k);
                    }

                    if (j > 0 && _markers(i, j - 1, k) == kBoundary) {
                        _system.b(i, j, k) += c.y * f(i, j - 1, k);
                    }

                    if (k + 1 < size.z && _markers(i, j, k + 1) == kBoundary) {
                        _system.b(i, j, k) += c.z * f(i, j, k + 1);
                    }

                    if (k > 0 && _markers(i, j, k - 1) == kBoundary) {
                        _system.b(i, j, k) += c.z * f(i, j, k - 1);
                    }
                }
            });
    }

    void GridBackwardEulerDiffusionSolver3::buildVectors(
        const ConstArrayAccessor3<Vector3D>& f,
        const Vector3D& c,
        size_t component) {
        Size3 size = f.size();

        _system.x.resize(size, 0.0);
        _system.b.resize(size, 0.0);

        // Build linear system
        _system.x.parallelForEachIndex(
            [&](size_t i, size_t j, size_t k) {
                _system.b(i, j, k) = _system.x(i, j, k) = f(i, j, k)[component];

                if (_boundaryType == Dirichlet && _markers(i, j, k) == backward_euler_kfluid) {
                    if (i + 1 < size.x && _markers(i + 1, j, k) == kBoundary) {
                        _system.b(i, j, k) += c.x * f(i + 1, j, k)[component];
                    }

                    if (i > 0 && _markers(i - 1, j, k) == kBoundary) {
                        _system.b(i, j, k) += c.x * f(i - 1, j, k)[component];
                    }

                    if (j + 1 < size.y && _markers(i, j + 1, k) == kBoundary) {
                        _system.b(i, j, k) += c.y * f(i, j + 1, k)[component];
                    }

                    if (j > 0 && _markers(i, j - 1, k) == kBoundary) {
                        _system.b(i, j, k) += c.y * f(i, j - 1, k)[component];
                    }

                    if (k + 1 < size.z && _markers(i, j, k + 1) == kBoundary) {
                        _system.b(i, j, k) += c.z * f(i, j, k + 1)[component];
                    }

                    if (k > 0 && _markers(i, j, k - 1) == kBoundary) {
                        _system.b(i, j, k) += c.z * f(i, j, k - 1)[component];
                    }
                }
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_BACKWARD_EULER_DIFFUSION_SOLVER3_H_