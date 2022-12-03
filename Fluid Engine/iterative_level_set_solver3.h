#ifndef INCLUDE_JET_ITERATIVE_LEVEL_SET_SOLVER3_H_
#define INCLUDE_JET_ITERATIVE_LEVEL_SET_SOLVER3_H_

#include "level_set_solver3.h"

#include "pch.h"
#include "array_utils.h"
#include "fdm_utils.h"
#include "iterative_level_set_solver3.h"
#include "parallel.h"

#include <algorithm>
#include <limits>
#include <utility>

namespace jet {

    //!
    //! \brief Abstract base class for 3-D PDE-based iterative level set solver.
    //!
    //! This class provides infrastructure for 3-D PDE-based iterative level set
    //! solver. Internally, the class implements upwind-style wave propagation and
    //! the inheriting classes must provide a way to compute the derivatives for
    //! given grid points.
    //!
    //! \see Osher, Stanley, and Ronald Fedkiw. Level set methods and dynamic
    //!     implicit surfaces. Vol. 153. Springer Science & Business Media, 2006.
    //!
    class IterativeLevelSetSolver3 : public LevelSetSolver3 {
    public:
        //! Default constructor.
        IterativeLevelSetSolver3();

        //! Default destructor.
        virtual ~IterativeLevelSetSolver3();

        //!
        //! Reinitializes given scalar field to signed-distance field.
        //!
        //! \param inputSdf Input signed-distance field which can be distorted.
        //! \param maxDistance Max range of reinitialization.
        //! \param outputSdf Output signed-distance field.
        //!
        void reinitialize(const ScalarGrid3& inputSdf, double maxDistance,
            ScalarGrid3* outputSdf) override;

        //!
        //! Extrapolates given scalar field from negative to positive SDF region.
        //!
        //! \param input Input scalar field to be extrapolated.
        //! \param sdf Reference signed-distance field.
        //! \param maxDistance Max range of extrapolation.
        //! \param output Output scalar field.
        //!
        void extrapolate(const ScalarGrid3& input, const ScalarField3& sdf,
            double maxDistance, ScalarGrid3* output) override;

        //!
        //! Extrapolates given collocated vector field from negative to positive SDF
        //! region.
        //!
        //! \param input Input collocated vector field to be extrapolated.
        //! \param sdf Reference signed-distance field.
        //! \param maxDistance Max range of extrapolation.
        //! \param output Output collocated vector field.
        //!
        void extrapolate(const CollocatedVectorGrid3& input,
            const ScalarField3& sdf, double maxDistance,
            CollocatedVectorGrid3* output) override;

        //!
        //! Extrapolates given face-centered vector field from negative to positive
        //! SDF region.
        //!
        //! \param input Input face-centered field to be extrapolated.
        //! \param sdf Reference signed-distance field.
        //! \param maxDistance Max range of extrapolation.
        //! \param output Output face-centered vector field.
        //!
        void extrapolate(const FaceCenteredGrid3& input, const ScalarField3& sdf,
            double maxDistance, FaceCenteredGrid3* output) override;

        //! Returns the maximum CFL limit.
        double maxCfl() const;

        //!
        //! \brief Sets the maximum CFL limit.
        //!
        //! This function sets the maximum CFL limit for the internal upwind-style
        //! PDE calculation. The negative input will be clamped to 0.
        //!
        void setMaxCfl(double newMaxCfl);

    protected:
        //! Computes the derivatives for given grid point.
        virtual void getDerivatives(ConstArrayAccessor3<double> grid,
            const Vector3D& gridSpacing, size_t i, size_t j,
            size_t k, std::array<double, 2>* dx,
            std::array<double, 2>* dy,
            std::array<double, 2>* dz) const = 0;

    private:
        double _maxCfl = 0.5;

        void extrapolate(const ConstArrayAccessor3<double>& input,
            const ConstArrayAccessor3<double>& sdf,
            const Vector3D& gridSpacing, double maxDistance,
            ArrayAccessor3<double> output);

        static unsigned int distanceToNumberOfIterations(double distance,
            double dtau);

        static double sign(const ConstArrayAccessor3<double>& sdf,
            const Vector3D& gridSpacing, size_t i, size_t j,
            size_t k);

        double pseudoTimeStep(ConstArrayAccessor3<double> sdf,
            const Vector3D& gridSpacing);
    };

    typedef std::shared_ptr<IterativeLevelSetSolver3> IterativeLevelSetSolver3Ptr;

    IterativeLevelSetSolver3::IterativeLevelSetSolver3() {
    }

    IterativeLevelSetSolver3::~IterativeLevelSetSolver3() {
    }

    void IterativeLevelSetSolver3::reinitialize(
        const ScalarGrid3& inputSdf,
        double maxDistance,
        ScalarGrid3* outputSdf) {
        const Size3 size = inputSdf.dataSize();
        const Vector3D gridSpacing = inputSdf.gridSpacing();

        JET_THROW_INVALID_ARG_IF(!inputSdf.hasSameShape(*outputSdf));

        ArrayAccessor3<double> outputAcc = outputSdf->dataAccessor();

        const double dtau = pseudoTimeStep(
            inputSdf.constDataAccessor(), gridSpacing);
        const unsigned int numberOfIterations
            = distanceToNumberOfIterations(maxDistance, dtau);

        copyRange3(
            inputSdf.constDataAccessor(), size.x, size.y, size.z, &outputAcc);

        Array3<double> temp(size);
        ArrayAccessor3<double> tempAcc = temp.accessor();

        JET_INFO << "Reinitializing with pseudoTimeStep: " << dtau
            << " numberOfIterations: " << numberOfIterations;

        for (unsigned int n = 0; n < numberOfIterations; ++n) {
            inputSdf.parallelForEachDataPointIndex(
                [&](size_t i, size_t j, size_t k) {
                    double s = sign(outputAcc, gridSpacing, i, j, k);

                    std::array<double, 2> dx, dy, dz;

                    getDerivatives(outputAcc, gridSpacing, i, j, k, &dx, &dy, &dz);

                    // Explicit Euler step
                    double val = outputAcc(i, j, k)
                        - dtau * std::max(s, 0.0)
                        * (std::sqrt(square(std::max(dx[0], 0.0))
                            + square(std::min(dx[1], 0.0))
                            + square(std::max(dy[0], 0.0))
                            + square(std::min(dy[1], 0.0))
                            + square(std::max(dz[0], 0.0))
                            + square(std::min(dz[1], 0.0))) - 1.0)
                        - dtau * std::min(s, 0.0)
                        * (std::sqrt(square(std::min(dx[0], 0.0))
                            + square(std::max(dx[1], 0.0))
                            + square(std::min(dy[0], 0.0))
                            + square(std::max(dy[1], 0.0))
                            + square(std::min(dz[0], 0.0))
                            + square(std::max(dz[1], 0.0))) - 1.0);
                    tempAcc(i, j, k) = val;
                });

            std::swap(tempAcc, outputAcc);
        }

        auto outputSdfAcc = outputSdf->dataAccessor();
        copyRange3(outputAcc, size.x, size.y, size.z, &outputSdfAcc);
    }

    void IterativeLevelSetSolver3::extrapolate(
        const ScalarGrid3& input,
        const ScalarField3& sdf,
        double maxDistance,
        ScalarGrid3* output) {
        JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));

        Array3<double> sdfGrid(input.dataSize());
        auto pos = input.dataPosition();
        sdfGrid.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            sdfGrid(i, j, k) = sdf.sample(pos(i, j, k));
            });

        extrapolate(
            input.constDataAccessor(),
            sdfGrid.constAccessor(),
            input.gridSpacing(),
            maxDistance,
            output->dataAccessor());
    }

    void IterativeLevelSetSolver3::extrapolate(
        const CollocatedVectorGrid3& input,
        const ScalarField3& sdf,
        double maxDistance,
        CollocatedVectorGrid3* output) {
        JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));

        Array3<double> sdfGrid(input.dataSize());
        auto pos = input.dataPosition();
        sdfGrid.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            sdfGrid(i, j, k) = sdf.sample(pos(i, j, k));
            });

        const Vector3D gridSpacing = input.gridSpacing();

        Array3<double> u(input.dataSize());
        Array3<double> u0(input.dataSize());
        Array3<double> v(input.dataSize());
        Array3<double> v0(input.dataSize());
        Array3<double> w(input.dataSize());
        Array3<double> w0(input.dataSize());

        input.parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            u(i, j, k) = input(i, j, k).x;
            v(i, j, k) = input(i, j, k).y;
            w(i, j, k) = input(i, j, k).z;
            });

        extrapolate(
            u,
            sdfGrid.constAccessor(),
            gridSpacing,
            maxDistance,
            u0);

        extrapolate(
            v,
            sdfGrid.constAccessor(),
            gridSpacing,
            maxDistance,
            v0);

        extrapolate(
            w,
            sdfGrid.constAccessor(),
            gridSpacing,
            maxDistance,
            w0);

        output->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            (*output)(i, j, k).x = u(i, j, k);
            (*output)(i, j, k).y = v(i, j, k);
            (*output)(i, j, k).z = w(i, j, k);
            });
    }

    void IterativeLevelSetSolver3::extrapolate(
        const FaceCenteredGrid3& input,
        const ScalarField3& sdf,
        double maxDistance,
        FaceCenteredGrid3* output) {
        JET_THROW_INVALID_ARG_IF(!input.hasSameShape(*output));

        const Vector3D gridSpacing = input.gridSpacing();

        auto u = input.uConstAccessor();
        auto uPos = input.uPosition();
        Array3<double> sdfAtU(u.size());
        input.parallelForEachUIndex([&](size_t i, size_t j, size_t k) {
            sdfAtU(i, j, k) = sdf.sample(uPos(i, j, k));
            });

        extrapolate(
            u,
            sdfAtU,
            gridSpacing,
            maxDistance,
            output->uAccessor());

        auto v = input.vConstAccessor();
        auto vPos = input.vPosition();
        Array3<double> sdfAtV(v.size());
        input.parallelForEachVIndex([&](size_t i, size_t j, size_t k) {
            sdfAtV(i, j, k) = sdf.sample(vPos(i, j, k));
            });

        extrapolate(
            v,
            sdfAtV,
            gridSpacing,
            maxDistance,
            output->vAccessor());

        auto w = input.wConstAccessor();
        auto wPos = input.wPosition();
        Array3<double> sdfAtW(w.size());
        input.parallelForEachWIndex([&](size_t i, size_t j, size_t k) {
            sdfAtW(i, j, k) = sdf.sample(wPos(i, j, k));
            });

        extrapolate(
            w,
            sdfAtW,
            gridSpacing,
            maxDistance,
            output->wAccessor());
    }

    void IterativeLevelSetSolver3::extrapolate(
        const ConstArrayAccessor3<double>& input,
        const ConstArrayAccessor3<double>& sdf,
        const Vector3D& gridSpacing,
        double maxDistance,
        ArrayAccessor3<double> output) {
        const Size3 size = input.size();

        ArrayAccessor3<double> outputAcc = output;

        const double dtau = pseudoTimeStep(sdf, gridSpacing);
        const unsigned int numberOfIterations
            = distanceToNumberOfIterations(maxDistance, dtau);

        copyRange3(input, size.x, size.y, size.z, &outputAcc);

        Array3<double> temp(size);
        ArrayAccessor3<double> tempAcc = temp.accessor();

        for (unsigned int n = 0; n < numberOfIterations; ++n) {
            parallelFor(
                kZeroSize, size.x, kZeroSize, size.y, kZeroSize, size.z,
                [&](size_t i, size_t j, size_t k) {
                    if (sdf(i, j, k) >= 0) {
                        std::array<double, 2> dx, dy, dz;
                        Vector3D grad = gradient3(sdf, gridSpacing, i, j, k);

                        getDerivatives(
                            outputAcc, gridSpacing, i, j, k, &dx, &dy, &dz);

                        tempAcc(i, j, k) = outputAcc(i, j, k)
                            - dtau * (std::max(grad.x, 0.0) * dx[0]
                                + std::min(grad.x, 0.0) * dx[1]
                                + std::max(grad.y, 0.0) * dy[0]
                                + std::min(grad.y, 0.0) * dy[1]
                                + std::max(grad.z, 0.0) * dz[0]
                                + std::min(grad.z, 0.0) * dz[1]);
                    }
                    else {
                        tempAcc(i, j, k) = outputAcc(i, j, k);
                    }
                });

            std::swap(tempAcc, outputAcc);
        }

        copyRange3(outputAcc, size.x, size.y, size.z, &output);
    }

    double IterativeLevelSetSolver3::maxCfl() const {
        return _maxCfl;
    }

    void IterativeLevelSetSolver3::setMaxCfl(double newMaxCfl) {
        _maxCfl = std::max(newMaxCfl, 0.0);
    }

    unsigned int IterativeLevelSetSolver3::distanceToNumberOfIterations(
        double distance,
        double dtau) {
        return static_cast<unsigned int>(std::ceil(distance / dtau));
    }

    double IterativeLevelSetSolver3::sign(
        const ConstArrayAccessor3<double>& sdf,
        const Vector3D& gridSpacing,
        size_t i,
        size_t j,
        size_t k) {
        double d = sdf(i, j, k);
        double e = min3(gridSpacing.x, gridSpacing.y, gridSpacing.z);
        return d / std::sqrt(d * d + e * e);
    }

    double IterativeLevelSetSolver3::pseudoTimeStep(
        ConstArrayAccessor3<double> sdf,
        const Vector3D& gridSpacing) {
        const Size3 size = sdf.size();

        const double h = max3(gridSpacing.x, gridSpacing.y, gridSpacing.z);

        double maxS = -std::numeric_limits<double>::max();
        double dtau = _maxCfl * h;

        for (size_t k = 0; k < size.z; ++k) {
            for (size_t j = 0; j < size.y; ++j) {
                for (size_t i = 0; i < size.x; ++i) {
                    double s = sign(sdf, gridSpacing, i, j, k);
                    maxS = std::max(s, maxS);
                }
            }
        }

        while (dtau * maxS / h > _maxCfl) {
            dtau *= 0.5;
        }

        return dtau;
    }

}  // namespace jet

#endif  // INCLUDE_JET_ITERATIVE_LEVEL_SET_SOLVER3_H_