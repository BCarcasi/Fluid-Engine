#ifndef INCLUDE_JET_GRID_SMOKE_SOLVER3_H_
#define INCLUDE_JET_GRID_SMOKE_SOLVER3_H_

#include "pch.h"

#include "grid_fluid_solver3.h"

#include <algorithm>
namespace jet {

    //!
    //! \brief      3-D grid-based smoke solver.
    //!
    //! This class extends GridFluidSolver3 to implement smoke simulation solver.
    //! It adds smoke density and temperature fields to define the smoke and uses
    //! buoyancy force to simulate hot rising smoke.
    //!
    //! \see Fedkiw, Ronald, Jos Stam, and Henrik Wann Jensen.
    //!     "Visual simulation of smoke." Proceedings of the 28th annual conference
    //!     on Computer graphics and interactive techniques. ACM, 2001.
    //!
    class GridSmokeSolver3 : public GridFluidSolver3 {
    public:
        class Builder;

        //! Default constructor.
        GridSmokeSolver3();

        //! Constructs solver with initial grid size.
        GridSmokeSolver3(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin);

        //! Destructor.
        virtual ~GridSmokeSolver3();

        //! Returns smoke diffusion coefficient.
        double smokeDiffusionCoefficient() const;

        //! Sets smoke diffusion coefficient.
        void setSmokeDiffusionCoefficient(double newValue);

        //! Returns temperature diffusion coefficient.
        double temperatureDiffusionCoefficient() const;

        //! Sets temperature diffusion coefficient.
        void setTemperatureDiffusionCoefficient(double newValue);

        //!
        //! \brief      Returns the buoyancy factor which will be multiplied to the
        //!     smoke density.
        //!
        //! This class computes buoyancy by looking up the value of smoke density
        //! and temperature, compare them to the average values, and apply
        //! multiplier factor to the diff between the value and the average. That
        //! multiplier is defined for each smoke density and temperature separately.
        //! For example, negative smoke density buoyancy factor means a heavier
        //! smoke should sink.
        //!
        //! \return     The buoyance factor for the smoke density.
        //!
        double buoyancySmokeDensityFactor() const;

        //!
        //! \brief          Sets the buoyancy factor which will be multiplied to the
        //!     smoke density.
        //!
        //! This class computes buoyancy by looking up the value of smoke density
        //! and temperature, compare them to the average values, and apply
        //! multiplier factor to the diff between the value and the average. That
        //! multiplier is defined for each smoke density and temperature separately.
        //! For example, negative smoke density buoyancy factor means a heavier
        //! smoke should sink.
        //!
        //! \param newValue The new buoyancy factor for smoke density.
        //!
        void setBuoyancySmokeDensityFactor(double newValue);

        //!
        //! \brief      Returns the buoyancy factor which will be multiplied to the
        //!     temperature.
        //!
        //! This class computes buoyancy by looking up the value of smoke density
        //! and temperature, compare them to the average values, and apply
        //! multiplier factor to the diff between the value and the average. That
        //! multiplier is defined for each smoke density and temperature separately.
        //! For example, negative smoke density buoyancy factor means a heavier
        //! smoke should sink.
        //!
        //! \return     The buoyance factor for the temperature.
        //!
        double buoyancyTemperatureFactor() const;

        //!
        //! \brief          Sets the buoyancy factor which will be multiplied to the
        //!     temperature.
        //!
        //! This class computes buoyancy by looking up the value of smoke density
        //! and temperature, compare them to the average values, and apply
        //! multiplier factor to the diff between the value and the average. That
        //! multiplier is defined for each smoke density and temperature separately.
        //! For example, negative smoke density buoyancy factor means a heavier
        //! smoke should sink.
        //!
        //! \param newValue The new buoyancy factor for temperature.
        //!
        void setBuoyancyTemperatureFactor(double newValue);

        //!
        //! \brief      Returns smoke decay factor.
        //!
        //! In addition to the diffusion, the smoke also can fade-out over time by
        //! setting the decay factor between 0 and 1.
        //!
        //! \return     The decay factor for smoke density.
        //!
        double smokeDecayFactor() const;

        //!
        //! \brief      Sets the smoke decay factor.
        //!
        //! In addition to the diffusion, the smoke also can fade-out over time by
        //! setting the decay factor between 0 and 1.
        //!
        //! \param[in]  newValue The new decay factor.
        //!
        void setSmokeDecayFactor(double newValue);

        //!
        //! \brief      Returns temperature decay factor.
        //!
        //! In addition to the diffusion, the smoke also can fade-out over time by
        //! setting the decay factor between 0 and 1.
        //!
        //! \return     The decay factor for smoke temperature.
        //!
        double smokeTemperatureDecayFactor() const;

        //!
        //! \brief      Sets the temperature decay factor.
        //!
        //! In addition to the diffusion, the temperature also can fade-out over
        //! time by setting the decay factor between 0 and 1.
        //!
        //! \param[in]  newValue The new decay factor.
        //!
        void setTemperatureDecayFactor(double newValue);

        //! Returns smoke density field.
        ScalarGrid3Ptr smokeDensity() const;

        //! Returns temperature field.
        ScalarGrid3Ptr temperature() const;

        //! Returns builder fox GridSmokeSolver3.
        static Builder builder();

    protected:
        void onEndAdvanceTimeStep(double timeIntervalInSeconds) override;

        void computeExternalForces(double timeIntervalInSeconds) override;

    private:
        size_t _smokeDensityDataId;
        size_t _temperatureDataId;
        double _smokeDiffusionCoefficient = 0.0;
        double _temperatureDiffusionCoefficient = 0.0;
        double _buoyancySmokeDensityFactor = -0.000625;
        double _buoyancyTemperatureFactor = 5.0;
        double _smokeDecayFactor = 0.001;
        double _temperatureDecayFactor = 0.001;

        void computeDiffusion(double timeIntervalInSeconds);

        void computeBuoyancyForce(double timeIntervalInSeconds);
    };

    //! Shared pointer type for the GridSmokeSolver3.
    typedef std::shared_ptr<GridSmokeSolver3> GridSmokeSolver3Ptr;


    //!
    //! \brief Front-end to create GridSmokeSolver3 objects step by step.
    //!
    class GridSmokeSolver3::Builder final
        : public GridFluidSolverBuilderBase3<GridSmokeSolver3::Builder> {
    public:
        //! Builds GridSmokeSolver3.
        GridSmokeSolver3 build() const;

        //! Builds shared pointer of GridSmokeSolver3 instance.
        GridSmokeSolver3Ptr makeShared() const;
    };

    GridSmokeSolver3::GridSmokeSolver3()
        : GridSmokeSolver3({ 1, 1, 1 }, { 1, 1, 1 }, { 0, 0, 0 }) {}

    GridSmokeSolver3::GridSmokeSolver3(const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& gridOrigin)
        : GridFluidSolver3(resolution, gridSpacing, gridOrigin) {
        auto grids = gridSystemData();
        _smokeDensityDataId = grids->addAdvectableScalarData(
            std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);
        _temperatureDataId = grids->addAdvectableScalarData(
            std::make_shared<CellCenteredScalarGrid3::Builder>(), 0.0);
    }

    GridSmokeSolver3::~GridSmokeSolver3() {}

    double GridSmokeSolver3::smokeDiffusionCoefficient() const {
        return _smokeDiffusionCoefficient;
    }

    void GridSmokeSolver3::setSmokeDiffusionCoefficient(double newValue) {
        _smokeDiffusionCoefficient = std::max(newValue, 0.0);
    }

    double GridSmokeSolver3::temperatureDiffusionCoefficient() const {
        return _temperatureDiffusionCoefficient;
    }

    void GridSmokeSolver3::setTemperatureDiffusionCoefficient(double newValue) {
        _temperatureDiffusionCoefficient = std::max(newValue, 0.0);
    }

    double GridSmokeSolver3::buoyancySmokeDensityFactor() const {
        return _buoyancySmokeDensityFactor;
    }

    void GridSmokeSolver3::setBuoyancySmokeDensityFactor(double newValue) {
        _buoyancySmokeDensityFactor = newValue;
    }

    double GridSmokeSolver3::buoyancyTemperatureFactor() const {
        return _buoyancyTemperatureFactor;
    }

    void GridSmokeSolver3::setBuoyancyTemperatureFactor(double newValue) {
        _buoyancyTemperatureFactor = newValue;
    }

    double GridSmokeSolver3::smokeDecayFactor() const { return _smokeDecayFactor; }

    void GridSmokeSolver3::setSmokeDecayFactor(double newValue) {
        _smokeDecayFactor = clamp(newValue, 0.0, 1.0);
    }

    double GridSmokeSolver3::smokeTemperatureDecayFactor() const {
        return _temperatureDecayFactor;
    }

    void GridSmokeSolver3::setTemperatureDecayFactor(double newValue) {
        _temperatureDecayFactor = clamp(newValue, 0.0, 1.0);
    }

    ScalarGrid3Ptr GridSmokeSolver3::smokeDensity() const {
        return gridSystemData()->advectableScalarDataAt(_smokeDensityDataId);
    }

    ScalarGrid3Ptr GridSmokeSolver3::temperature() const {
        return gridSystemData()->advectableScalarDataAt(_temperatureDataId);
    }

    void GridSmokeSolver3::onEndAdvanceTimeStep(double timeIntervalInSeconds) {
        computeDiffusion(timeIntervalInSeconds);
    }

    void GridSmokeSolver3::computeExternalForces(double timeIntervalInSeconds) {
        computeBuoyancyForce(timeIntervalInSeconds);
    }

    void GridSmokeSolver3::computeDiffusion(double timeIntervalInSeconds) {
        if (diffusionSolver() != nullptr) {
            if (_smokeDiffusionCoefficient > kEpsilonD) {
                auto den = smokeDensity();
                auto den0 = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(
                    den->clone());

                diffusionSolver()->solve(*den0, _smokeDiffusionCoefficient,
                    timeIntervalInSeconds, den.get(),
                    *colliderSdf());
                extrapolateIntoCollider(den.get());
            }

            if (_temperatureDiffusionCoefficient > kEpsilonD) {
                auto temp = temperature();
                auto temp0 = std::dynamic_pointer_cast<CellCenteredScalarGrid3>(
                    temp->clone());

                diffusionSolver()->solve(*temp0, _temperatureDiffusionCoefficient,
                    timeIntervalInSeconds, temp.get(),
                    *colliderSdf());
                extrapolateIntoCollider(temp.get());
            }
        }

        auto den = smokeDensity();
        den->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            (*den)(i, j, k) *= 1.0 - _smokeDecayFactor;
            });
        auto temp = temperature();
        temp->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            (*temp)(i, j, k) *= 1.0 - _temperatureDecayFactor;
            });
    }

    void GridSmokeSolver3::computeBuoyancyForce(double timeIntervalInSeconds) {
        auto grids = gridSystemData();
        auto vel = grids->velocity();

        Vector3D up(0, 1, 0);
        if (gravity().lengthSquared() > kEpsilonD) {
            up = -gravity().normalized();
        }

        if (std::abs(_buoyancySmokeDensityFactor) > kEpsilonD ||
            std::abs(_buoyancyTemperatureFactor) > kEpsilonD) {
            auto den = smokeDensity();
            auto temp = temperature();

            double tAmb = 0.0;
            temp->forEachCellIndex(
                [&](size_t i, size_t j, size_t k) { tAmb += (*temp)(i, j, k); });
            tAmb /= static_cast<double>(
                temp->resolution().x * temp->resolution().y * temp->resolution().z);

            auto u = vel->uAccessor();
            auto v = vel->vAccessor();
            auto w = vel->wAccessor();
            auto uPos = vel->uPosition();
            auto vPos = vel->vPosition();
            auto wPos = vel->wPosition();

            if (std::abs(up.x) > kEpsilonD) {
                vel->parallelForEachUIndex([&](size_t i, size_t j, size_t k) {
                    Vector3D pt = uPos(i, j, k);
                    double fBuoy =
                        _buoyancySmokeDensityFactor * den->sample(pt) +
                        _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                    u(i, j, k) += timeIntervalInSeconds * fBuoy * up.x;
                    });
            }

            if (std::abs(up.y) > kEpsilonD) {
                vel->parallelForEachVIndex([&](size_t i, size_t j, size_t k) {
                    Vector3D pt = vPos(i, j, k);
                    double fBuoy =
                        _buoyancySmokeDensityFactor * den->sample(pt) +
                        _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                    v(i, j, k) += timeIntervalInSeconds * fBuoy * up.y;
                    });
            }

            if (std::abs(up.z) > kEpsilonD) {
                vel->parallelForEachWIndex([&](size_t i, size_t j, size_t k) {
                    Vector3D pt = wPos(i, j, k);
                    double fBuoy =
                        _buoyancySmokeDensityFactor * den->sample(pt) +
                        _buoyancyTemperatureFactor * (temp->sample(pt) - tAmb);
                    w(i, j, k) += timeIntervalInSeconds * fBuoy * up.z;
                    });
            }

            applyBoundaryCondition();
        }
    }

    GridSmokeSolver3::Builder GridSmokeSolver3::builder() { return Builder(); }

    GridSmokeSolver3 GridSmokeSolver3::Builder::build() const {
        return GridSmokeSolver3(_resolution, getGridSpacing(), _gridOrigin);
    }

    GridSmokeSolver3Ptr GridSmokeSolver3::Builder::makeShared() const {
        return std::shared_ptr<GridSmokeSolver3>(
            new GridSmokeSolver3(_resolution, getGridSpacing(), _gridOrigin),
            [](GridSmokeSolver3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_SMOKE_SOLVER3_H_