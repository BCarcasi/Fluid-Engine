#ifndef INCLUDE_JET_GRID_FLUID_SOLVER3_H_
#define INCLUDE_JET_GRID_FLUID_SOLVER3_H_

#include "advection_solver3.h"
#include "cell_centered_scalar_grid3.h"
#include "collider3.h"
#include "face_centered_grid3.h"
#include "grid_boundary_condition_solver3.h"
#include "grid_diffusion_solver3.h"
#include "grid_emitter3.h"
#include "grid_pressure_solver3.h"
#include "grid_system_data3.h"
#include "physics_animation.h"

#include "pch.h"

#include "array_utils.h"
#include "constant_scalar_field3.h"
#include "constants.h"
#include "cubic_semi_lagrangian3.h"
#include "grid_backward_euler_diffusion_solver3.h"
#include "grid_blocked_boundary_condition_solver3.h"
#include "grid_fluid_solver3.h"
#include "grid_fractional_single_phase_pressure_solver3.h"
#include "level_set_utils.h"
#include "surface_to_implicit3.h"
#include "timer.h"

#include <algorithm>

namespace jet {

    //!
    //! \brief Abstract base class for grid-based 3-D fluid solver.
    //!
    //! This is an abstract base class for grid-based 3-D fluid solver based on
    //! Jos Stam's famous 1999 paper - "Stable Fluids". This solver takes fractional
    //! step method as its foundation which is consisted of independent advection,
    //! diffusion, external forces, and pressure projection steps. Each step is
    //! configurable so that a custom step can be implemented. For example, if a
    //! user wants to change the advection solver to her/his own implementation,
    //! simply call GridFluidSolver3::setAdvectionSolver(newSolver).
    //!
    class GridFluidSolver3 : public PhysicsAnimation {
    public:
        class Builder;

        //! Default constructor.
        GridFluidSolver3();

        //! Constructs solver with initial grid size.
        GridFluidSolver3(const Size3& resolution, const Vector3D& gridSpacing,
            const Vector3D& gridOrigin);

        //! Default destructor.
        virtual ~GridFluidSolver3();

        //! Returns the gravity vector of the system.
        const Vector3D& gravity() const;

        //! Sets the gravity of the system.
        void setGravity(const Vector3D& newGravity);

        //! Returns the viscosity coefficient.
        double viscosityCoefficient() const;

        //!
        //! \brief Sets the viscosity coefficient.
        //!
        //! This function sets the viscosity coefficient. Non-positive input will be
        //! clamped to zero.
        //!
        //! \param[in] newValue The new viscosity coefficient value.
        //!
        void setViscosityCoefficient(double newValue);

        //!
        //! \brief Returns the CFL number from the current velocity field for given
        //!     time interval.
        //!
        //! \param[in] timeIntervalInSeconds The time interval in seconds.
        //!
        double cfl(double timeIntervalInSeconds) const;

        //! Returns the max allowed CFL number.
        double maxCfl() const;

        //! Sets the max allowed CFL number.
        void setMaxCfl(double newCfl);

        //! Returns true if the solver is using compressed linear system.
        bool useCompressedLinearSystem() const;

        //! Sets whether the solver should use compressed linear system.
        void setUseCompressedLinearSystem(bool onoff);

        //! Returns the advection solver instance.
        const AdvectionSolver3Ptr& advectionSolver() const;

        //! Sets the advection solver.
        void setAdvectionSolver(const AdvectionSolver3Ptr& newSolver);

        //! Returns the diffusion solver instance.
        const GridDiffusionSolver3Ptr& diffusionSolver() const;

        //! Sets the diffusion solver.
        void setDiffusionSolver(const GridDiffusionSolver3Ptr& newSolver);

        //! Returns the pressure solver instance.
        const GridPressureSolver3Ptr& pressureSolver() const;

        //! Sets the pressure solver.
        void setPressureSolver(const GridPressureSolver3Ptr& newSolver);

        //! Returns the closed domain boundary flag.
        int closedDomainBoundaryFlag() const;

        //! Sets the closed domain boundary flag.
        void setClosedDomainBoundaryFlag(int flag);

        //!
        //! \brief Returns the grid system data.
        //!
        //! This function returns the grid system data. The grid system data stores
        //! the core fluid flow fields such as velocity. By default, the data
        //! instance has velocity field only.
        //!
        //! \see GridSystemData3
        //!
        const GridSystemData3Ptr& gridSystemData() const;

        //!
        //! \brief Resizes grid system data.
        //!
        //! This function resizes grid system data. You can also resize the grid by
        //! calling resize function directly from
        //! GridFluidSolver3::gridSystemData(), but this function provides a
        //! shortcut for the same operation.
        //!
        //! \param[in] newSize        The new size.
        //! \param[in] newGridSpacing The new grid spacing.
        //! \param[in] newGridOrigin  The new grid origin.
        //!
        void resizeGrid(const Size3& newSize, const Vector3D& newGridSpacing,
            const Vector3D& newGridOrigin);

        //!
        //! \brief Returns the resolution of the grid system data.
        //!
        //! This function returns the resolution of the grid system data. This is
        //! equivalent to calling gridSystemData()->resolution(), but provides a
        //! shortcut.
        //!
        Size3 resolution() const;

        //!
        //! \brief Returns the grid spacing of the grid system data.
        //!
        //! This function returns the resolution of the grid system data. This is
        //! equivalent to calling gridSystemData()->gridSpacing(), but provides a
        //! shortcut.
        //!
        Vector3D gridSpacing() const;

        //!
        //! \brief Returns the origin of the grid system data.
        //!
        //! This function returns the resolution of the grid system data. This is
        //! equivalent to calling gridSystemData()->origin(), but provides a
        //! shortcut.
        //!
        Vector3D gridOrigin() const;

        //!
        //! \brief Returns the velocity field.
        //!
        //! This function returns the velocity field from the grid system data.
        //! It is just a shortcut to the most commonly accessed data chunk.
        //!
        const FaceCenteredGrid3Ptr& velocity() const;

        //! Returns the collider.
        const Collider3Ptr& collider() const;

        //! Sets the collider.
        void setCollider(const Collider3Ptr& newCollider);

        //! Returns the emitter.
        const GridEmitter3Ptr& emitter() const;

        //! Sets the emitter.
        void setEmitter(const GridEmitter3Ptr& newEmitter);

        //! Returns builder fox GridFluidSolver3.
        static Builder builder();

    protected:
        //! Called when it needs to setup initial condition.
        void onInitialize() override;

        //! Called when advancing a single time-step.
        void onAdvanceTimeStep(double timeIntervalInSeconds) override;

        //!
        //! \brief Returns the required sub-time-steps for given time interval.
        //!
        //! This function returns the required sub-time-steps for given time
        //! interval based on the max allowed CFL number. If the time interval is
        //! too large so that it makes the CFL number greater than the max value,
        //! This function will return a numebr that is greater than 1.
        //!
        //! \see GridFluidSolver3::maxCfl
        //!
        unsigned int numberOfSubTimeSteps(
            double timeIntervalInSeconds) const override;

        //! Called at the beginning of a time-step.
        virtual void onBeginAdvanceTimeStep(double timeIntervalInSeconds);

        //! Called at the end of a time-step.
        virtual void onEndAdvanceTimeStep(double timeIntervalInSeconds);

        //!
        //! \brief Computes the external force terms.
        //!
        //! This function computes the external force applied for given time
        //! interval. By default, it only computes the gravity.
        //!
        //! \see GridFluidSolver3::computeGravity
        //!
        virtual void computeExternalForces(double timeIntervalInSeconds);

        //! Computes the viscosity term using the diffusion solver.
        virtual void computeViscosity(double timeIntervalInSeconds);

        //! Computes the pressure term using the pressure solver.
        virtual void computePressure(double timeIntervalInSeconds);

        //! Computes the advection term using the advection solver.
        virtual void computeAdvection(double timeIntervalInSeconds);

        //!
        //! \breif Returns the signed-distance representation of the fluid.
        //!
        //! This function returns the signed-distance representation of the fluid.
        //! Positive sign area is considered to be atmosphere and won't be included
        //! for computing the dynamics. By default, this will return constant scalar
        //! field of -kMaxD, meaning that the entire volume is occupied with fluid.
        //!
        virtual ScalarField3Ptr fluidSdf() const;

        //! Computes the gravity term.
        void computeGravity(double timeIntervalInSeconds);

        //!
        //! \brief Applies the boundary condition to the velocity field.
        //!
        //! This function applies the boundary condition to the velocity field by
        //! constraining the flow based on the boundary condition solver.
        //!
        void applyBoundaryCondition();

        //! Extrapolates given field into the collider-occupied region.
        void extrapolateIntoCollider(ScalarGrid3* grid);

        //! Extrapolates given field into the collider-occupied region.
        void extrapolateIntoCollider(CollocatedVectorGrid3* grid);

        //! Extrapolates given field into the collider-occupied region.
        void extrapolateIntoCollider(FaceCenteredGrid3* grid);

        //! Returns the signed-distance field representation of the collider.
        ScalarField3Ptr colliderSdf() const;

        //! Returns the velocity field of the collider.
        VectorField3Ptr colliderVelocityField() const;

    private:
        Vector3D _gravity = Vector3D(0.0, -9.8, 0.0);
        double _viscosityCoefficient = 0.0;
        double _maxCfl = 5.0;
        bool _useCompressedLinearSys = false;
        int _closedDomainBoundaryFlag = kDirectionAll;

        GridSystemData3Ptr _grids;
        Collider3Ptr _collider;
        GridEmitter3Ptr _emitter;

        AdvectionSolver3Ptr _advectionSolver;
        GridDiffusionSolver3Ptr _diffusionSolver;
        GridPressureSolver3Ptr _pressureSolver;
        GridBoundaryConditionSolver3Ptr _boundaryConditionSolver;

        void beginAdvanceTimeStep(double timeIntervalInSeconds);

        void endAdvanceTimeStep(double timeIntervalInSeconds);

        void updateCollider(double timeIntervalInSeconds);

        void updateEmitter(double timeIntervalInSeconds);
    };

    //! Shared pointer type for the GridFluidSolver3.
    typedef std::shared_ptr<GridFluidSolver3> GridFluidSolver3Ptr;

    //!
    //! \brief Base class for grid-based fluid solver builder.
    //!
    template <typename DerivedBuilder>
    class GridFluidSolverBuilderBase3 {
    public:
        //! Returns builder with grid resolution.
        DerivedBuilder& withResolution(const Size3& resolution);

        //! Returns builder with grid spacing.
        DerivedBuilder& withGridSpacing(const Vector3D& gridSpacing);

        //! Returns builder with grid spacing.
        DerivedBuilder& withGridSpacing(double gridSpacing);

        //!
        //! \brief Returns builder with domain size in x-direction.
        //!
        //! To build a solver, one can use either grid spacing directly or domain
        //! size in x-direction to set the final grid spacing.
        //!
        DerivedBuilder& withDomainSizeX(double domainSizeX);

        //! Returns builder with grid origin
        DerivedBuilder& withOrigin(const Vector3D& gridOrigin);

    protected:
        Size3 _resolution{ 1, 1, 1 };
        Vector3D _gridSpacing{ 1, 1, 1 };
        Vector3D _gridOrigin{ 0, 0, 0 };
        double _domainSizeX = 1.0;
        bool _useDomainSize = false;

        Vector3D getGridSpacing() const;
    };

    template <typename T>
    T& GridFluidSolverBuilderBase3<T>::withResolution(const Size3& resolution) {
        _resolution = resolution;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& GridFluidSolverBuilderBase3<T>::withGridSpacing(
        const Vector3D& gridSpacing) {
        _gridSpacing = gridSpacing;
        _useDomainSize = false;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& GridFluidSolverBuilderBase3<T>::withGridSpacing(double gridSpacing) {
        _gridSpacing.x = gridSpacing;
        _gridSpacing.y = gridSpacing;
        _gridSpacing.z = gridSpacing;
        _useDomainSize = false;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& GridFluidSolverBuilderBase3<T>::withDomainSizeX(double domainSizeX) {
        _domainSizeX = domainSizeX;
        _useDomainSize = true;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& GridFluidSolverBuilderBase3<T>::withOrigin(const Vector3D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return static_cast<T&>(*this);
    }

    template <typename T>
    Vector3D GridFluidSolverBuilderBase3<T>::getGridSpacing() const {
        Vector3D gridSpacing = _gridSpacing;
        if (_useDomainSize) {
            gridSpacing.set(_domainSizeX / static_cast<double>(_resolution.x));
        }
        return gridSpacing;
    }

    //!
    //! \brief Front-end to create GridFluidSolver3 objects step by step.
    //!
    class GridFluidSolver3::Builder final
        : public GridFluidSolverBuilderBase3<GridFluidSolver3::Builder> {
    public:
        //! Builds GridFluidSolver3.
        GridFluidSolver3 build() const;

        //! Builds shared pointer of GridFluidSolver3 instance.
        GridFluidSolver3Ptr makeShared() const {
            return std::make_shared<GridFluidSolver3>(_resolution, getGridSpacing(),
                _gridOrigin);
        }
    };

    GridFluidSolver3::GridFluidSolver3()
        : GridFluidSolver3({ 1, 1, 1 }, { 1, 1, 1 }, { 0, 0, 0 }) {}

    GridFluidSolver3::GridFluidSolver3(const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& gridOrigin) {
        _grids = std::make_shared<GridSystemData3>();
        _grids->resize(resolution, gridSpacing, gridOrigin);

        setAdvectionSolver(std::make_shared<CubicSemiLagrangian3>());
        setDiffusionSolver(std::make_shared<GridBackwardEulerDiffusionSolver3>());
        setPressureSolver(
            std::make_shared<GridFractionalSinglePhasePressureSolver3>());
        setIsUsingFixedSubTimeSteps(false);
    }

    GridFluidSolver3::~GridFluidSolver3() {}

    const Vector3D& GridFluidSolver3::gravity() const { return _gravity; }

    void GridFluidSolver3::setGravity(const Vector3D& newGravity) {
        _gravity = newGravity;
    }

    double GridFluidSolver3::viscosityCoefficient() const {
        return _viscosityCoefficient;
    }

    void GridFluidSolver3::setViscosityCoefficient(double newValue) {
        _viscosityCoefficient = std::max(newValue, 0.0);
    }

    double GridFluidSolver3::cfl(double timeIntervalInSeconds) const {
        auto vel = _grids->velocity();
        double maxVel = 0.0;
        vel->forEachCellIndex([&](size_t i, size_t j, size_t k) {
            Vector3D v =
                vel->valueAtCellCenter(i, j, k) + timeIntervalInSeconds * _gravity;
            maxVel = std::max(maxVel, v.x);
            maxVel = std::max(maxVel, v.y);
            maxVel = std::max(maxVel, v.z);
            });

        Vector3D gridSpacing = _grids->gridSpacing();
        double minGridSize = min3(gridSpacing.x, gridSpacing.y, gridSpacing.z);

        return maxVel * timeIntervalInSeconds / minGridSize;
    }

    double GridFluidSolver3::maxCfl() const { return _maxCfl; }

    void GridFluidSolver3::setMaxCfl(double newCfl) {
        _maxCfl = std::max(newCfl, kEpsilonD);
    }

    bool GridFluidSolver3::useCompressedLinearSystem() const {
        return _useCompressedLinearSys;
    }

    void GridFluidSolver3::setUseCompressedLinearSystem(bool onoff) {
        _useCompressedLinearSys = onoff;
    }

    const AdvectionSolver3Ptr& GridFluidSolver3::advectionSolver() const {
        return _advectionSolver;
    }

    void GridFluidSolver3::setAdvectionSolver(
        const AdvectionSolver3Ptr& newSolver) {
        _advectionSolver = newSolver;
    }

    const GridDiffusionSolver3Ptr& GridFluidSolver3::diffusionSolver() const {
        return _diffusionSolver;
    }

    void GridFluidSolver3::setDiffusionSolver(
        const GridDiffusionSolver3Ptr& newSolver) {
        _diffusionSolver = newSolver;
    }

    const GridPressureSolver3Ptr& GridFluidSolver3::pressureSolver() const {
        return _pressureSolver;
    }

    void GridFluidSolver3::setPressureSolver(
        const GridPressureSolver3Ptr& newSolver) {
        _pressureSolver = newSolver;
        if (_pressureSolver != nullptr) {
            _boundaryConditionSolver =
                _pressureSolver->suggestedBoundaryConditionSolver();

            // Apply domain boundary flag
            _boundaryConditionSolver->setClosedDomainBoundaryFlag(
                _closedDomainBoundaryFlag);
        }
    }

    int GridFluidSolver3::closedDomainBoundaryFlag() const {
        return _closedDomainBoundaryFlag;
    }

    void GridFluidSolver3::setClosedDomainBoundaryFlag(int flag) {
        _closedDomainBoundaryFlag = flag;
        _boundaryConditionSolver->setClosedDomainBoundaryFlag(
            _closedDomainBoundaryFlag);
    }

    const GridSystemData3Ptr& GridFluidSolver3::gridSystemData() const {
        return _grids;
    }

    void GridFluidSolver3::resizeGrid(const Size3& newSize,
        const Vector3D& newGridSpacing,
        const Vector3D& newGridOrigin) {
        _grids->resize(newSize, newGridSpacing, newGridOrigin);
    }

    Size3 GridFluidSolver3::resolution() const { return _grids->resolution(); }

    Vector3D GridFluidSolver3::gridSpacing() const { return _grids->gridSpacing(); }

    Vector3D GridFluidSolver3::gridOrigin() const { return _grids->origin(); }

    const FaceCenteredGrid3Ptr& GridFluidSolver3::velocity() const {
        return _grids->velocity();
    }

    const Collider3Ptr& GridFluidSolver3::collider() const { return _collider; }

    void GridFluidSolver3::setCollider(const Collider3Ptr& newCollider) {
        _collider = newCollider;
    }

    const GridEmitter3Ptr& GridFluidSolver3::emitter() const { return _emitter; }

    void GridFluidSolver3::setEmitter(const GridEmitter3Ptr& newEmitter) {
        _emitter = newEmitter;
    }

    void GridFluidSolver3::onInitialize() {
        // When initializing the solver, update the collider and emitter state as
        // well since they also affects the initial condition of the simulation.
        Timer timer;
        updateCollider(0.0);
        JET_INFO << "Update collider took " << timer.durationInSeconds()
            << " seconds";

        timer.reset();
        updateEmitter(0.0);
        JET_INFO << "Update emitter took " << timer.durationInSeconds()
            << " seconds";
    }

    void GridFluidSolver3::onAdvanceTimeStep(double timeIntervalInSeconds) {
        // The minimum grid resolution is 1x1.
        if (_grids->resolution().x == 0 || _grids->resolution().y == 0 ||
            _grids->resolution().z == 0) {
            JET_WARN << "Empty grid. Skipping the simulation.";
            return;
        }

        beginAdvanceTimeStep(timeIntervalInSeconds);

        Timer timer;
        computeExternalForces(timeIntervalInSeconds);
        JET_INFO << "Computing external force took " << timer.durationInSeconds()
            << " seconds";

        timer.reset();
        computeViscosity(timeIntervalInSeconds);
        JET_INFO << "Computing viscosity force took " << timer.durationInSeconds()
            << " seconds";

        timer.reset();
        computePressure(timeIntervalInSeconds);
        JET_INFO << "Computing pressure force took " << timer.durationInSeconds()
            << " seconds";

        timer.reset();
        computeAdvection(timeIntervalInSeconds);
        JET_INFO << "Computing advection force took " << timer.durationInSeconds()
            << " seconds";

        endAdvanceTimeStep(timeIntervalInSeconds);
    }

    unsigned int GridFluidSolver3::numberOfSubTimeSteps(
        double timeIntervalInSeconds) const {
        double currentCfl = cfl(timeIntervalInSeconds);
        return static_cast<unsigned int>(
            std::max(std::ceil(currentCfl / _maxCfl), 1.0));
    }

    void GridFluidSolver3::onBeginAdvanceTimeStep(double timeIntervalInSeconds) {
        UNUSED_VARIABLE(timeIntervalInSeconds);
    }

    void GridFluidSolver3::onEndAdvanceTimeStep(double timeIntervalInSeconds) {
        UNUSED_VARIABLE(timeIntervalInSeconds);
    }

    void GridFluidSolver3::computeExternalForces(double timeIntervalInSeconds) {
        computeGravity(timeIntervalInSeconds);
    }

    void GridFluidSolver3::computeViscosity(double timeIntervalInSeconds) {
        if (_diffusionSolver != nullptr && _viscosityCoefficient > kEpsilonD) {
            auto vel = velocity();
            auto vel0 = std::dynamic_pointer_cast<FaceCenteredGrid3>(vel->clone());

            _diffusionSolver->solve(*vel0, _viscosityCoefficient,
                timeIntervalInSeconds, vel.get(),
                *colliderSdf(), *fluidSdf());
            applyBoundaryCondition();
        }
    }

    void GridFluidSolver3::computePressure(double timeIntervalInSeconds) {
        if (_pressureSolver != nullptr) {
            auto vel = velocity();
            auto vel0 = std::dynamic_pointer_cast<FaceCenteredGrid3>(vel->clone());

            _pressureSolver->solve(*vel0, timeIntervalInSeconds, vel.get(),
                *colliderSdf(), *colliderVelocityField(),
                *fluidSdf(), _useCompressedLinearSys);
            applyBoundaryCondition();
        }
    }

    void GridFluidSolver3::computeAdvection(double timeIntervalInSeconds) {
        auto vel = velocity();
        if (_advectionSolver != nullptr) {
            // Solve advections for custom scalar fields
            size_t n = _grids->numberOfAdvectableScalarData();
            for (size_t i = 0; i < n; ++i) {
                auto grid = _grids->advectableScalarDataAt(i);
                auto grid0 = grid->clone();
                _advectionSolver->advect(*grid0, *vel, timeIntervalInSeconds,
                    grid.get(), *colliderSdf());
                extrapolateIntoCollider(grid.get());
            }

            // Solve advections for custom vector fields
            n = _grids->numberOfAdvectableVectorData();
            size_t velIdx = _grids->velocityIndex();
            for (size_t i = 0; i < n; ++i) {
                // Handle velocity layer separately.
                if (i == velIdx) {
                    continue;
                }

                auto grid = _grids->advectableVectorDataAt(i);
                auto grid0 = grid->clone();

                auto collocated =
                    std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid);
                auto collocated0 =
                    std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid0);
                if (collocated != nullptr) {
                    _advectionSolver->advect(*collocated0, *vel,
                        timeIntervalInSeconds,
                        collocated.get(), *colliderSdf());
                    extrapolateIntoCollider(collocated.get());
                    continue;
                }

                auto faceCentered =
                    std::dynamic_pointer_cast<FaceCenteredGrid3>(grid);
                auto faceCentered0 =
                    std::dynamic_pointer_cast<FaceCenteredGrid3>(grid0);
                if (faceCentered != nullptr && faceCentered0 != nullptr) {
                    _advectionSolver->advect(*faceCentered0, *vel,
                        timeIntervalInSeconds,
                        faceCentered.get(), *colliderSdf());
                    extrapolateIntoCollider(faceCentered.get());
                    continue;
                }
            }

            // Solve velocity advection
            auto vel0 = std::dynamic_pointer_cast<FaceCenteredGrid3>(vel->clone());
            _advectionSolver->advect(*vel0, *vel0, timeIntervalInSeconds, vel.get(),
                *colliderSdf());
            applyBoundaryCondition();
        }
    }

    ScalarField3Ptr GridFluidSolver3::fluidSdf() const {
        return std::make_shared<ConstantScalarField3>(-kMaxD);
    }

    void GridFluidSolver3::computeGravity(double timeIntervalInSeconds) {
        if (_gravity.lengthSquared() > kEpsilonD) {
            auto vel = _grids->velocity();
            auto u = vel->uAccessor();
            auto v = vel->vAccessor();
            auto w = vel->wAccessor();

            if (std::abs(_gravity.x) > kEpsilonD) {
                vel->forEachUIndex([&](size_t i, size_t j, size_t k) {
                    u(i, j, k) += timeIntervalInSeconds * _gravity.x;
                    });
            }

            if (std::abs(_gravity.y) > kEpsilonD) {
                vel->forEachVIndex([&](size_t i, size_t j, size_t k) {
                    v(i, j, k) += timeIntervalInSeconds * _gravity.y;
                    });
            }

            if (std::abs(_gravity.z) > kEpsilonD) {
                vel->forEachWIndex([&](size_t i, size_t j, size_t k) {
                    w(i, j, k) += timeIntervalInSeconds * _gravity.z;
                    });
            }

            applyBoundaryCondition();
        }
    }

    void GridFluidSolver3::applyBoundaryCondition() {
        auto vel = _grids->velocity();

        if (vel != nullptr && _boundaryConditionSolver != nullptr) {
            unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
            _boundaryConditionSolver->constrainVelocity(vel.get(), depth);
        }
    }

    void GridFluidSolver3::extrapolateIntoCollider(ScalarGrid3* grid) {
        Array3<char> marker(grid->dataSize());
        auto pos = grid->dataPosition();
        marker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(colliderSdf()->sample(pos(i, j, k)))) {
                marker(i, j, k) = 0;
            }
            else {
                marker(i, j, k) = 1;
            }
            });

        unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
        extrapolateToRegion(grid->constDataAccessor(), marker, depth,
            grid->dataAccessor());
    }

    void GridFluidSolver3::extrapolateIntoCollider(CollocatedVectorGrid3* grid) {
        Array3<char> marker(grid->dataSize());
        auto pos = grid->dataPosition();
        marker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(colliderSdf()->sample(pos(i, j, k)))) {
                marker(i, j, k) = 0;
            }
            else {
                marker(i, j, k) = 1;
            }
            });

        unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
        extrapolateToRegion(grid->constDataAccessor(), marker, depth,
            grid->dataAccessor());
    }

    void GridFluidSolver3::extrapolateIntoCollider(FaceCenteredGrid3* grid) {
        auto u = grid->uAccessor();
        auto v = grid->vAccessor();
        auto w = grid->wAccessor();
        auto uPos = grid->uPosition();
        auto vPos = grid->vPosition();
        auto wPos = grid->wPosition();

        Array3<char> uMarker(u.size());
        Array3<char> vMarker(v.size());
        Array3<char> wMarker(w.size());

        uMarker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(colliderSdf()->sample(uPos(i, j, k)))) {
                uMarker(i, j, k) = 0;
            }
            else {
                uMarker(i, j, k) = 1;
            }
            });

        vMarker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(colliderSdf()->sample(vPos(i, j, k)))) {
                vMarker(i, j, k) = 0;
            }
            else {
                vMarker(i, j, k) = 1;
            }
            });

        wMarker.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (isInsideSdf(colliderSdf()->sample(wPos(i, j, k)))) {
                wMarker(i, j, k) = 0;
            }
            else {
                wMarker(i, j, k) = 1;
            }
            });

        unsigned int depth = static_cast<unsigned int>(std::ceil(_maxCfl));
        extrapolateToRegion(grid->uConstAccessor(), uMarker, depth, u);
        extrapolateToRegion(grid->vConstAccessor(), vMarker, depth, v);
        extrapolateToRegion(grid->wConstAccessor(), wMarker, depth, w);
    }

    ScalarField3Ptr GridFluidSolver3::colliderSdf() const {
        return _boundaryConditionSolver->colliderSdf();
    }

    VectorField3Ptr GridFluidSolver3::colliderVelocityField() const {
        return _boundaryConditionSolver->colliderVelocityField();
    }

    void GridFluidSolver3::beginAdvanceTimeStep(double timeIntervalInSeconds) {
        // Update collider and emitter
        Timer timer;
        updateCollider(timeIntervalInSeconds);
        JET_INFO << "Update collider took " << timer.durationInSeconds()
            << " seconds";

        timer.reset();
        updateEmitter(timeIntervalInSeconds);
        JET_INFO << "Update emitter took " << timer.durationInSeconds()
            << " seconds";

        // Update boundary condition solver
        if (_boundaryConditionSolver != nullptr) {
            _boundaryConditionSolver->updateCollider(
                _collider, _grids->resolution(), _grids->gridSpacing(),
                _grids->origin());
        }

        // Apply boundary condition to the velocity field in case the field got
        // updated externally.
        applyBoundaryCondition();

        // Invoke callback
        onBeginAdvanceTimeStep(timeIntervalInSeconds);
    }

    void GridFluidSolver3::endAdvanceTimeStep(double timeIntervalInSeconds) {
        // Invoke callback
        onEndAdvanceTimeStep(timeIntervalInSeconds);
    }

    void GridFluidSolver3::updateCollider(double timeIntervalInSeconds) {
        if (_collider != nullptr) {
            _collider->update(currentTimeInSeconds(), timeIntervalInSeconds);
        }
    }

    void GridFluidSolver3::updateEmitter(double timeIntervalInSeconds) {
        if (_emitter != nullptr) {
            _emitter->update(currentTimeInSeconds(), timeIntervalInSeconds);
        }
    }

    GridFluidSolver3::Builder GridFluidSolver3::builder() { return Builder(); }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_FLUID_SOLVER3_H_