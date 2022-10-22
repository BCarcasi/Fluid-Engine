#ifndef INCLUDE_JET_PARTICLE_SYSTEM_SOLVER3_H_
#define INCLUDE_JET_PARTICLE_SYSTEM_SOLVER3_H_

#include "collider3.h"
#include "constants.h"
#include "vector_field3.h"
#include "particle_emitter3.h"
#include "particle_system_data3.h"
#include "physics_animation.h"

#include "pch.h"

#include "array_utils.h"
#include "constant_vector_field3.h"
#include "parallel.h"
#include "timer.h"

#include <algorithm>


namespace jet {

    //!
    //! \brief      Basic 3-D particle system solver.
    //!
    //! This class implements basic particle system solver. It includes gravity,
    //! air drag, and collision. But it does not compute particle-to-particle
    //! interaction. Thus, this solver is suitable for performing simple spray-like
    //! simulations with low computational cost. This class can be further extend
    //! to add more sophisticated simulations, such as SPH, to handle
    //! particle-to-particle intersection.
    //!
    //! \see        SphSolver3
    //!
    class ParticleSystemSolver3 : public PhysicsAnimation {
    public:
        class Builder;

        //! Constructs an empty solver.
        ParticleSystemSolver3();

        //! Constructs a solver with particle parameters.
        ParticleSystemSolver3(
            double radius,
            double mass);

        //! Destructor.
        virtual ~ParticleSystemSolver3();

        //! Returns the drag coefficient.
        double dragCoefficient() const;

        //!
        //! \brief      Sets the drag coefficient.
        //!
        //! The drag coefficient controls the amount of air-drag. The coefficient
        //! should be a positive number and 0 means no drag force.
        //!
        //! \param[in]  newDragCoefficient The new drag coefficient.
        //!
        void setDragCoefficient(double newDragCoefficient);

        //! Sets the restitution coefficient.
        double restitutionCoefficient() const;

        //!
        //! \brief      Sets the restitution coefficient.
        //!
        //! The restitution coefficient controls the bouncy-ness of a particle when
        //! it hits a collider surface. The range of the coefficient should be 0 to
        //! 1 -- 0 means no bounce back and 1 means perfect reflection.
        //!
        //! \param[in]  newRestitutionCoefficient The new restitution coefficient.
        //!
        void setRestitutionCoefficient(double newRestitutionCoefficient);

        //! Returns the gravity.
        const Vector3D& gravity() const;

        //! Sets the gravity.
        void setGravity(const Vector3D& newGravity);

        //!
        //! \brief      Returns the particle system data.
        //!
        //! This function returns the particle system data. The data is created when
        //! this solver is constructed and also owned by the solver.
        //!
        //! \return     The particle system data.
        //!
        const ParticleSystemData3Ptr& particleSystemData() const;

        //! Returns the collider.
        const Collider3Ptr& collider() const;

        //! Sets the collider.
        void setCollider(const Collider3Ptr& newCollider);

        //! Returns the emitter.
        const ParticleEmitter3Ptr& emitter() const;

        //! Sets the emitter.
        void setEmitter(const ParticleEmitter3Ptr& newEmitter);

        //! Returns the wind field.
        const VectorField3Ptr& wind() const;

        //!
        //! \brief      Sets the wind.
        //!
        //! Wind can be applied to the particle system by setting a vector field to
        //! the solver.
        //!
        //! \param[in]  newWind The new wind.
        //!
        void setWind(const VectorField3Ptr& newWind);

        //! Returns builder fox ParticleSystemSolver3.
        static Builder builder();

    protected:
        //! Initializes the simulator.
        void onInitialize() override;

        //! Called to advane a single time-step.
        void onAdvanceTimeStep(double timeStepInSeconds) override;

        //! Accumulates forces applied to the particles.
        virtual void accumulateForces(double timeStepInSeconds);

        //! Called when a time-step is about to begin.
        virtual void onBeginAdvanceTimeStep(double timeStepInSeconds);

        //! Called after a time-step is completed.
        virtual void onEndAdvanceTimeStep(double timeStepInSeconds);

        //! Resolves any collisions occured by the particles.
        void resolveCollision();

        //! Resolves any collisions occured by the particles where the particle
        //! state is given by the position and velocity arrays.
        void resolveCollision(
            ArrayAccessor1<Vector3D> newPositions,
            ArrayAccessor1<Vector3D> newVelocities);

        //! Assign a new particle system data.
        void setParticleSystemData(const ParticleSystemData3Ptr& newParticles);

    private:
        double _dragCoefficient = 1e-4;
        double _restitutionCoefficient = 0.0;
        Vector3D _gravity = Vector3D(0.0, kGravity, 0.0);

        ParticleSystemData3Ptr _particleSystemData;
        ParticleSystemData3::VectorData _newPositions;
        ParticleSystemData3::VectorData _newVelocities;
        Collider3Ptr _collider;
        ParticleEmitter3Ptr _emitter;
        VectorField3Ptr _wind;

        void beginAdvanceTimeStep(double timeStepInSeconds);

        void endAdvanceTimeStep(double timeStepInSeconds);

        void accumulateExternalForces();

        void timeIntegration(double timeStepInSeconds);

        void updateCollider(double timeStepInSeconds);

        void updateEmitter(double timeStepInSeconds);
    };

    //! Shared pointer type for the ParticleSystemSolver3.
    typedef std::shared_ptr<ParticleSystemSolver3> ParticleSystemSolver3Ptr;


    //!
    //! \brief Base class for particle-based solver builder.
    //!
    template <typename DerivedBuilder>
    class ParticleSystemSolverBuilderBase3 {
    public:
        //! Returns builder with particle radius.
        DerivedBuilder& withRadius(double radius);

        //! Returns builder with mass per particle.
        DerivedBuilder& withMass(double mass);

    protected:
        double _radius = 1e-3;
        double _mass = 1e-3;
    };

    template <typename T>
    T& ParticleSystemSolverBuilderBase3<T>::withRadius(double radius) {
        _radius = radius;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& ParticleSystemSolverBuilderBase3<T>::withMass(double mass) {
        _mass = mass;
        return static_cast<T&>(*this);
    }

    //!
    //! \brief Front-end to create ParticleSystemSolver3 objects step by step.
    //!
    class ParticleSystemSolver3::Builder final
        : public ParticleSystemSolverBuilderBase3<ParticleSystemSolver3::Builder> {
    public:
        //! Builds ParticleSystemSolver3.
        ParticleSystemSolver3 build() const;

        //! Builds shared pointer of ParticleSystemSolver3 instance.
        ParticleSystemSolver3Ptr makeShared() const;
    };

    ParticleSystemSolver3::ParticleSystemSolver3()
        : ParticleSystemSolver3(1e-3, 1e-3) {
    }

    ParticleSystemSolver3::ParticleSystemSolver3(
        double radius,
        double mass) {
        _particleSystemData = std::make_shared<ParticleSystemData3>();
        _particleSystemData->setRadius(radius);
        _particleSystemData->setMass(mass);
        _wind = std::make_shared<ConstantVectorField3>(Vector3D());
    }

    ParticleSystemSolver3::~ParticleSystemSolver3() {
    }

    double ParticleSystemSolver3::dragCoefficient() const {
        return _dragCoefficient;
    }

    void ParticleSystemSolver3::setDragCoefficient(double newDragCoefficient) {
        _dragCoefficient = std::max(newDragCoefficient, 0.0);
    }

    double ParticleSystemSolver3::restitutionCoefficient() const {
        return _restitutionCoefficient;
    }

    void ParticleSystemSolver3::setRestitutionCoefficient(
        double newRestitutionCoefficient) {
        _restitutionCoefficient = clamp(newRestitutionCoefficient, 0.0, 1.0);
    }

    const Vector3D& ParticleSystemSolver3::gravity() const {
        return _gravity;
    }

    void ParticleSystemSolver3::setGravity(const Vector3D& newGravity) {
        _gravity = newGravity;
    }

    const ParticleSystemData3Ptr&
        ParticleSystemSolver3::particleSystemData() const {
        return _particleSystemData;
    }

    const Collider3Ptr& ParticleSystemSolver3::collider() const {
        return _collider;
    }

    void ParticleSystemSolver3::setCollider(
        const Collider3Ptr& newCollider) {
        _collider = newCollider;
    }

    const ParticleEmitter3Ptr& ParticleSystemSolver3::emitter() const {
        return _emitter;
    }

    void ParticleSystemSolver3::setEmitter(
        const ParticleEmitter3Ptr& newEmitter) {
        _emitter = newEmitter;
        newEmitter->setTarget(_particleSystemData);
    }

    const VectorField3Ptr& ParticleSystemSolver3::wind() const {
        return _wind;
    }

    void ParticleSystemSolver3::setWind(const VectorField3Ptr& newWind) {
        _wind = newWind;
    }

    void ParticleSystemSolver3::onInitialize() {
        // When initializing the solver, update the collider and emitter state as
        // well since they also affects the initial condition of the simulation.
        Timer timer;
        updateCollider(0.0);
        JET_INFO << "Update collider took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        updateEmitter(0.0);
        JET_INFO << "Update emitter took "
            << timer.durationInSeconds() << " seconds";
    }

    void ParticleSystemSolver3::onAdvanceTimeStep(double timeStepInSeconds) {
        beginAdvanceTimeStep(timeStepInSeconds);

        Timer timer;
        accumulateForces(timeStepInSeconds);
        JET_INFO << "Accumulating forces took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        timeIntegration(timeStepInSeconds);
        JET_INFO << "Time integration took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        resolveCollision();
        JET_INFO << "Resolving collision took "
            << timer.durationInSeconds() << " seconds";

        endAdvanceTimeStep(timeStepInSeconds);
    }

    void ParticleSystemSolver3::accumulateForces(double timeStepInSeconds) {
        UNUSED_VARIABLE(timeStepInSeconds);

        // Add external forces
        accumulateExternalForces();
    }

    void ParticleSystemSolver3::beginAdvanceTimeStep(double timeStepInSeconds) {
        // Clear forces
        auto forces = _particleSystemData->forces();
        setRange1(forces.size(), Vector3D(), &forces);

        // Update collider and emitter
        Timer timer;
        updateCollider(timeStepInSeconds);
        JET_INFO << "Update collider took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        updateEmitter(timeStepInSeconds);
        JET_INFO << "Update emitter took "
            << timer.durationInSeconds() << " seconds";

        // Allocate buffers
        size_t n = _particleSystemData->numberOfParticles();
        _newPositions.resize(n);
        _newVelocities.resize(n);

        onBeginAdvanceTimeStep(timeStepInSeconds);
    }

    void ParticleSystemSolver3::endAdvanceTimeStep(double timeStepInSeconds) {
        // Update data
        size_t n = _particleSystemData->numberOfParticles();
        auto positions = _particleSystemData->positions();
        auto velocities = _particleSystemData->velocities();
        parallelFor(
            kZeroSize,
            n,
            [&](size_t i) {
                positions[i] = _newPositions[i];
                velocities[i] = _newVelocities[i];
            });

        onEndAdvanceTimeStep(timeStepInSeconds);
    }

    void ParticleSystemSolver3::onBeginAdvanceTimeStep(double timeStepInSeconds) {
        UNUSED_VARIABLE(timeStepInSeconds);
    }

    void ParticleSystemSolver3::onEndAdvanceTimeStep(double timeStepInSeconds) {
        UNUSED_VARIABLE(timeStepInSeconds);
    }

    void ParticleSystemSolver3::resolveCollision() {
        resolveCollision(
            _newPositions.accessor(),
            _newVelocities.accessor());
    }

    void ParticleSystemSolver3::resolveCollision(
        ArrayAccessor1<Vector3D> newPositions,
        ArrayAccessor1<Vector3D> newVelocities) {
        if (_collider != nullptr) {
            size_t numberOfParticles = _particleSystemData->numberOfParticles();
            const double radius = _particleSystemData->radius();

            parallelFor(
                kZeroSize,
                numberOfParticles,
                [&](size_t i) {
                    _collider->resolveCollision(
                        radius,
                        _restitutionCoefficient,
                        &newPositions[i],
                        &newVelocities[i]);
                });
        }
    }

    void ParticleSystemSolver3::setParticleSystemData(
        const ParticleSystemData3Ptr& newParticles) {
        _particleSystemData = newParticles;
    }

    void ParticleSystemSolver3::accumulateExternalForces() {
        size_t n = _particleSystemData->numberOfParticles();
        auto forces = _particleSystemData->forces();
        auto velocities = _particleSystemData->velocities();
        auto positions = _particleSystemData->positions();
        const double mass = _particleSystemData->mass();

        parallelFor(
            kZeroSize,
            n,
            [&](size_t i) {
                // Gravity
                Vector3D force = mass * _gravity;

                // Wind forces
                Vector3D relativeVel = velocities[i] - _wind->sample(positions[i]);
                force += -_dragCoefficient * relativeVel;

                forces[i] += force;
            });
    }

    void ParticleSystemSolver3::timeIntegration(double timeStepInSeconds) {
        size_t n = _particleSystemData->numberOfParticles();
        auto forces = _particleSystemData->forces();
        auto velocities = _particleSystemData->velocities();
        auto positions = _particleSystemData->positions();
        const double mass = _particleSystemData->mass();

        parallelFor(
            kZeroSize,
            n,
            [&](size_t i) {
                // Integrate velocity first
                Vector3D& newVelocity = _newVelocities[i];
                newVelocity = velocities[i]
                    + timeStepInSeconds * forces[i] / mass;

                // Integrate position.
                Vector3D& newPosition = _newPositions[i];
                newPosition = positions[i] + timeStepInSeconds * newVelocity;
            });
    }

    void ParticleSystemSolver3::updateCollider(double timeStepInSeconds) {
        if (_collider != nullptr) {
            _collider->update(currentTimeInSeconds(), timeStepInSeconds);
        }
    }

    void ParticleSystemSolver3::updateEmitter(double timeStepInSeconds) {
        if (_emitter != nullptr) {
            _emitter->update(currentTimeInSeconds(), timeStepInSeconds);
        }
    }

    ParticleSystemSolver3::Builder ParticleSystemSolver3::builder() {
        return Builder();
    }


    ParticleSystemSolver3 ParticleSystemSolver3::Builder::build() const {
        return ParticleSystemSolver3(_radius, _mass);
    }

    ParticleSystemSolver3Ptr ParticleSystemSolver3::Builder::makeShared() const {
        return std::shared_ptr<ParticleSystemSolver3>(
            new ParticleSystemSolver3(_radius, _mass),
            [](ParticleSystemSolver3* obj) {
                delete obj;
            });
    }


}  // namespace jet

#endif  // INCLUDE_JET_PARTICLE_SYSTEM_SOLVER3_H_