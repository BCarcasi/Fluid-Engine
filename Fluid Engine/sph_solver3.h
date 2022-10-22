#ifndef INCLUDE_JET_SPH_SOLVER3_H_
#define INCLUDE_JET_SPH_SOLVER3_H_

#include "constants.h"
#include "particle_system_solver3.h"
#include "sph_system_data3.h"

#include "pch.h"
#include "physics_helpers.h"
#include "parallel.h"
#include "sph_kernels3.h"
#include "timer.h"

#include <algorithm>

namespace jet {

    //!
    //! \brief 3-D SPH solver.
    //!
    //! This class implements 3-D SPH solver. The main pressure solver is based on
    //! equation-of-state (EOS).
    //!
    //! \see M{\"u}ller et al., Particle-based fluid simulation for interactive
    //!      applications, SCA 2003.
    //! \see M. Becker and M. Teschner, Weakly compressible SPH for free surface
    //!      flows, SCA 2007.
    //! \see Adams and Wicke, Meshless approximation methods and applications in
    //!      physics based modeling and animation, Eurographics tutorials 2009.
    //!
    class SphSolver3 : public ParticleSystemSolver3 {
    public:
        class Builder;

        //! Constructs a solver with empty particle set.
        SphSolver3();

        //! Constructs a solver with target density, spacing, and relative kernel
        //! radius.
        SphSolver3(
            double targetDensity,
            double targetSpacing,
            double relativeKernelRadius);

        virtual ~SphSolver3();

        //! Returns the exponent part of the equation-of-state.
        double eosExponent() const;

        //!
        //! \brief Sets the exponent part of the equation-of-state.
        //!
        //! This function sets the exponent part of the equation-of-state.
        //! The value must be greater than 1.0, and smaller inputs will be clamped.
        //! Default is 7.
        //!
        void setEosExponent(double newEosExponent);

        //! Returns the negative pressure scale.
        double negativePressureScale() const;

        //!
        //! \brief Sets the negative pressure scale.
        //!
        //! This function sets the negative pressure scale. By setting the number
        //! between 0 and 1, the solver will scale the effect of negative pressure
        //! which can prevent the clumping of the particles near the surface. Input
        //! value outside 0 and 1 will be clamped within the range. Default is 0.
        //!
        void setNegativePressureScale(double newNegativePressureScale);

        //! Returns the viscosity coefficient.
        double viscosityCoefficient() const;

        //! Sets the viscosity coefficient.
        void setViscosityCoefficient(double newViscosityCoefficient);

        //! Returns the pseudo viscosity coefficient.
        double pseudoViscosityCoefficient() const;

        //!
        //! \brief Sets the pseudo viscosity coefficient.
        //!
        //! This function sets the pseudo viscosity coefficient which applies
        //! additional pseudo-physical damping to the system. Default is 10.
        //!
        void setPseudoViscosityCoefficient(double newPseudoViscosityCoefficient);

        //! Returns the speed of sound.
        double speedOfSound() const;

        //!
        //! \brief Sets the speed of sound.
        //!
        //! This function sets the speed of sound which affects the stiffness of the
        //! EOS and the time-step size. Higher value will make EOS stiffer and the
        //! time-step smaller. The input value must be higher than 0.0.
        //!
        void setSpeedOfSound(double newSpeedOfSound);

        //!
        //! \brief Multiplier that scales the max allowed time-step.
        //!
        //! This function returns the multiplier that scales the max allowed
        //! time-step. When the scale is 1.0, the time-step is bounded by the speed
        //! of sound and max acceleration.
        //!
        double timeStepLimitScale() const;

        //!
        //! \brief Sets the multiplier that scales the max allowed time-step.
        //!
        //! This function sets the multiplier that scales the max allowed
        //! time-step. When the scale is 1.0, the time-step is bounded by the speed
        //! of sound and max acceleration.
        //!
        void setTimeStepLimitScale(double newScale);

        //! Returns the SPH system data.
        SphSystemData3Ptr sphSystemData() const;

        //! Returns builder fox SphSolver3.
        static Builder builder();

    protected:
        //! Returns the number of sub-time-steps.
        unsigned int numberOfSubTimeSteps(
            double timeIntervalInSeconds) const override;

        //! Accumulates the force to the forces array in the particle system.
        void accumulateForces(double timeStepInSeconds) override;

        //! Performs pre-processing step before the simulation.
        void onBeginAdvanceTimeStep(double timeStepInSeconds) override;

        //! Performs post-processing step before the simulation.
        void onEndAdvanceTimeStep(double timeStepInSeconds) override;

        //! Accumulates the non-pressure forces to the forces array in the particle
        //! system.
        virtual void accumulateNonPressureForces(double timeStepInSeconds);

        //! Accumulates the pressure force to the forces array in the particle
        //! system.
        virtual void accumulatePressureForce(double timeStepInSeconds);

        //! Computes the pressure.
        void computePressure();

        //! Accumulates the pressure force to the given \p pressureForces array.
        void accumulatePressureForce(
            const ConstArrayAccessor1<Vector3D>& positions,
            const ConstArrayAccessor1<double>& densities,
            const ConstArrayAccessor1<double>& pressures,
            ArrayAccessor1<Vector3D> pressureForces);

        //! Accumulates the viscosity force to the forces array in the particle
        //! system.
        void accumulateViscosityForce();

        //! Computes pseudo viscosity.
        void computePseudoViscosity(double timeStepInSeconds);

    private:
        //! Exponent component of equation-of-state (or Tait's equation).
        double _eosExponent = 7.0;

        //! Negative pressure scaling factor.
        //! Zero means clamping. One means do nothing.
        double _negativePressureScale = 0.0;

        //! Viscosity coefficient.
        double _viscosityCoefficient = 0.01;

        //! Pseudo-viscosity coefficient velocity filtering.
        //! This is a minimum "safety-net" for SPH solver which is quite
        //! sensitive to the parameters.
        double _pseudoViscosityCoefficient = 10.0;

        //! Speed of sound in medium to determin the stiffness of the system.
        //! Ideally, it should be the actual speed of sound in the fluid, but in
        //! practice, use lower value to trace-off performance and compressibility.
        double _speedOfSound = 100.0;

        //! Scales the max allowed time-step.
        double _timeStepLimitScale = 1.0;
    };

    //! Shared pointer type for the SphSolver3.
    typedef std::shared_ptr<SphSolver3> SphSolver3Ptr;


    //!
    //! \brief Base class for SPH-based fluid solver builder.
    //!
    template <typename DerivedBuilder>
    class SphSolverBuilderBase3 {
    public:
        //! Returns builder with target density.
        DerivedBuilder& withTargetDensity(double targetDensity);

        //! Returns builder with target spacing.
        DerivedBuilder& withTargetSpacing(double targetSpacing);

        //! Returns builder with relative kernel radius.
        DerivedBuilder& withRelativeKernelRadius(double relativeKernelRadius);

    protected:
        double _targetDensity = kWaterDensity;
        double _targetSpacing = 0.1;
        double _relativeKernelRadius = 1.8;
    };

    template <typename T>
    T& SphSolverBuilderBase3<T>::withTargetDensity(double targetDensity) {
        _targetDensity = targetDensity;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& SphSolverBuilderBase3<T>::withTargetSpacing(double targetSpacing) {
        _targetSpacing = targetSpacing;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& SphSolverBuilderBase3<T>::withRelativeKernelRadius(
        double relativeKernelRadius) {
        _relativeKernelRadius = relativeKernelRadius;
        return static_cast<T&>(*this);
    }

    //!
    //! \brief Front-end to create SphSolver3 objects step by step.
    //!
    class SphSolver3::Builder final
        : public SphSolverBuilderBase3<SphSolver3::Builder> {
    public:
        //! Builds SphSolver3.
        SphSolver3 build() const;

        //! Builds shared pointer of SphSolver3 instance.
        SphSolver3Ptr makeShared() const;
    };

    static double kTimeStepLimitBySpeedFactor = 0.4;
    static double kTimeStepLimitByForceFactor = 0.25;

    SphSolver3::SphSolver3() {
        setParticleSystemData(std::make_shared<SphSystemData3>());
        setIsUsingFixedSubTimeSteps(false);
    }

    SphSolver3::SphSolver3(
        double targetDensity,
        double targetSpacing,
        double relativeKernelRadius) {
        auto sphParticles = std::make_shared<SphSystemData3>();
        setParticleSystemData(sphParticles);
        sphParticles->setTargetDensity(targetDensity);
        sphParticles->setTargetSpacing(targetSpacing);
        sphParticles->setRelativeKernelRadius(relativeKernelRadius);
        setIsUsingFixedSubTimeSteps(false);
    }

    SphSolver3::~SphSolver3() {
    }

    double SphSolver3::eosExponent() const {
        return _eosExponent;
    }

    void SphSolver3::setEosExponent(double newEosExponent) {
        _eosExponent = std::max(newEosExponent, 1.0);
    }

    double SphSolver3::negativePressureScale() const {
        return _negativePressureScale;
    }

    void SphSolver3::setNegativePressureScale(
        double newNegativePressureScale) {
        _negativePressureScale = clamp(newNegativePressureScale, 0.0, 1.0);
    }

    double SphSolver3::viscosityCoefficient() const {
        return _viscosityCoefficient;
    }

    void SphSolver3::setViscosityCoefficient(double newViscosityCoefficient) {
        _viscosityCoefficient = std::max(newViscosityCoefficient, 0.0);
    }

    double SphSolver3::pseudoViscosityCoefficient() const {
        return _pseudoViscosityCoefficient;
    }

    void SphSolver3::setPseudoViscosityCoefficient(
        double newPseudoViscosityCoefficient) {
        _pseudoViscosityCoefficient
            = std::max(newPseudoViscosityCoefficient, 0.0);
    }

    double SphSolver3::speedOfSound() const {
        return _speedOfSound;
    }

    void SphSolver3::setSpeedOfSound(double newSpeedOfSound) {
        _speedOfSound = std::max(newSpeedOfSound, kEpsilonD);
    }

    double SphSolver3::timeStepLimitScale() const {
        return _timeStepLimitScale;
    }

    void SphSolver3::setTimeStepLimitScale(double newScale) {
        _timeStepLimitScale = std::max(newScale, 0.0);
    }

    SphSystemData3Ptr SphSolver3::sphSystemData() const {
        return std::dynamic_pointer_cast<SphSystemData3>(particleSystemData());
    }

    unsigned int SphSolver3::numberOfSubTimeSteps(
        double timeIntervalInSeconds) const {
        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();
        auto f = particles->forces();

        const double kernelRadius = particles->kernelRadius();
        const double mass = particles->mass();

        double maxForceMagnitude = 0.0;

        for (size_t i = 0; i < numberOfParticles; ++i) {
            maxForceMagnitude = std::max(maxForceMagnitude, f[i].length());
        }

        double timeStepLimitBySpeed
            = kTimeStepLimitBySpeedFactor * kernelRadius / _speedOfSound;
        double timeStepLimitByForce
            = kTimeStepLimitByForceFactor
            * std::sqrt(kernelRadius * mass / maxForceMagnitude);

        double desiredTimeStep
            = _timeStepLimitScale
            * std::min(timeStepLimitBySpeed, timeStepLimitByForce);

        return static_cast<unsigned int>(
            std::ceil(timeIntervalInSeconds / desiredTimeStep));
    }

    void SphSolver3::accumulateForces(double timeStepInSeconds) {
        accumulateNonPressureForces(timeStepInSeconds);
        accumulatePressureForce(timeStepInSeconds);
    }

    void SphSolver3::onBeginAdvanceTimeStep(double timeStepInSeconds) {
        UNUSED_VARIABLE(timeStepInSeconds);

        auto particles = sphSystemData();

        Timer timer;
        particles->buildNeighborSearcher();
        particles->buildNeighborLists();
        particles->updateDensities();

        JET_INFO << "Building neighbor lists and updating densities took "
            << timer.durationInSeconds()
            << " seconds";
    }

    void SphSolver3::onEndAdvanceTimeStep(double timeStepInSeconds) {
        computePseudoViscosity(timeStepInSeconds);

        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();
        auto densities = particles->densities();

        double maxDensity = 0.0;
        for (size_t i = 0; i < numberOfParticles; ++i) {
            maxDensity = std::max(maxDensity, densities[i]);
        }

        JET_INFO << "Max density: " << maxDensity << " "
            << "Max density / target density ratio: "
            << maxDensity / particles->targetDensity();
    }

    void SphSolver3::accumulateNonPressureForces(double timeStepInSeconds) {
        ParticleSystemSolver3::accumulateForces(timeStepInSeconds);
        accumulateViscosityForce();
    }

    void SphSolver3::accumulatePressureForce(double timeStepInSeconds) {
        UNUSED_VARIABLE(timeStepInSeconds);

        auto particles = sphSystemData();
        auto x = particles->positions();
        auto d = particles->densities();
        auto p = particles->pressures();
        auto f = particles->forces();

        computePressure();
        accumulatePressureForce(x, d, p, f);
    }

    void SphSolver3::computePressure() {
        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();
        auto d = particles->densities();
        auto p = particles->pressures();

        // See Murnaghan-Tait equation of state from
        // https://en.wikipedia.org/wiki/Tait_equation
        const double targetDensity = particles->targetDensity();
        const double eosScale = targetDensity * square(_speedOfSound);

        parallelFor(
            kZeroSize,
            numberOfParticles,
            [&](size_t i) {
                p[i] = computePressureFromEos(
                    d[i],
                    targetDensity,
                    eosScale,
                    eosExponent(),
                    negativePressureScale());
            });
    }

    void SphSolver3::accumulatePressureForce(
        const ConstArrayAccessor1<Vector3D>& positions,
        const ConstArrayAccessor1<double>& densities,
        const ConstArrayAccessor1<double>& pressures,
        ArrayAccessor1<Vector3D> pressureForces) {
        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();

        const double massSquared = square(particles->mass());
        const SphSpikyKernel3 kernel(particles->kernelRadius());

        parallelFor(
            kZeroSize,
            numberOfParticles,
            [&](size_t i) {
                const auto& neighbors = particles->neighborLists()[i];
                for (size_t j : neighbors) {
                    double dist = positions[i].distanceTo(positions[j]);

                    if (dist > 0.0) {
                        Vector3D dir = (positions[j] - positions[i]) / dist;
                        pressureForces[i] -= massSquared
                            * (pressures[i] / (densities[i] * densities[i])
                                + pressures[j] / (densities[j] * densities[j]))
                            * kernel.gradient(dist, dir);
                    }
                }
            });
    }


    void SphSolver3::accumulateViscosityForce() {
        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();
        auto x = particles->positions();
        auto v = particles->velocities();
        auto d = particles->densities();
        auto f = particles->forces();

        const double massSquared = square(particles->mass());
        const SphSpikyKernel3 kernel(particles->kernelRadius());

        parallelFor(
            kZeroSize,
            numberOfParticles,
            [&](size_t i) {
                const auto& neighbors = particles->neighborLists()[i];
                for (size_t j : neighbors) {
                    double dist = x[i].distanceTo(x[j]);

                    f[i] += viscosityCoefficient() * massSquared
                        * (v[j] - v[i]) / d[j]
                        * kernel.secondDerivative(dist);
                }
            });
    }

    void SphSolver3::computePseudoViscosity(double timeStepInSeconds) {
        auto particles = sphSystemData();
        size_t numberOfParticles = particles->numberOfParticles();
        auto x = particles->positions();
        auto v = particles->velocities();
        auto d = particles->densities();

        const double mass = particles->mass();
        const SphSpikyKernel3 kernel(particles->kernelRadius());

        Array1<Vector3D> smoothedVelocities(numberOfParticles);

        parallelFor(
            kZeroSize,
            numberOfParticles,
            [&](size_t i) {
                double weightSum = 0.0;
                Vector3D smoothedVelocity;

                const auto& neighbors = particles->neighborLists()[i];
                for (size_t j : neighbors) {
                    double dist = x[i].distanceTo(x[j]);
                    double wj = mass / d[j] * kernel(dist);
                    weightSum += wj;
                    smoothedVelocity += wj * v[j];
                }

                double wi = mass / d[i];
                weightSum += wi;
                smoothedVelocity += wi * v[i];

                if (weightSum > 0.0) {
                    smoothedVelocity /= weightSum;
                }

                smoothedVelocities[i] = smoothedVelocity;
            });

        double factor = timeStepInSeconds * _pseudoViscosityCoefficient;
        factor = clamp(factor, 0.0, 1.0);

        parallelFor(
            kZeroSize,
            numberOfParticles,
            [&](size_t i) {
                v[i] = lerp(
                    v[i], smoothedVelocities[i], factor);
            });
    }

    SphSolver3::Builder SphSolver3::builder() {
        return Builder();
    }

    SphSolver3 SphSolver3::Builder::build() const {
        return SphSolver3(
            _targetDensity,
            _targetSpacing,
            _relativeKernelRadius);
    }

    SphSolver3Ptr SphSolver3::Builder::makeShared() const {
        return std::shared_ptr<SphSolver3>(
            new SphSolver3(
                _targetDensity,
                _targetSpacing,
                _relativeKernelRadius),
            [](SphSolver3* obj) {
                delete obj;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_SPH_SOLVER3_H_