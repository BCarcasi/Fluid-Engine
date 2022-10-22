#ifndef INCLUDE_JET_PIC_SOLVER3_H_
#define INCLUDE_JET_PIC_SOLVER3_H_

#include "grid_fluid_solver3.h"
#include "particle_emitter3.h"
#include "particle_system_data3.h"

#include "pch.h"
#include "array_utils.h"
#include "level_set_utils.h"
#include "pic_solver3.h"
#include "timer.h"
#include <algorithm>

namespace jet {

    //!
    //! \brief 3-D Particle-in-Cell (PIC) implementation.
    //!
    //! This class implements 3-D Particle-in-Cell (PIC) method by inheriting
    //! GridFluidSolver3. Since it is a grid-particle hybrid method, the solver
    //! also has a particle system to track fluid particles.
    //!
    //! \see Zhu, Yongning, and Robert Bridson. "Animating sand as a fluid."
    //!     ACM Transactions on Graphics (TOG). Vol. 34. No. 3. ACM, 3005.
    //!
    class PicSolver3 : public GridFluidSolver3 {
    public:
        class Builder;

        //! Default constructor.
        PicSolver3();

        //! Constructs solver with initial grid size.
        PicSolver3(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin);

        //! Default destructor.
        virtual ~PicSolver3();

        //! Returns the signed-distance field of particles.
        ScalarGrid3Ptr signedDistanceField() const;

        //! Returns the particle system data.
        const ParticleSystemData3Ptr& particleSystemData() const;

        //! Returns the particle emitter.
        const ParticleEmitter3Ptr& particleEmitter() const;

        //! Sets the particle emitter.
        void setParticleEmitter(const ParticleEmitter3Ptr& newEmitter);

        //! Returns builder fox PicSolver3.
        static Builder builder();

    protected:
        Array3<char> _uMarkers;
        Array3<char> _vMarkers;
        Array3<char> _wMarkers;

        //! Initializes the simulator.
        void onInitialize() override;

        //! Invoked before a simulation time-step begins.
        void onBeginAdvanceTimeStep(double timeIntervalInSeconds) override;

        //! Computes the advection term of the fluid solver.
        void computeAdvection(double timeIntervalInSeconds) override;

        //! Returns the signed-distance field of the fluid.
        ScalarField3Ptr fluidSdf() const override;

        //! Transfers velocity field from particles to grids.
        virtual void transferFromParticlesToGrids();

        //! Transfers velocity field from grids to particles.
        virtual void transferFromGridsToParticles();

        //! Moves particles.
        virtual void moveParticles(double timeIntervalInSeconds);

    private:
        size_t _signedDistanceFieldId;
        ParticleSystemData3Ptr _particles;
        ParticleEmitter3Ptr _particleEmitter;

        void extrapolateVelocityToAir();

        void buildSignedDistanceField();

        void updateParticleEmitter(double timeIntervalInSeconds);
    };

    //! Shared pointer type for the PicSolver3.
    typedef std::shared_ptr<PicSolver3> PicSolver3Ptr;


    //!
    //! \brief Front-end to create PicSolver3 objects step by step.
    //!
    class PicSolver3::Builder final
        : public GridFluidSolverBuilderBase3<PicSolver3::Builder> {
    public:
        //! Builds PicSolver3.
        PicSolver3 build() const;

        //! Builds shared pointer of PicSolver3 instance.
        PicSolver3Ptr makeShared() const;
    };

    PicSolver3::PicSolver3() : PicSolver3({ 1, 1, 1 }, { 1, 1, 1 }, { 0, 0, 0 }) {
    }

    PicSolver3::PicSolver3(
        const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& gridOrigin)
        : GridFluidSolver3(resolution, gridSpacing, gridOrigin) {
        auto grids = gridSystemData();
        _signedDistanceFieldId = grids->addScalarData(
            std::make_shared<CellCenteredScalarGrid3::Builder>(), kMaxD);
        _particles = std::make_shared<ParticleSystemData3>();
    }

    PicSolver3::~PicSolver3() {
    }

    ScalarGrid3Ptr PicSolver3::signedDistanceField() const {
        return gridSystemData()->scalarDataAt(_signedDistanceFieldId);
    }

    const ParticleSystemData3Ptr& PicSolver3::particleSystemData() const {
        return _particles;
    }

    const ParticleEmitter3Ptr& PicSolver3::particleEmitter() const {
        return _particleEmitter;
    }

    void PicSolver3::setParticleEmitter(const ParticleEmitter3Ptr& newEmitter) {
        _particleEmitter = newEmitter;
        newEmitter->setTarget(_particles);
    }

    void PicSolver3::onInitialize() {
        GridFluidSolver3::onInitialize();

        Timer timer;
        updateParticleEmitter(0.0);
        JET_INFO << "Update particle emitter took "
            << timer.durationInSeconds() << " seconds";
    }

    void PicSolver3::onBeginAdvanceTimeStep(double timeIntervalInSeconds) {
        UNUSED_VARIABLE(timeIntervalInSeconds);

        JET_INFO << "Number of PIC-type particles: "
            << _particles->numberOfParticles();

        Timer timer;
        updateParticleEmitter(timeIntervalInSeconds);
        JET_INFO << "Update particle emitter took "
            << timer.durationInSeconds() << " seconds";

        JET_INFO << "Number of PIC-type particles: "
            << _particles->numberOfParticles();

        timer.reset();
        transferFromParticlesToGrids();
        JET_INFO << "transferFromParticlesToGrids took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        buildSignedDistanceField();
        JET_INFO << "buildSignedDistanceField took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        extrapolateVelocityToAir();
        JET_INFO << "extrapolateVelocityToAir took "
            << timer.durationInSeconds() << " seconds";

        applyBoundaryCondition();
    }

    void PicSolver3::computeAdvection(double timeIntervalInSeconds) {
        Timer timer;
        extrapolateVelocityToAir();
        JET_INFO << "extrapolateVelocityToAir took "
            << timer.durationInSeconds() << " seconds";

        applyBoundaryCondition();

        timer.reset();
        transferFromGridsToParticles();
        JET_INFO << "transferFromGridsToParticles took "
            << timer.durationInSeconds() << " seconds";

        timer.reset();
        moveParticles(timeIntervalInSeconds);
        JET_INFO << "moveParticles took "
            << timer.durationInSeconds() << " seconds";
    }

    ScalarField3Ptr PicSolver3::fluidSdf() const {
        return signedDistanceField();
    }

    void PicSolver3::transferFromParticlesToGrids() {
        auto flow = gridSystemData()->velocity();
        auto positions = _particles->positions();
        auto velocities = _particles->velocities();
        size_t numberOfParticles = _particles->numberOfParticles();

        // Clear velocity to zero
        flow->fill(Vector3D());

        // Weighted-average velocity
        auto u = flow->uAccessor();
        auto v = flow->vAccessor();
        auto w = flow->wAccessor();
        Array3<double> uWeight(u.size());
        Array3<double> vWeight(v.size());
        Array3<double> wWeight(w.size());
        _uMarkers.resize(u.size());
        _vMarkers.resize(v.size());
        _wMarkers.resize(w.size());
        _uMarkers.set(0);
        _vMarkers.set(0);
        _wMarkers.set(0);
        LinearArraySampler3<double, double> uSampler(
            flow->uConstAccessor(),
            flow->gridSpacing(),
            flow->uOrigin());
        LinearArraySampler3<double, double> vSampler(
            flow->vConstAccessor(),
            flow->gridSpacing(),
            flow->vOrigin());
        LinearArraySampler3<double, double> wSampler(
            flow->wConstAccessor(),
            flow->gridSpacing(),
            flow->wOrigin());
        for (size_t i = 0; i < numberOfParticles; ++i) {
            std::array<Point3UI, 8> indices;
            std::array<double, 8> weights;

            uSampler.getCoordinatesAndWeights(positions[i], &indices, &weights);
            for (int j = 0; j < 8; ++j) {
                u(indices[j]) += velocities[i].x * weights[j];
                uWeight(indices[j]) += weights[j];
                _uMarkers(indices[j]) = 1;
            }

            vSampler.getCoordinatesAndWeights(positions[i], &indices, &weights);
            for (int j = 0; j < 8; ++j) {
                v(indices[j]) += velocities[i].y * weights[j];
                vWeight(indices[j]) += weights[j];
                _vMarkers(indices[j]) = 1;
            }

            wSampler.getCoordinatesAndWeights(positions[i], &indices, &weights);
            for (int j = 0; j < 8; ++j) {
                w(indices[j]) += velocities[i].z * weights[j];
                wWeight(indices[j]) += weights[j];
                _wMarkers(indices[j]) = 1;
            }
        }

        uWeight.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (uWeight(i, j, k) > 0.0) {
                u(i, j, k) /= uWeight(i, j, k);
            }
            });
        vWeight.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (vWeight(i, j, k) > 0.0) {
                v(i, j, k) /= vWeight(i, j, k);
            }
            });
        wWeight.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            if (wWeight(i, j, k) > 0.0) {
                w(i, j, k) /= wWeight(i, j, k);
            }
            });
    }

    void PicSolver3::transferFromGridsToParticles() {
        auto flow = gridSystemData()->velocity();
        auto positions = _particles->positions();
        auto velocities = _particles->velocities();
        size_t numberOfParticles = _particles->numberOfParticles();

        parallelFor(kZeroSize, numberOfParticles, [&](size_t i) {
            velocities[i] = flow->sample(positions[i]);
            });
    }

    void PicSolver3::moveParticles(double timeIntervalInSeconds) {
        auto flow = gridSystemData()->velocity();
        auto positions = _particles->positions();
        auto velocities = _particles->velocities();
        size_t numberOfParticles = _particles->numberOfParticles();
        int domainBoundaryFlag = closedDomainBoundaryFlag();
        BoundingBox3D boundingBox = flow->boundingBox();

        parallelFor(kZeroSize, numberOfParticles, [&](size_t i) {
            Vector3D pt0 = positions[i];
            Vector3D pt1 = pt0;
            Vector3D vel = velocities[i];

            // Adaptive time-stepping
            unsigned int numSubSteps
                = static_cast<unsigned int>(std::max(maxCfl(), 1.0));
            double dt = timeIntervalInSeconds / numSubSteps;
            for (unsigned int t = 0; t < numSubSteps; ++t) {
                Vector3D vel0 = flow->sample(pt0);

                // Mid-point rule
                Vector3D midPt = pt0 + 0.5 * dt * vel0;
                Vector3D midVel = flow->sample(midPt);
                pt1 = pt0 + dt * midVel;

                pt0 = pt1;
            }

            if ((domainBoundaryFlag & kDirectionLeft)
                && pt1.x <= boundingBox.lowerCorner.x) {
                pt1.x = boundingBox.lowerCorner.x;
                vel.x = 0.0;
            }
            if ((domainBoundaryFlag & kDirectionRight)
                && pt1.x >= boundingBox.upperCorner.x) {
                pt1.x = boundingBox.upperCorner.x;
                vel.x = 0.0;
            }
            if ((domainBoundaryFlag & kDirectionDown)
                && pt1.y <= boundingBox.lowerCorner.y) {
                pt1.y = boundingBox.lowerCorner.y;
                vel.y = 0.0;
            }
            if ((domainBoundaryFlag & kDirectionUp)
                && pt1.y >= boundingBox.upperCorner.y) {
                pt1.y = boundingBox.upperCorner.y;
                vel.y = 0.0;
            }
            if ((domainBoundaryFlag & kDirectionBack)
                && pt1.z <= boundingBox.lowerCorner.z) {
                pt1.z = boundingBox.lowerCorner.z;
                vel.z = 0.0;
            }
            if ((domainBoundaryFlag & kDirectionFront)
                && pt1.z >= boundingBox.upperCorner.z) {
                pt1.z = boundingBox.upperCorner.z;
                vel.z = 0.0;
            }

            positions[i] = pt1;
            velocities[i] = vel;
            });

        Collider3Ptr col = collider();
        if (col != nullptr) {
            parallelFor(
                kZeroSize,
                numberOfParticles,
                [&](size_t i) {
                    col->resolveCollision(
                        0.0,
                        0.0,
                        &positions[i],
                        &velocities[i]);
                });
        }
    }

    void PicSolver3::extrapolateVelocityToAir() {
        auto vel = gridSystemData()->velocity();
        auto u = vel->uAccessor();
        auto v = vel->vAccessor();
        auto w = vel->wAccessor();

        unsigned int depth = static_cast<unsigned int>(std::ceil(maxCfl()));
        extrapolateToRegion(vel->uConstAccessor(), _uMarkers, depth, u);
        extrapolateToRegion(vel->vConstAccessor(), _vMarkers, depth, v);
        extrapolateToRegion(vel->wConstAccessor(), _wMarkers, depth, w);
    }

    void PicSolver3::buildSignedDistanceField() {
        auto sdf = signedDistanceField();
        auto sdfPos = sdf->dataPosition();
        double maxH = max3(
            sdf->gridSpacing().x, sdf->gridSpacing().y, sdf->gridSpacing().z);
        double radius = 1.2 * maxH / std::sqrt(2.0);
        double sdfBandRadius = 2.0 * radius;

        _particles->buildNeighborSearcher(2 * radius);
        auto searcher = _particles->neighborSearcher();
        sdf->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            Vector3D pt = sdfPos(i, j, k);
            double minDist = sdfBandRadius;
            searcher->forEachNearbyPoint(
                pt, sdfBandRadius, [&](size_t, const Vector3D& x) {
                    minDist = std::min(minDist, pt.distanceTo(x));
                });
            (*sdf)(i, j, k) = minDist - radius;
            });

        extrapolateIntoCollider(sdf.get());
    }

    void PicSolver3::updateParticleEmitter(double timeIntervalInSeconds) {
        if (_particleEmitter != nullptr) {
            _particleEmitter->update(currentTimeInSeconds(), timeIntervalInSeconds);
        }
    }

    PicSolver3::Builder PicSolver3::builder() {
        return Builder();
    }


    PicSolver3 PicSolver3::Builder::build() const {
        return PicSolver3(_resolution, getGridSpacing(), _gridOrigin);
    }

    PicSolver3Ptr PicSolver3::Builder::makeShared() const {
        return std::shared_ptr<PicSolver3>(
            new PicSolver3(
                _resolution,
                getGridSpacing(),
                _gridOrigin),
            [](PicSolver3* obj) {
                delete obj;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_PIC_SOLVER3_H_