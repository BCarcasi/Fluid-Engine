#ifndef INCLUDE_JET_VOLUME_PARTICLE_EMITTER3_H_
#define INCLUDE_JET_VOLUME_PARTICLE_EMITTER3_H_

#include "bounding_box3.h"
#include "implicit_surface3.h"
#include "particle_emitter3.h"
#include "point_generator3.h"

#include <limits>
#include <memory>
#include <random>

#include "pch.h"

#include "bcc_lattice_point_generator.h"
#include "point_hash_grid_searcher3.h"
#include "samplers.h"
#include "surface_to_implicit3.h"


namespace jet {

    //!
    //! \brief 3-D volumetric particle emitter.
    //!
    //! This class emits particles from volumetric geometry.
    //!
    class VolumeParticleEmitter3 final : public ParticleEmitter3 {
    public:
        class Builder;

        //!
        //! Constructs an emitter that spawns particles from given implicit surface
        //! which defines the volumetric geometry. Provided bounding box limits
        //! the particle generation region.
        //!
        //! \param[in]  implicitSurface         The implicit surface.
        //! \param[in]  maxRegion               The max region.
        //! \param[in]  spacing                 The spacing between particles.
        //! \param[in]  initialVel              The initial velocity.
        //! \param[in]  linearVel               The linear velocity of the emitter.
        //! \param[in]  angularVel              The angular velocity of the emitter.
        //! \param[in]  maxNumberOfParticles    The max number of particles to be
        //!                                     emitted.
        //! \param[in]  jitter                  The jitter amount between 0 and 1.
        //! \param[in]  isOneShot               True if emitter gets disabled after one shot.
        //! \param[in]  allowOverlapping        True if particles can be overlapped.
        //! \param[in]  seed                    The random seed.
        //!
        VolumeParticleEmitter3(
            const ImplicitSurface3Ptr& implicitSurface,
            const BoundingBox3D& maxRegion,
            double spacing,
            const Vector3D& initialVel = Vector3D(),
            const Vector3D& linearVel = Vector3D(),
            const Vector3D& angularVel = Vector3D(),
            size_t maxNumberOfParticles = kMaxSize,
            double jitter = 0.0,
            bool isOneShot = true,
            bool allowOverlapping = false,
            uint32_t seed = 0);

        //!
        //! \brief      Sets the point generator.
        //!
        //! This function sets the point generator that defines the pattern of the
        //! point distribution within the volume.
        //!
        //! \param[in]  newPointsGen The new points generator.
        //!
        void setPointGenerator(const PointGenerator3Ptr& newPointsGen);

        //! Returns source surface.
        const ImplicitSurface3Ptr& surface() const;

        //! Sets the source surface.
        void setSurface(const ImplicitSurface3Ptr& newSurface);

        //! Returns max particle gen region.
        const BoundingBox3D& maxRegion() const;

        //! Sets the max particle gen region.
        void setMaxRegion(const BoundingBox3D& newBox);

        //! Returns jitter amount.
        double jitter() const;

        //! Sets jitter amount between 0 and 1.
        void setJitter(double newJitter);

        //! Returns true if particles should be emitted just once.
        bool isOneShot() const;

        //!
        //! \brief      Sets the flag to true if particles are emitted just once.
        //!
        //! If true is set, the emitter will generate particles only once even after
        //! multiple emit calls. If false, it will keep generating particles from
        //! the volumetric geometry. Default value is true.
        //!
        //! \param[in]  newValue True if particles should be emitted just once.
        //!
        void setIsOneShot(bool newValue);

        //! Returns true if particles can be overlapped.
        bool allowOverlapping() const;

        //!
        //! \brief      Sets the flag to true if particles can overlap each other.
        //!
        //! If true is set, the emitter will generate particles even if the new
        //! particles can find existing nearby particles within the particle
        //! spacing.
        //!
        //! \param[in]  newValue True if particles can be overlapped.
        //!
        void setAllowOverlapping(bool newValue);

        //! Returns max number of particles to be emitted.
        size_t maxNumberOfParticles() const;

        //! Sets the max number of particles to be emitted.
        void setMaxNumberOfParticles(size_t newMaxNumberOfParticles);

        //! Returns the spacing between particles.
        double spacing() const;

        //! Sets the spacing between particles.
        void setSpacing(double newSpacing);

        //! Sets the initial velocity of the particles.
        Vector3D initialVelocity() const;

        //! Returns the initial velocity of the particles.
        void setInitialVelocity(const Vector3D& newInitialVel);

        //! Returns the linear velocity of the emitter.
        Vector3D linearVelocity() const;

        //! Sets the linear velocity of the emitter.
        void setLinearVelocity(const Vector3D& newLinearVel);

        //! Returns the angular velocity of the emitter.
        Vector3D angularVelocity() const;

        //! Sets the linear velocity of the emitter.
        void setAngularVelocity(const Vector3D& newAngularVel);

        //! Returns builder fox VolumeParticleEmitter3.
        static Builder builder();

    private:
        std::mt19937 _rng;

        ImplicitSurface3Ptr _implicitSurface;
        BoundingBox3D _bounds;
        double _spacing;
        Vector3D _initialVel;
        Vector3D _linearVel;
        Vector3D _angularVel;
        PointGenerator3Ptr _pointsGen;

        size_t _maxNumberOfParticles = kMaxSize;
        size_t _numberOfEmittedParticles = 0;

        double _jitter = 0.0;
        bool _isOneShot = true;
        bool _allowOverlapping = false;

        //!
        //! \brief      Emits particles to the particle system data.
        //!
        //! \param[in]  currentTimeInSeconds    Current simulation time.
        //! \param[in]  timeIntervalInSeconds   The time-step interval.
        //!
        void onUpdate(
            double currentTimeInSeconds,
            double timeIntervalInSeconds) override;

        void emit(
            const ParticleSystemData3Ptr& particles,
            Array1<Vector3D>* newPositions,
            Array1<Vector3D>* newVelocities);

        double random();

        Vector3D velocityAt(const Vector3D& point) const;
    };

    //! Shared pointer for the VolumeParticleEmitter3 type.
    typedef std::shared_ptr<VolumeParticleEmitter3> VolumeParticleEmitter3Ptr;


    //!
    //! \brief Front-end to create VolumeParticleEmitter3 objects step by step.
    //!
    class VolumeParticleEmitter3::Builder final {
    public:
        //! Returns builder with implicit surface defining volume shape.
        Builder& withImplicitSurface(const ImplicitSurface3Ptr& implicitSurface);

        //! Returns builder with surface defining volume shape.
        Builder& withSurface(const Surface3Ptr& surface);

        //! Returns builder with max region.
        Builder& withMaxRegion(const BoundingBox3D& bounds);

        //! Returns builder with spacing.
        Builder& withSpacing(double spacing);

        //! Returns builder with initial velocity.
        Builder& withInitialVelocity(const Vector3D& initialVel);

        //! Returns builder with linear velocity.
        Builder& withLinearVelocity(const Vector3D& linearVel);

        //! Returns builder with angular velocity.
        Builder& withAngularVelocity(const Vector3D& angularVel);

        //! Returns builder with max number of particles.
        Builder& withMaxNumberOfParticles(size_t maxNumberOfParticles);

        //! Returns builder with jitter amount.
        Builder& withJitter(double jitter);

        //! Returns builder with one-shot flag.
        Builder& withIsOneShot(bool isOneShot);

        //! Returns builder with overlapping flag.
        Builder& withAllowOverlapping(bool allowOverlapping);

        //! Returns builder with random seed.
        Builder& withRandomSeed(uint32_t seed);

        //! Builds VolumeParticleEmitter3.
        VolumeParticleEmitter3 build() const;

        //! Builds shared pointer of VolumeParticleEmitter3 instance.
        VolumeParticleEmitter3Ptr makeShared() const;

    private:
        ImplicitSurface3Ptr _implicitSurface;
        bool _isBoundSet = false;
        BoundingBox3D _bounds;
        double _spacing = 0.1;
        Vector3D _initialVel;
        Vector3D _linearVel;
        Vector3D _angularVel;
        size_t _maxNumberOfParticles = kMaxSize;
        double _jitter = 0.0;
        bool _isOneShot = true;
        bool _allowOverlapping = false;
        uint32_t _seed = 0;
    };

    static const size_t kDefaultHashGridResolutionParticleEmitter3 = 64;

    VolumeParticleEmitter3::VolumeParticleEmitter3(
        const ImplicitSurface3Ptr& implicitSurface, const BoundingBox3D& maxRegion,
        double spacing, const Vector3D& initialVel, const Vector3D& linearVel,
        const Vector3D& angularVel, size_t maxNumberOfParticles, double jitter,
        bool isOneShot, bool allowOverlapping, uint32_t seed)
        : _rng(seed),
        _implicitSurface(implicitSurface),
        _bounds(maxRegion),
        _spacing(spacing),
        _initialVel(initialVel),
        _linearVel(linearVel),
        _angularVel(angularVel),
        _maxNumberOfParticles(maxNumberOfParticles),
        _jitter(jitter),
        _isOneShot(isOneShot),
        _allowOverlapping(allowOverlapping) {
        _pointsGen = std::make_shared<BccLatticePointGenerator>();
    }

    void VolumeParticleEmitter3::onUpdate(double currentTimeInSeconds,
        double timeIntervalInSeconds) {
        UNUSED_VARIABLE(currentTimeInSeconds);
        UNUSED_VARIABLE(timeIntervalInSeconds);

        auto particles = target();

        if (particles == nullptr) {
            return;
        }

        if (!isEnabled()) {
            return;
        }

        Array1<Vector3D> newPositions;
        Array1<Vector3D> newVelocities;

        emit(particles, &newPositions, &newVelocities);

        particles->addParticles(newPositions, newVelocities);

        if (_isOneShot) {
            setIsEnabled(false);
        }
    }

    void VolumeParticleEmitter3::emit(const ParticleSystemData3Ptr& particles,
        Array1<Vector3D>* newPositions,
        Array1<Vector3D>* newVelocities) {
        if (!_implicitSurface) {
            return;
        }

        _implicitSurface->updateQueryEngine();

        BoundingBox3D region = _bounds;
        if (_implicitSurface->isBounded()) {
            BoundingBox3D surfaceBBox = _implicitSurface->boundingBox();
            region.lowerCorner = max(region.lowerCorner, surfaceBBox.lowerCorner);
            region.upperCorner = min(region.upperCorner, surfaceBBox.upperCorner);
        }

        // Reserving more space for jittering
        const double j = jitter();
        const double maxJitterDist = 0.5 * j * _spacing;
        size_t numNewParticles = 0;

        if (_allowOverlapping || _isOneShot) {
            _pointsGen->forEachPoint(region, _spacing, [&](const Vector3D& point) {
                Vector3D randomDir = uniformSampleSphere(random(), random());
                Vector3D offset = maxJitterDist * randomDir;
                Vector3D candidate = point + offset;
                if (_implicitSurface->signedDistance(candidate) <= 0.0) {
                    if (_numberOfEmittedParticles < _maxNumberOfParticles) {
                        newPositions->append(candidate);
                        ++_numberOfEmittedParticles;
                        ++numNewParticles;
                    }
                    else {
                        return false;
                    }
                }

                return true;
                });
        }
        else {
            // Use serial hash grid searcher for continuous update.
            PointHashGridSearcher3 neighborSearcher(
                Size3(kDefaultHashGridResolutionParticleEmitter3, kDefaultHashGridResolutionParticleEmitter3,
                    kDefaultHashGridResolutionParticleEmitter3),
                2.0 * _spacing);
            if (!_allowOverlapping) {
                neighborSearcher.build(particles->positions());
            }

            _pointsGen->forEachPoint(region, _spacing, [&](const Vector3D& point) {
                Vector3D randomDir = uniformSampleSphere(random(), random());
                Vector3D offset = maxJitterDist * randomDir;
                Vector3D candidate = point + offset;
                if (_implicitSurface->isInside(candidate) &&
                    (!_allowOverlapping &&
                        !neighborSearcher.hasNearbyPoint(candidate, _spacing))) {
                    if (_numberOfEmittedParticles < _maxNumberOfParticles) {
                        newPositions->append(candidate);
                        neighborSearcher.add(candidate);
                        ++_numberOfEmittedParticles;
                        ++numNewParticles;
                    }
                    else {
                        return false;
                    }
                }

                return true;
                });
        }

        JET_INFO << "Number of newly generated particles: " << numNewParticles;
        JET_INFO << "Number of total generated particles: "
            << _numberOfEmittedParticles;

        newVelocities->resize(newPositions->size());
        newVelocities->parallelForEachIndex([&](size_t i) {
            (*newVelocities)[i] = velocityAt((*newPositions)[i]);
            });
    }

    void VolumeParticleEmitter3::setPointGenerator(
        const PointGenerator3Ptr& newPointsGen) {
        _pointsGen = newPointsGen;
    }

    const ImplicitSurface3Ptr& VolumeParticleEmitter3::surface() const {
        return _implicitSurface;
    }

    void VolumeParticleEmitter3::setSurface(const ImplicitSurface3Ptr& newSurface) {
        _implicitSurface = newSurface;
    }

    const BoundingBox3D& VolumeParticleEmitter3::maxRegion() const {
        return _bounds;
    }

    void VolumeParticleEmitter3::setMaxRegion(const BoundingBox3D& newMaxRegion) {
        _bounds = newMaxRegion;
    }

    double VolumeParticleEmitter3::jitter() const { return _jitter; }

    void VolumeParticleEmitter3::setJitter(double newJitter) {
        _jitter = clamp(newJitter, 0.0, 1.0);
    }

    bool VolumeParticleEmitter3::isOneShot() const { return _isOneShot; }

    void VolumeParticleEmitter3::setIsOneShot(bool newValue) {
        _isOneShot = newValue;
    }

    bool VolumeParticleEmitter3::allowOverlapping() const {
        return _allowOverlapping;
    }

    void VolumeParticleEmitter3::setAllowOverlapping(bool newValue) {
        _allowOverlapping = newValue;
    }

    size_t VolumeParticleEmitter3::maxNumberOfParticles() const {
        return _maxNumberOfParticles;
    }

    void VolumeParticleEmitter3::setMaxNumberOfParticles(
        size_t newMaxNumberOfParticles) {
        _maxNumberOfParticles = newMaxNumberOfParticles;
    }

    double VolumeParticleEmitter3::spacing() const { return _spacing; }

    void VolumeParticleEmitter3::setSpacing(double newSpacing) {
        _spacing = newSpacing;
    }

    Vector3D VolumeParticleEmitter3::initialVelocity() const { return _initialVel; }

    void VolumeParticleEmitter3::setInitialVelocity(const Vector3D& newInitialVel) {
        _initialVel = newInitialVel;
    }

    Vector3D VolumeParticleEmitter3::linearVelocity() const { return _linearVel; }

    void VolumeParticleEmitter3::setLinearVelocity(const Vector3D& newLinearVel) {
        _linearVel = newLinearVel;
    }

    Vector3D VolumeParticleEmitter3::angularVelocity() const { return _angularVel; }

    void VolumeParticleEmitter3::setAngularVelocity(const Vector3D& newAngularVel) {
        _angularVel = newAngularVel;
    }

    double VolumeParticleEmitter3::random() {
        std::uniform_real_distribution<> d(0.0, 1.0);
        return d(_rng);
    }

    Vector3D VolumeParticleEmitter3::velocityAt(const Vector3D& point) const {
        Vector3D r = point - _implicitSurface->transform.translation();
        return _linearVel + _angularVel.cross(r) + _initialVel;
    }

    VolumeParticleEmitter3::Builder VolumeParticleEmitter3::builder() {
        return Builder();
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withImplicitSurface(
            const ImplicitSurface3Ptr& implicitSurface) {
        _implicitSurface = implicitSurface;
        if (!_isBoundSet) {
            _bounds = _implicitSurface->boundingBox();
        }
        return *this;
    }

    VolumeParticleEmitter3::Builder& VolumeParticleEmitter3::Builder::withSurface(
        const Surface3Ptr& surface) {
        _implicitSurface = std::make_shared<SurfaceToImplicit3>(surface);
        if (!_isBoundSet) {
            _bounds = surface->boundingBox();
        }
        return *this;
    }

    VolumeParticleEmitter3::Builder& VolumeParticleEmitter3::Builder::withMaxRegion(
        const BoundingBox3D& bounds) {
        _bounds = bounds;
        _isBoundSet = true;
        return *this;
    }

    VolumeParticleEmitter3::Builder& VolumeParticleEmitter3::Builder::withSpacing(
        double spacing) {
        _spacing = spacing;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withInitialVelocity(
            const Vector3D& initialVel) {
        _initialVel = initialVel;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withLinearVelocity(const Vector3D& linearVel) {
        _linearVel = linearVel;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withAngularVelocity(
            const Vector3D& angularVel) {
        _angularVel = angularVel;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withMaxNumberOfParticles(
            size_t maxNumberOfParticles) {
        _maxNumberOfParticles = maxNumberOfParticles;
        return *this;
    }

    VolumeParticleEmitter3::Builder& VolumeParticleEmitter3::Builder::withJitter(
        double jitter) {
        _jitter = jitter;
        return *this;
    }

    VolumeParticleEmitter3::Builder& VolumeParticleEmitter3::Builder::withIsOneShot(
        bool isOneShot) {
        _isOneShot = isOneShot;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withAllowOverlapping(bool allowOverlapping) {
        _allowOverlapping = allowOverlapping;
        return *this;
    }

    VolumeParticleEmitter3::Builder&
        VolumeParticleEmitter3::Builder::withRandomSeed(uint32_t seed) {
        _seed = seed;
        return *this;
    }

    VolumeParticleEmitter3 VolumeParticleEmitter3::Builder::build() const {
        return VolumeParticleEmitter3(_implicitSurface, _bounds, _spacing,
            _initialVel, _linearVel, _angularVel,
            _maxNumberOfParticles, _jitter, _isOneShot,
            _allowOverlapping, _seed);
    }

    VolumeParticleEmitter3Ptr VolumeParticleEmitter3::Builder::makeShared() const {
        return std::shared_ptr<VolumeParticleEmitter3>(
            new VolumeParticleEmitter3(_implicitSurface, _bounds, _spacing,
                _initialVel, _linearVel, _angularVel,
                _maxNumberOfParticles, _jitter, _isOneShot,
                _allowOverlapping),
            [](VolumeParticleEmitter3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VOLUME_PARTICLE_EMITTER3_H_