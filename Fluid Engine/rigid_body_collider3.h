#ifndef INCLUDE_JET_RIGID_BODY_COLLIDER3_H_
#define INCLUDE_JET_RIGID_BODY_COLLIDER3_H_

#include "collider3.h"
#include "quaternion.h"

#include "pch.h"

namespace jet {

    //!
    //! \brief 3-D rigid body collider class.
    //!
    //! This class implements 3-D rigid body collider. The collider can only take
    //! rigid body motion with linear and rotational velocities.
    //!
    class RigidBodyCollider3 final : public Collider3 {
    public:
        class Builder;

        //! Linear velocity of the rigid body.
        Vector3D linearVelocity;

        //! Angular velocity of the rigid body.
        Vector3D angularVelocity;

        //! Constructs a collider with a surface.
        explicit RigidBodyCollider3(const Surface3Ptr& surface);

        //! Constructs a collider with a surface and other parameters.
        RigidBodyCollider3(
            const Surface3Ptr& surface,
            const Vector3D& linearVelocity,
            const Vector3D& angularVelocity);

        //! Returns the velocity of the collider at given \p point.
        Vector3D velocityAt(const Vector3D& point) const override;

        //! Returns builder fox RigidBodyCollider3.
        static Builder builder();
    };

    //! Shared pointer for the RigidBodyCollider3 type.
    typedef std::shared_ptr<RigidBodyCollider3> RigidBodyCollider3Ptr;


    //!
    //! \brief Front-end to create RigidBodyCollider3 objects step by step.
    //!
    class RigidBodyCollider3::Builder final {
    public:
        //! Returns builder with surface.
        Builder& withSurface(const Surface3Ptr& surface);

        //! Returns builder with linear velocity.
        Builder& withLinearVelocity(const Vector3D& linearVelocity);

        //! Returns builder with angular velocity.
        Builder& withAngularVelocity(const Vector3D& angularVelocity);

        //! Builds RigidBodyCollider3.
        RigidBodyCollider3 build() const;

        //! Builds shared pointer of RigidBodyCollider3 instance.
        RigidBodyCollider3Ptr makeShared() const;

    private:
        Surface3Ptr _surface;
        Vector3D _linearVelocity{ 0, 0, 0 };
        Vector3D _angularVelocity{ 0, 0, 0 };
    };
    RigidBodyCollider3::RigidBodyCollider3(const Surface3Ptr& surface) {
        setSurface(surface);
    }

    RigidBodyCollider3::RigidBodyCollider3(
        const Surface3Ptr& surface,
        const Vector3D& linearVelocity_,
        const Vector3D& angularVelocity_)
        : linearVelocity(linearVelocity_)
        , angularVelocity(angularVelocity_) {
        setSurface(surface);
    }

    Vector3D RigidBodyCollider3::velocityAt(const Vector3D& point) const {
        Vector3D r = point - surface()->transform.translation();
        return linearVelocity + angularVelocity.cross(r);
    }

    RigidBodyCollider3::Builder RigidBodyCollider3::builder() {
        return Builder();
    }

    RigidBodyCollider3::Builder&
        RigidBodyCollider3::Builder::withSurface(const Surface3Ptr& surface) {
        _surface = surface;
        return *this;
    }

    RigidBodyCollider3::Builder&
        RigidBodyCollider3::Builder::withLinearVelocity(
            const Vector3D& linearVelocity) {
        _linearVelocity = linearVelocity;
        return *this;
    }

    RigidBodyCollider3::Builder&
        RigidBodyCollider3::Builder::withAngularVelocity(
            const Vector3D& angularVelocity) {
        _angularVelocity = angularVelocity;
        return *this;
    }

    RigidBodyCollider3 RigidBodyCollider3::Builder::build() const {
        return RigidBodyCollider3(
            _surface,
            _linearVelocity,
            _angularVelocity);
    }

    RigidBodyCollider3Ptr RigidBodyCollider3::Builder::makeShared() const {
        return std::shared_ptr<RigidBodyCollider3>(
            new RigidBodyCollider3(
                _surface,
                _linearVelocity,
                _angularVelocity),
            [](RigidBodyCollider3* obj) {
                delete obj;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_RIGID_BODY_COLLIDER3_H_