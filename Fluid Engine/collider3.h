#ifndef INCLUDE_JET_COLLIDER3_H_
#define INCLUDE_JET_COLLIDER3_H_

#include "surface3.h"
#include <functional>

#include "pch.h"


#include <algorithm>

namespace jet {

    //!
    //! \brief Abstract base class for generic collider object.
    //!
    //! This class contains basic interfaces for colliders. Most of the
    //! functionalities are implemented within this class, except the member
    //! function Collider3::velocityAt. This class also let the subclasses to
    //! provide a Surface3 instance to define collider surface using
    //! Collider3::setSurface function.
    //!
    class Collider3 {
    public:
        //!
        //! \brief Callback function type for update calls.
        //!
        //! This type of callback function will take the collider pointer, current
        //! time, and time interval in seconds.
        //!
        typedef std::function<void(Collider3*, double, double)>
            OnBeginUpdateCallback;

        //! Default constructor.
        Collider3();

        //! Default destructor.
        virtual ~Collider3();

        //! Returns the velocity of the collider at given \p point.
        virtual Vector3D velocityAt(const Vector3D& point) const = 0;

        //!
        //! Resolves collision for given point.
        //!
        //! \param radius Radius of the colliding point.
        //! \param restitutionCoefficient Defines the restitution effect.
        //! \param position Input and output position of the point.
        //! \param position Input and output velocity of the point.
        //!
        void resolveCollision(
            double radius,
            double restitutionCoefficient,
            Vector3D* position,
            Vector3D* velocity);

        //! Returns friction coefficent.
        double frictionCoefficient() const;

        //!
        //! \brief Sets the friction coefficient.
        //!
        //! This function assigns the friction coefficient to the collider. Any
        //! negative inputs will be clamped to zero.
        //!
        void setFrictionCoefficient(double newFrictionCoeffient);

        //! Returns the surface instance.
        const Surface3Ptr& surface() const;

        //! Updates the collider state.
        void update(double currentTimeInSeconds, double timeIntervalInSeconds);

        //!
        //! \brief      Sets the callback function to be called when
        //!             Collider2::update function is invoked.
        //!
        //! The callback function takes current simulation time in seconds unit. Use
        //! this callback to track any motion or state changes related to this
        //! collider.
        //!
        //! \param[in]  callback The callback function.
        //!
        void setOnBeginUpdateCallback(const OnBeginUpdateCallback& callback);

    protected:
        //! Internal query result structure.
        struct ColliderQueryResult final {
            double distance;
            Vector3D point;
            Vector3D normal;
            Vector3D velocity;
        };

        //! Assigns the surface instance from the subclass.
        void setSurface(const Surface3Ptr& newSurface);

        //! Outputs closest point's information.
        void getClosestPoint(
            const Surface3Ptr& surface,
            const Vector3D& queryPoint,
            ColliderQueryResult* result) const;

        //! Returns true if given point is in the opposite side of the surface.
        bool isPenetrating(
            const ColliderQueryResult& colliderPoint,
            const Vector3D& position,
            double radius);

    private:
        Surface3Ptr _surface;
        double _frictionCoeffient = 0.0;
        OnBeginUpdateCallback _onUpdateCallback;
    };

    //! Shared pointer type for the Collider2.
    typedef std::shared_ptr<Collider3> Collider3Ptr;

    Collider3::Collider3() {}

    Collider3::~Collider3() {}

    void Collider3::resolveCollision(double radius, double restitutionCoefficient,
        Vector3D* newPosition, Vector3D* newVelocity) {
        JET_ASSERT(_surface);

        if (!_surface->isValidGeometry()) {
            return;
        }

        ColliderQueryResult colliderPoint;

        getClosestPoint(_surface, *newPosition, &colliderPoint);

        // Check if the new position is penetrating the surface
        if (isPenetrating(colliderPoint, *newPosition, radius)) {
            // Target point is the closest non-penetrating position from the
            // new position.
            Vector3D targetNormal = colliderPoint.normal;
            Vector3D targetPoint = colliderPoint.point + radius * targetNormal;
            Vector3D colliderVelAtTargetPoint = colliderPoint.velocity;

            // Get new candidate relative velocity from the target point.
            Vector3D relativeVel = *newVelocity - colliderVelAtTargetPoint;
            double normalDotRelativeVel = targetNormal.dot(relativeVel);
            Vector3D relativeVelN = normalDotRelativeVel * targetNormal;
            Vector3D relativeVelT = relativeVel - relativeVelN;

            // Check if the velocity is facing opposite direction of the surface
            // normal
            if (normalDotRelativeVel < 0.0) {
                // Apply restitution coefficient to the surface normal component of
                // the velocity
                Vector3D deltaRelativeVelN =
                    (-restitutionCoefficient - 1.0) * relativeVelN;
                relativeVelN *= -restitutionCoefficient;

                // Apply friction to the tangential component of the velocity
                // From Bridson et al., Robust Treatment of Collisions, Contact and
                // Friction for Cloth Animation, 2002
                // http://graphics.stanford.edu/papers/cloth-sig02/cloth.pdf
                if (relativeVelT.lengthSquared() > 0.0) {
                    double frictionScale = std::max(
                        1.0 - _frictionCoeffient * deltaRelativeVelN.length() /
                        relativeVelT.length(),
                        0.0);
                    relativeVelT *= frictionScale;
                }

                // Reassemble the components
                *newVelocity =
                    relativeVelN + relativeVelT + colliderVelAtTargetPoint;
            }

            // Geometric fix
            *newPosition = targetPoint;
        }
    }

    double Collider3::frictionCoefficient() const { return _frictionCoeffient; }

    void Collider3::setFrictionCoefficient(double newFrictionCoeffient) {
        _frictionCoeffient = std::max(newFrictionCoeffient, 0.0);
    }

    const Surface3Ptr& Collider3::surface() const { return _surface; }

    void Collider3::setSurface(const Surface3Ptr& newSurface) {
        _surface = newSurface;
    }

    void Collider3::getClosestPoint(const Surface3Ptr& surface,
        const Vector3D& queryPoint,
        ColliderQueryResult* result) const {
        result->distance = surface->closestDistance(queryPoint);
        result->point = surface->closestPoint(queryPoint);
        result->normal = surface->closestNormal(queryPoint);
        result->velocity = velocityAt(queryPoint);
    }

    bool Collider3::isPenetrating(const ColliderQueryResult& colliderPoint,
        const Vector3D& position, double radius) {
        // If the new candidate position of the particle is inside
        // the volume defined by the surface OR the new distance to the surface is
        // less than the particle's radius, this particle is in colliding state.
        return _surface->isInside(position) || colliderPoint.distance < radius;
    }

    void Collider3::update(double currentTimeInSeconds,
        double timeIntervalInSeconds) {
        JET_ASSERT(_surface);

        if (!_surface->isValidGeometry()) {
            return;
        }

        _surface->updateQueryEngine();

        if (_onUpdateCallback) {
            _onUpdateCallback(this, currentTimeInSeconds, timeIntervalInSeconds);
        }
    }

    void Collider3::setOnBeginUpdateCallback(
        const OnBeginUpdateCallback& callback) {
        _onUpdateCallback = callback;
    }

}  // namespace jet

#endif  // INCLUDE_JET_COLLIDER3_H_