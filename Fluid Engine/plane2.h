#ifndef INCLUDE_JET_PLANE2_H_
#define INCLUDE_JET_PLANE2_H_

#include "surface2.h"

#include "pch.h"
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief 2-D plane geometry.
    //!
    //! This class represents 2-D plane geometry which extends Surface2 by
    //! overriding surface-related queries.
    //!
    class Plane2 final : public Surface2 {
    public:
        class Builder;

        //! Plane normal.
        Vector2D normal = Vector2D(0, 1);

        //! Point that lies on the plane.
        Vector2D point;

        //! Constructs a plane that crosses (0, 0) with surface normal (0, 1).
        Plane2(
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        //! Constructs a plane that cross \p point with surface normal \p normal.
        Plane2(
            const Vector2D& normal,
            const Vector2D& point,
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        //! Copy constructor.
        Plane2(const Plane2& other);

        //! Returns true if bounding box can be defined.
        bool isBounded() const override;

        //! Returns builder fox Plane2.
        static Builder builder();

    private:
        Vector2D closestPointLocal(const Vector2D& otherPoint) const override;

        bool intersectsLocal(const Ray2D& ray) const override;

        BoundingBox2D boundingBoxLocal() const override;

        Vector2D closestNormalLocal(const Vector2D& otherPoint) const override;

        SurfaceRayIntersection2 closestIntersectionLocal(
            const Ray2D& ray) const override;
    };

    //! Shared pointer for the Plane2 type.
    typedef std::shared_ptr<Plane2> Plane2Ptr;


    //!
    //! \brief Front-end to create Plane2 objects step by step.
    //!
    class Plane2::Builder final : public SurfaceBuilderBase2<Plane2::Builder> {
    public:
        //! Returns builder with plane normal.
        Builder& withNormal(const Vector2D& normal);

        //! Returns builder with point on the plane.
        Builder& withPoint(const Vector2D& point);

        //! Builds Plane2.
        Plane2 build() const;

        //! Builds shared pointer of Plane2 instance.
        Plane2Ptr makeShared() const;

    private:
        Vector2D _normal{ 0, 1 };
        Vector2D _point{ 0, 0 };
    };
    Plane2::Plane2(const Transform2& transform_, bool isNormalFlipped_)
        : Surface2(transform_, isNormalFlipped_) {}

    Plane2::Plane2(const Vector2D& normal_, const Vector2D& point_,
        const Transform2& transform_, bool isNormalFlipped_)
        : Surface2(transform_, isNormalFlipped_), normal(normal_), point(point_) {}

    Plane2::Plane2(const Plane2& other)
        : Surface2(other), normal(other.normal), point(other.point) {}

    bool Plane2::isBounded() const {
        return false;
    }

    Vector2D Plane2::closestPointLocal(const Vector2D& otherPoint) const {
        Vector2D r = otherPoint - point;
        return r - normal.dot(r) * normal + point;
    }

    Vector2D Plane2::closestNormalLocal(const Vector2D& otherPoint) const {
        UNUSED_VARIABLE(otherPoint);
        return normal;
    }

    bool Plane2::intersectsLocal(const Ray2D& ray) const {
        return std::fabs(ray.direction.dot(normal)) > 0;
    }

    SurfaceRayIntersection2 Plane2::closestIntersectionLocal(
        const Ray2D& ray) const {
        SurfaceRayIntersection2 intersection;
        double dDotN = ray.direction.dot(normal);

        // Check if not parallel
        if (std::fabs(dDotN) > 0) {
            double t = normal.dot(point - ray.origin) / dDotN;
            if (t >= 0.0) {
                intersection.isIntersecting = true;
                intersection.distance = t;
                intersection.point = ray.pointAt(t);
                intersection.normal = normal;
            }
        }

        return intersection;
    }

    BoundingBox2D Plane2::boundingBoxLocal() const {
        if (std::fabs(normal.dot(Vector2D(1, 0)) - 1.0) < kEpsilonD) {
            return BoundingBox2D(point - Vector2D(0, kMaxD),
                point + Vector2D(0, kMaxD));
        }
        else if (std::fabs(normal.dot(Vector2D(0, 1)) - 1.0) < kEpsilonD) {
            return BoundingBox2D(point - Vector2D(kMaxD, 0),
                point + Vector2D(kMaxD, 0));
        }
        else {
            return BoundingBox2D(Vector2D(kMaxD, kMaxD), Vector2D(kMaxD, kMaxD));
        }
    }

    Plane2::Builder Plane2::builder() { return Builder(); }

    Plane2::Builder& Plane2::Builder::withNormal(const Vector2D& normal) {
        _normal = normal;
        return *this;
    }

    Plane2::Builder& Plane2::Builder::withPoint(const Vector2D& point) {
        _point = point;
        return *this;
    }

    Plane2 Plane2::Builder::build() const {
        return Plane2(_normal, _point, _transform, _isNormalFlipped);
    }

    Plane2Ptr Plane2::Builder::makeShared() const {
        return std::shared_ptr<Plane2>(
            new Plane2(_normal, _point, _transform, _isNormalFlipped),
            [](Plane2* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_PLANE2_H_