#ifndef INCLUDE_JET_PLANE3_H_
#define INCLUDE_JET_PLANE3_H_

#include "surface3.h"

#include "pch.h"
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief 3-D plane geometry.
    //!
    //! This class represents 3-D plane geometry which extends Surface3 by
    //! overriding surface-related queries.
    //!
    class Plane3 final : public Surface3 {
    public:
        class Builder;

        //! Plane normal.
        Vector3D normal = Vector3D(0, 1, 0);

        //! Point that lies on the plane.
        Vector3D point;

        //! Constructs a plane that crosses (0, 0, 0) with surface normal (0, 1, 0).
        Plane3(
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        //! Constructs a plane that cross \p point with surface normal \p normal.
        Plane3(
            const Vector3D& normal,
            const Vector3D& point,
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        //! Constructs a plane with three points on the surface. The normal will be
        //! set using the counter clockwise direction.
        Plane3(
            const Vector3D& point0,
            const Vector3D& point1,
            const Vector3D& point2,
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        //! Copy constructor.
        Plane3(const Plane3& other);

        //! Returns true if bounding box can be defined.
        bool isBounded() const override;

        //! Returns builder fox Plane3.
        static Builder builder();

    protected:
        Vector3D closestPointLocal(const Vector3D& otherPoint) const override;

        bool intersectsLocal(const Ray3D& ray) const override;

        BoundingBox3D boundingBoxLocal() const override;

        Vector3D closestNormalLocal(const Vector3D& otherPoint) const override;

        SurfaceRayIntersection3 closestIntersectionLocal(
            const Ray3D& ray) const override;
    };

    //! Shared pointer for the Plane3 type.
    typedef std::shared_ptr<Plane3> Plane3Ptr;


    //!
    //! \brief Front-end to create Plane3 objects step by step.
    //!
    class Plane3::Builder final : public SurfaceBuilderBase3<Plane3::Builder> {
    public:
        //! Returns builder with plane normal.
        Builder& withNormal(const Vector3D& normal);

        //! Returns builder with point on the plane.
        Builder& withPoint(const Vector3D& point);

        //! Builds Plane3.
        Plane3 build() const;

        //! Builds shared pointer of Plane3 instance.
        Plane3Ptr makeShared() const;

    private:
        Vector3D _normal{ 0, 1, 0 };
        Vector3D _point{ 0, 0, 0 };
    };

    Plane3::Plane3(const Transform3& transform_, bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_) {}

    Plane3::Plane3(const Vector3D& normal, const Vector3D& point,
        const Transform3& transform_, bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_), normal(normal), point(point) {}

    Plane3::Plane3(const Vector3D& point0, const Vector3D& point1,
        const Vector3D& point2, const Transform3& transform_,
        bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_) {
        normal = (point1 - point0).cross(point2 - point0).normalized();
        point = point0;
    }

    Plane3::Plane3(const Plane3& other)
        : Surface3(other), normal(other.normal), point(other.point) {}

    bool Plane3::isBounded() const {
        return false;
    }

    Vector3D Plane3::closestPointLocal(const Vector3D& otherPoint) const {
        Vector3D r = otherPoint - point;
        return r - normal.dot(r) * normal + point;
    }

    Vector3D Plane3::closestNormalLocal(const Vector3D& otherPoint) const {
        UNUSED_VARIABLE(otherPoint);
        return normal;
    }

    bool Plane3::intersectsLocal(const Ray3D& ray) const {
        return std::fabs(ray.direction.dot(normal)) > 0;
    }

    SurfaceRayIntersection3 Plane3::closestIntersectionLocal(
        const Ray3D& ray) const {
        SurfaceRayIntersection3 intersection;

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

    BoundingBox3D Plane3::boundingBoxLocal() const {
        static const double eps = std::numeric_limits<double>::epsilon();
        static const double dmax = std::numeric_limits<double>::max();

        if (std::fabs(normal.dot(Vector3D(1, 0, 0)) - 1.0) < eps) {
            return BoundingBox3D(point - Vector3D(0, dmax, dmax),
                point + Vector3D(0, dmax, dmax));
        }
        else if (std::fabs(normal.dot(Vector3D(0, 1, 0)) - 1.0) < eps) {
            return BoundingBox3D(point - Vector3D(dmax, 0, dmax),
                point + Vector3D(dmax, 0, dmax));
        }
        else if (std::fabs(normal.dot(Vector3D(0, 0, 1)) - 1.0) < eps) {
            return BoundingBox3D(point - Vector3D(dmax, dmax, 0),
                point + Vector3D(dmax, dmax, 0));
        }
        else {
            return BoundingBox3D(Vector3D(dmax, dmax, dmax),
                Vector3D(dmax, dmax, dmax));
        }
    }

    Plane3::Builder Plane3::builder() { return Builder(); }

    Plane3::Builder& Plane3::Builder::withNormal(const Vector3D& normal) {
        _normal = normal;
        return *this;
    }

    Plane3::Builder& Plane3::Builder::withPoint(const Vector3D& point) {
        _point = point;
        return *this;
    }

    Plane3 Plane3::Builder::build() const {
        return Plane3(_normal, _point, _transform, _isNormalFlipped);
    }

    Plane3Ptr Plane3::Builder::makeShared() const {
        return std::shared_ptr<Plane3>(
            new Plane3(_normal, _point, _transform, _isNormalFlipped),
            [](Plane3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_PLANE3_H_