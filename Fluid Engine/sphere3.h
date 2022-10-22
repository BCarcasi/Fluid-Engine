#ifndef INCLUDE_JET_SPHERE3_H_
#define INCLUDE_JET_SPHERE3_H_

#include "surface3.h"
#include "bounding_box3.h"

#include "pch.h"

namespace jet {

    //!
    //! \brief 3-D sphere geometry.
    //!
    //! This class represents 3-D sphere geometry which extends Surface3 by
    //! overriding surface-related queries.
    //!
    class Sphere3 final : public Surface3 {
    public:
        class Builder;

        //! Center of the sphere.
        Vector3D center;

        //! Radius of the sphere.
        double radius = 1.0;

        //! Constructs a sphere with center at (0, 0, 0) and radius of 1.
        Sphere3(
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        //! Constructs a sphere with \p center and \p radius.
        Sphere3(
            const Vector3D& center,
            double radius,
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        //! Copy constructor.
        Sphere3(const Sphere3& other);

        //! Returns builder fox Sphere3.
        static Builder builder();

    private:
        Vector3D closestPointLocal(const Vector3D& otherPoint) const override;

        double closestDistanceLocal(const Vector3D& otherPoint) const override;

        bool intersectsLocal(const Ray3D& ray) const override;

        BoundingBox3D boundingBoxLocal() const override;

        Vector3D closestNormalLocal(const Vector3D& otherPoint) const override;

        SurfaceRayIntersection3 closestIntersectionLocal(
            const Ray3D& ray) const override;
    };

    //! Shared pointer for the Sphere3 type.
    typedef std::shared_ptr<Sphere3> Sphere3Ptr;

    //!
    //! \brief Front-end to create Sphere3 objects step by step.
    //!
    class Sphere3::Builder final : public SurfaceBuilderBase3<Sphere3::Builder> {
    public:
        //! Returns builder with sphere center.
        Builder& withCenter(const Vector3D& center);

        //! Returns builder with sphere radius.
        Builder& withRadius(double radius);

        //! Builds Sphere3.
        Sphere3 build() const;

        //! Builds shared pointer of Sphere3 instance.
        Sphere3Ptr makeShared() const;

    private:
        Vector3D _center{ 0, 0, 0 };
        double _radius = 0.0;
    };

    Sphere3::Sphere3(const Transform3& transform_, bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_) {}

    Sphere3::Sphere3(const Vector3D& center_, double radius_,
        const Transform3& transform_, bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_),
        center(center_),
        radius(radius_) {}

    Sphere3::Sphere3(const Sphere3& other)
        : Surface3(other), center(other.center), radius(other.radius) {}

    Vector3D Sphere3::closestPointLocal(const Vector3D& otherPoint) const {
        return radius * closestNormalLocal(otherPoint) + center;
    }

    double Sphere3::closestDistanceLocal(const Vector3D& otherPoint) const {
        return std::fabs(center.distanceTo(otherPoint) - radius);
    }

    Vector3D Sphere3::closestNormalLocal(const Vector3D& otherPoint) const {
        if (center.isSimilar(otherPoint)) {
            return Vector3D(1, 0, 0);
        }
        else {
            return (otherPoint - center).normalized();
        }
    }

    bool Sphere3::intersectsLocal(const Ray3D& ray) const {
        Vector3D r = ray.origin - center;
        double b = ray.direction.dot(r);
        double c = r.lengthSquared() - square(radius);
        double d = b * b - c;

        if (d > 0.) {
            d = std::sqrt(d);
            double tMin = -b - d;
            double tMax = -b + d;

            if (tMin < 0.0) {
                tMin = tMax;
            }

            if (tMin >= 0.0) {
                return true;
            }
        }

        return false;
    }

    SurfaceRayIntersection3 Sphere3::closestIntersectionLocal(
        const Ray3D& ray) const {
        SurfaceRayIntersection3 intersection;
        Vector3D r = ray.origin - center;
        double b = ray.direction.dot(r);
        double c = r.lengthSquared() - square(radius);
        double d = b * b - c;

        if (d > 0.) {
            d = std::sqrt(d);
            double tMin = -b - d;
            double tMax = -b + d;

            if (tMin < 0.0) {
                tMin = tMax;
            }

            if (tMin < 0.0) {
                intersection.isIntersecting = false;
            }
            else {
                intersection.isIntersecting = true;
                intersection.distance = tMin;
                intersection.point = ray.origin + tMin * ray.direction;
                intersection.normal = (intersection.point - center).normalized();
            }
        }
        else {
            intersection.isIntersecting = false;
        }

        return intersection;
    }

    BoundingBox3D Sphere3::boundingBoxLocal() const {
        Vector3D r(radius, radius, radius);
        return BoundingBox3D(center - r, center + r);
    }

    Sphere3::Builder Sphere3::builder() { return Builder(); }

    Sphere3::Builder& Sphere3::Builder::withCenter(const Vector3D& center) {
        _center = center;
        return *this;
    }

    Sphere3::Builder& Sphere3::Builder::withRadius(double radius) {
        _radius = radius;
        return *this;
    }

    Sphere3 Sphere3::Builder::build() const {
        return Sphere3(_center, _radius, _transform, _isNormalFlipped);
    }

    Sphere3Ptr Sphere3::Builder::makeShared() const {
        return std::shared_ptr<Sphere3>(
            new Sphere3(_center, _radius, _transform, _isNormalFlipped),
            [](Sphere3* obj) { delete obj; });
    }

}  // namespace jet


#endif  // INCLUDE_JET_SPHERE3_H_