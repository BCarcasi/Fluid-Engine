#ifndef INCLUDE_JET_SURFACE2_H_
#define INCLUDE_JET_SURFACE2_H_

#include "bounding_box2.h"
#include "constants.h"
#include "ray2.h"
#include "transform2.h"
#include <memory>

#include "pch.h"

#include <algorithm>


namespace jet {

    // Struct that represents ray-surface intersection point.
    struct SurfaceRayIntersection2 {
        bool isIntersecting = false;
        double distance = kMaxD;
        Vector2D point;
        Vector2D normal;
    };

    // Abstract base class for 2-D surface.
    class Surface2 {
    public:
        // Local-to-world transform.
        Transform2 transform;

        // Flips normal.
        bool isNormalFlipped = false;

        // Constructs a surface with normal direction.
        Surface2(const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        // Copy constructor.
        Surface2(const Surface2& other);

        // Default destructor.
        virtual ~Surface2();

        // Returns the closest point from the given point otherPoint to the
        // surface.
        Vector2D closestPoint(const Vector2D& otherPoint) const;

        // Returns the bounding box of this surface object.
        BoundingBox2D boundingBox() const;

        // Returns true if the given ray intersects with this surface object.
        bool intersects(const Ray2D& ray) const;

        // Returns the closest distance from the given point otherPoint to the
        // point on the surface.
        double closestDistance(const Vector2D& otherPoint) const;

        // Returns the closest intersection point for given ray.
        SurfaceRayIntersection2 closestIntersection(const Ray2D& ray) const;

        // Returns the normal to the closest point on the surface from the given
        // point otherPoint.
        Vector2D closestNormal(const Vector2D& otherPoint) const;

        // Updates internal spatial query engine.
        virtual void updateQueryEngine();

        // Returns true if bounding box can be defined.
        virtual bool isBounded() const;

        // Returns true if the surface is a valid geometry.
        virtual bool isValidGeometry() const;

        // Returns true if otherPoint is inside the volume defined by the
        // surface.
        bool isInside(const Vector2D& otherPoint) const;

    protected:
        // Returns the closest point from the given point otherPoint to the
        // surface in local frame.
        virtual Vector2D closestPointLocal(const Vector2D& otherPoint) const = 0;

        // Returns the bounding box of this surface object in local frame.
        virtual BoundingBox2D boundingBoxLocal() const = 0;

        // Returns the closest intersection point for given ray in local frame.
        virtual SurfaceRayIntersection2 closestIntersectionLocal(
            const Ray2D& ray) const = 0;

        // Returns the normal to the closest point on the surface from the given
        // point otherPoint in local frame.
        virtual Vector2D closestNormalLocal(const Vector2D& otherPoint) const = 0;

        // Returns true if the given ray intersects with this surface object
        // in local frame.
        virtual bool intersectsLocal(const Ray2D& ray) const;

        // Returns the closest distance from the given point otherPoint to the
        // point on the surface in local frame.
        virtual double closestDistanceLocal(const Vector2D& otherPoint) const;

        // Returns true if otherPoint is inside by given depth the volume
        // defined by the surface in local frame.
        virtual bool isInsideLocal(const Vector2D& otherPoint) const;
    };

    // Shared pointer for the Surface2 type.
    typedef std::shared_ptr<Surface2> Surface2Ptr;

    //
    // Base class for 2-D surface builder.
    //
    template <typename DerivedBuilder>
    class SurfaceBuilderBase2 {
    public:
        // Returns builder with flipped normal flag.
        DerivedBuilder& withIsNormalFlipped(bool isNormalFlipped);

        // Returns builder with translation.
        DerivedBuilder& withTranslation(const Vector2D& translation);

        // Returns builder with orientation.
        DerivedBuilder& withOrientation(double orientation);

        // Returns builder with transform.
        DerivedBuilder& withTransform(const Transform2& transform);

    protected:
        bool _isNormalFlipped = false;
        Transform2 _transform;
    };

    template <typename T>
    T& SurfaceBuilderBase2<T>::withIsNormalFlipped(bool isNormalFlipped) {
        _isNormalFlipped = isNormalFlipped;
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& SurfaceBuilderBase2<T>::withTranslation(const Vector2D& translation) {
        _transform.setTranslation(translation);
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& SurfaceBuilderBase2<T>::withOrientation(double orientation) {
        _transform.setOrientation(orientation);
        return static_cast<T&>(*this);
    }

    template <typename T>
    T& SurfaceBuilderBase2<T>::withTransform(const Transform2& transform) {
        _transform = transform;
        return static_cast<T&>(*this);
    }

    Surface2::Surface2(const Transform2& transform_, bool isNormalFlipped_)
        : transform(transform_), isNormalFlipped(isNormalFlipped_) {}

    Surface2::Surface2(const Surface2& other)
        : transform(other.transform), isNormalFlipped(other.isNormalFlipped) {}

    Surface2::~Surface2() {}

    Vector2D Surface2::closestPoint(const Vector2D& otherPoint) const {
        return transform.toWorld(closestPointLocal(transform.toLocal(otherPoint)));
    }

    BoundingBox2D Surface2::boundingBox() const {
        return transform.toWorld(boundingBoxLocal());
    }

    bool Surface2::intersects(const Ray2D& ray) const {
        return intersectsLocal(transform.toLocal(ray));
    }

    double Surface2::closestDistance(const Vector2D& otherPoint) const {
        return closestDistanceLocal(transform.toLocal(otherPoint));
    }

    SurfaceRayIntersection2 Surface2::closestIntersection(const Ray2D& ray) const {
        auto result = closestIntersectionLocal(transform.toLocal(ray));
        result.point = transform.toWorld(result.point);
        result.normal = transform.toWorldDirection(result.normal);
        result.normal *= (isNormalFlipped) ? -1.0 : 1.0;
        return result;
    }

    Vector2D Surface2::closestNormal(const Vector2D& otherPoint) const {
        auto result = transform.toWorldDirection(
            closestNormalLocal(transform.toLocal(otherPoint)));
        result *= (isNormalFlipped) ? -1.0 : 1.0;
        return result;
    }

    void Surface2::updateQueryEngine() {
        // Do nothing
    }

    bool Surface2::isBounded() const { return true; }

    bool Surface2::isValidGeometry() const { return true; }

    bool Surface2::isInside(const Vector2D& otherPoint) const {
        return isNormalFlipped == !isInsideLocal(transform.toLocal(otherPoint));
    }

    bool Surface2::intersectsLocal(const Ray2D& rayLocal) const {
        auto result = closestIntersectionLocal(rayLocal);
        return result.isIntersecting;
    }

    double Surface2::closestDistanceLocal(const Vector2D& otherPointLocal) const {
        return otherPointLocal.distanceTo(closestPointLocal(otherPointLocal));
    }

    bool Surface2::isInsideLocal(const Vector2D& otherPointLocal) const {
        Vector2D cpLocal = closestPointLocal(otherPointLocal);
        Vector2D normalLocal = closestNormalLocal(otherPointLocal);
        return (otherPointLocal - cpLocal).dot(normalLocal) < 0.0;
    }

}  // namespace jet

#endif  // INCLUDE_JET_SURFACE2_H_