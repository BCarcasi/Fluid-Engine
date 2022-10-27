#ifndef INCLUDE_JET_BOX2_H_
#define INCLUDE_JET_BOX2_H_

#include "surface2.h"
#include "bounding_box2.h"

#include "pch.h"

#include "plane2.h"

namespace jet {

    //!
    //! \brief 2-D box geometry.
    //!
    //! This class represents 2-D box geometry which extends Surface2 by overriding
    //! surface-related queries. This box implementation is an axis-aligned box
    //! that wraps lower-level primitive type, BoundingBox2D.
    //!
    class Box2 final : public Surface2 {
    public:
        class Builder;

        //! Bounding box of this box.
        BoundingBox2D bound = BoundingBox2D(Vector2D(), Vector2D(1.0, 1.0));

        //! Constructs (0, 0) x (1, 1) box.
        Box2(
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        //! Constructs a box with given \p lowerCorner and \p upperCorner.
        Box2(
            const Vector2D& lowerCorner,
            const Vector2D& upperCorner,
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        //! Constructs a box with BoundingBox2D instance.
        Box2(
            const BoundingBox2D& boundingBox,
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        //! Copy constructor.
        Box2(const Box2& other);

        //! Returns builder fox Box2.
        static Builder builder();

    protected:
        // Surface2 implementations

        Vector2D closestPointLocal(const Vector2D& otherPoint) const override;

        bool intersectsLocal(const Ray2D& ray) const override;

        BoundingBox2D boundingBoxLocal() const override;

        Vector2D closestNormalLocal(const Vector2D& otherPoint) const override;

        SurfaceRayIntersection2 closestIntersectionLocal(
            const Ray2D& ray) const override;
    };

    //! Shared pointer type for the Box2.
    typedef std::shared_ptr<Box2> Box2Ptr;


    //!
    //! \brief Front-end to create Box2 objects step by step.
    //!
    class Box2::Builder final : public SurfaceBuilderBase2<Box2::Builder> {
    public:
        //! Returns builder with lower corner set.
        Builder& withLowerCorner(const Vector2D& pt);

        //! Returns builder with upper corner set.
        Builder& withUpperCorner(const Vector2D& pt);

        //! Returns builder with bounding box.
        Builder& withBoundingBox(const BoundingBox2D& bbox);

        //! Builds Box2.
        Box2 build() const;

        //! Builds shared pointer of Box2 instance.
        Box2Ptr makeShared() const;

    private:
        Vector2D _lowerCorner{ 0, 0 };
        Vector2D _upperCorner{ 1, 1 };
    };

    Box2::Box2(const Transform2& transform, bool isNormalFlipped)
        : Surface2(transform, isNormalFlipped) {}

    Box2::Box2(const Vector2D& lowerCorner, const Vector2D& upperCorner,
        const Transform2& transform, bool isNormalFlipped)
        : Box2(BoundingBox2D(lowerCorner, upperCorner), transform,
            isNormalFlipped) {}

    Box2::Box2(const BoundingBox2D& boundingBox, const Transform2& transform,
        bool isNormalFlipped)
        : Surface2(transform, isNormalFlipped), bound(boundingBox) {}

    Box2::Box2(const Box2& other) : Surface2(other), bound(other.bound) {}

    Vector2D Box2::closestPointLocal(const Vector2D& otherPoint) const {
        if (bound.contains(otherPoint)) {
            Plane2 planes[4] = { Plane2(Vector2D(1, 0), bound.upperCorner),
                                Plane2(Vector2D(0, 1), bound.upperCorner),
                                Plane2(Vector2D(-1, 0), bound.lowerCorner),
                                Plane2(Vector2D(0, -1), bound.lowerCorner) };

            Vector2D result = planes[0].closestPoint(otherPoint);
            double distanceSquared = result.distanceSquaredTo(otherPoint);

            for (int i = 1; i < 4; ++i) {
                Vector2D localResult = planes[i].closestPoint(otherPoint);
                double localDistanceSquared =
                    localResult.distanceSquaredTo(otherPoint);

                if (localDistanceSquared < distanceSquared) {
                    result = localResult;
                    distanceSquared = localDistanceSquared;
                }
            }

            return result;
        }
        else {
            return clamp(otherPoint, bound.lowerCorner, bound.upperCorner);
        }
    }

    Vector2D Box2::closestNormalLocal(const Vector2D& otherPoint) const {
        Plane2 planes[4] = { Plane2(Vector2D(1, 0), bound.upperCorner),
                            Plane2(Vector2D(0, 1), bound.upperCorner),
                            Plane2(Vector2D(-1, 0), bound.lowerCorner),
                            Plane2(Vector2D(0, -1), bound.lowerCorner) };

        if (bound.contains(otherPoint)) {
            Vector2D closestNormal = planes[0].normal;
            Vector2D closestPoint = planes[0].closestPoint(otherPoint);
            double minDistanceSquared = (closestPoint - otherPoint).lengthSquared();

            for (int i = 1; i < 4; ++i) {
                Vector2D localClosestPoint = planes[i].closestPoint(otherPoint);
                double localDistanceSquared =
                    (localClosestPoint - otherPoint).lengthSquared();

                if (localDistanceSquared < minDistanceSquared) {
                    closestNormal = planes[i].normal;
                    minDistanceSquared = localDistanceSquared;
                }
            }

            return closestNormal;
        }
        else {
            Vector2D closestPoint =
                clamp(otherPoint, bound.lowerCorner, bound.upperCorner);
            Vector2D closestPointToInputPoint = otherPoint - closestPoint;
            Vector2D closestNormal = planes[0].normal;
            double maxCosineAngle = closestNormal.dot(closestPointToInputPoint);

            for (int i = 1; i < 4; ++i) {
                double cosineAngle = planes[i].normal.dot(closestPointToInputPoint);

                if (cosineAngle > maxCosineAngle) {
                    closestNormal = planes[i].normal;
                    maxCosineAngle = cosineAngle;
                }
            }

            return closestNormal;
        }
    }

    bool Box2::intersectsLocal(const Ray2D& ray) const {
        return bound.intersects(ray);
    }

    SurfaceRayIntersection2 Box2::closestIntersectionLocal(const Ray2D& ray) const {
        SurfaceRayIntersection2 intersection;
        BoundingBoxRayIntersection2D bbRayIntersection =
            bound.closestIntersection(ray);
        intersection.isIntersecting = bbRayIntersection.isIntersecting;
        if (intersection.isIntersecting) {
            intersection.distance = bbRayIntersection.tNear;
            intersection.point = ray.pointAt(bbRayIntersection.tNear);
            intersection.normal = closestNormalLocal(intersection.point);
        }
        return intersection;
    }

    BoundingBox2D Box2::boundingBoxLocal() const { return bound; }

    Box2::Builder Box2::builder() { return Builder(); }

    Box2::Builder& Box2::Builder::withLowerCorner(const Vector2D& pt) {
        _lowerCorner = pt;
        return *this;
    }

    Box2::Builder& Box2::Builder::withUpperCorner(const Vector2D& pt) {
        _upperCorner = pt;
        return *this;
    }

    Box2::Builder& Box2::Builder::withBoundingBox(const BoundingBox2D& bbox) {
        _lowerCorner = bbox.lowerCorner;
        _upperCorner = bbox.upperCorner;
        return *this;
    }

    Box2 Box2::Builder::build() const {
        return Box2(_lowerCorner, _upperCorner, _transform, _isNormalFlipped);
    }

    Box2Ptr Box2::Builder::makeShared() const {
        return std::shared_ptr<Box2>(
            new Box2(_lowerCorner, _upperCorner, _transform, _isNormalFlipped),
            [](Box2* obj) { delete obj; });
    }

}  // namespace jet


#endif  // INCLUDE_JET_BOX2_H_