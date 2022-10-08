#ifndef INCLUDE_JET_SURFACE_TO_IMPLICIT2_H_
#define INCLUDE_JET_SURFACE_TO_IMPLICIT2_H_

#include "implicit_surface2.h"
#include <memory>

namespace jet {

    //
    // 2-D implicit surface wrapper for generic Surface2 instance.
    //
    // This class represents 2-D implicit surface that converts Surface2 instance
    // to an ImplicitSurface2 object. The conversion is made by evaluating closest
    // point and normal from a given point for the given (explicit) surface. Thus,
    // this conversion won't work for every single surfaces. Use this class only
    // for the basic primitives such as Sphere2 or Box2.
    //
    class SurfaceToImplicit2 final : public ImplicitSurface2 {
    public:
        class Builder;

        // Constructs an instance with generic Surface2 instance.
        SurfaceToImplicit2(
            const Surface2Ptr& surface,
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        // Copy constructor.
        SurfaceToImplicit2(const SurfaceToImplicit2& other);

        // Updates internal spatial query engine.
        void updateQueryEngine() override;

        // Returns true if bounding box can be defined.
        bool isBounded() const override;

        // Returns true if the surface is a valid geometry.
        bool isValidGeometry() const override;

        // Returns the raw surface instance.
        Surface2Ptr surface() const;

        // Returns builder fox SurfaceToImplicit2.
        static Builder builder();

    protected:
        Vector2D closestPointLocal(const Vector2D& otherPoint) const override;

        double closestDistanceLocal(const Vector2D& otherPoint) const override;

        bool intersectsLocal(const Ray2D& ray) const override;

        BoundingBox2D boundingBoxLocal() const override;

        Vector2D closestNormalLocal(
            const Vector2D& otherPoint) const override;

        double signedDistanceLocal(const Vector2D& otherPoint) const override;

        SurfaceRayIntersection2 closestIntersectionLocal(
            const Ray2D& ray) const override;

        bool isInsideLocal(const Vector2D& otherPoint) const override;

    private:
        Surface2Ptr _surface;
    };

    // Shared pointer for the SurfaceToImplicit2 type.
    typedef std::shared_ptr<SurfaceToImplicit2> SurfaceToImplicit2Ptr;


    //
    // Front-end to create SurfaceToImplicit2 objects step by step.
    //
    class SurfaceToImplicit2::Builder final
        : public SurfaceBuilderBase2<SurfaceToImplicit2::Builder> {
    public:
        // Returns builder with surface.
        Builder& withSurface(const Surface2Ptr& surface);

        // Builds SurfaceToImplicit2.
        SurfaceToImplicit2 build() const;

        // Builds shared pointer of SurfaceToImplicit2 instance.
        SurfaceToImplicit2Ptr makeShared() const;

    private:
        Surface2Ptr _surface;
    };

    SurfaceToImplicit2::SurfaceToImplicit2(const Surface2Ptr& surface,
        const Transform2& transform,
        bool isNormalFlipped)
        : ImplicitSurface2(transform, isNormalFlipped), _surface(surface) {}

    SurfaceToImplicit2::SurfaceToImplicit2(const SurfaceToImplicit2& other)
        : ImplicitSurface2(other), _surface(other._surface) {}

    bool SurfaceToImplicit2::isBounded() const { return _surface->isBounded(); }

    void SurfaceToImplicit2::updateQueryEngine() { _surface->updateQueryEngine(); }

    bool SurfaceToImplicit2::isValidGeometry() const {
        return _surface->isValidGeometry();
    }

    Surface2Ptr SurfaceToImplicit2::surface() const { return _surface; }

    SurfaceToImplicit2::Builder SurfaceToImplicit2::builder() { return Builder(); }

    Vector2D SurfaceToImplicit2::closestPointLocal(
        const Vector2D& otherPoint) const {
        return _surface->closestPoint(otherPoint);
    }

    Vector2D SurfaceToImplicit2::closestNormalLocal(
        const Vector2D& otherPoint) const {
        return _surface->closestNormal(otherPoint);
    }

    double SurfaceToImplicit2::closestDistanceLocal(
        const Vector2D& otherPoint) const {
        return _surface->closestDistance(otherPoint);
    }

    bool SurfaceToImplicit2::intersectsLocal(const Ray2D& ray) const {
        return _surface->intersects(ray);
    }

    SurfaceRayIntersection2 SurfaceToImplicit2::closestIntersectionLocal(
        const Ray2D& ray) const {
        return _surface->closestIntersection(ray);
    }

    BoundingBox2D SurfaceToImplicit2::boundingBoxLocal() const {
        return _surface->boundingBox();
    }

    bool SurfaceToImplicit2::isInsideLocal(const Vector2D& otherPoint) const {
        return _surface->isInside(otherPoint);
    }

    double SurfaceToImplicit2::signedDistanceLocal(
        const Vector2D& otherPoint) const {
        Vector2D x = _surface->closestPoint(otherPoint);
        bool inside = _surface->isInside(otherPoint);
        return (inside) ? -x.distanceTo(otherPoint) : x.distanceTo(otherPoint);
    }

    SurfaceToImplicit2::Builder& SurfaceToImplicit2::Builder::withSurface(
        const Surface2Ptr& surface) {
        _surface = surface;
        return *this;
    }

    SurfaceToImplicit2 SurfaceToImplicit2::Builder::build() const {
        return SurfaceToImplicit2(_surface, _transform, _isNormalFlipped);
    }

    SurfaceToImplicit2Ptr SurfaceToImplicit2::Builder::makeShared() const {
        return std::shared_ptr<SurfaceToImplicit2>(
            new SurfaceToImplicit2(_surface, _transform, _isNormalFlipped),
            [](SurfaceToImplicit2* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_SURFACE_TO_IMPLICIT2_H_