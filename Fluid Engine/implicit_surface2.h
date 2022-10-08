#ifndef INCLUDE_JET_IMPLICIT_SURFACE2_H_
#define INCLUDE_JET_IMPLICIT_SURFACE2_H_

#include "surface2.h"
#include "level_set_utils.h"

namespace jet {

    // Abstract base class for 2-D implicit surface.
    class ImplicitSurface2 : public Surface2 {
    public:
        // Default constructor.
        ImplicitSurface2(
            const Transform2& transform = Transform2(),
            bool isNormalFlipped = false);

        // Copy constructor.
        ImplicitSurface2(const ImplicitSurface2& other);

        // Default destructor.
        virtual ~ImplicitSurface2();

        // Returns signed distance from the given point otherPoint.
        double signedDistance(const Vector2D& otherPoint) const;

    protected:
        // Returns signed distance from the given point otherPoint in local
        // space.
        virtual double signedDistanceLocal(const Vector2D& otherPoint) const = 0;

    private:
        double closestDistanceLocal(const Vector2D& otherPoint) const override;

        bool isInsideLocal(const Vector2D& otherPoint) const override;
    };

    // Shared pointer type for the ImplicitSurface2.
    typedef std::shared_ptr<ImplicitSurface2> ImplicitSurface2Ptr;

    ImplicitSurface2::ImplicitSurface2(
        const Transform2& transform, bool isNormalFlipped)
        : Surface2(transform, isNormalFlipped) {
    }

    ImplicitSurface2::ImplicitSurface2(const ImplicitSurface2& other) :
        Surface2(other) {
    }

    ImplicitSurface2::~ImplicitSurface2() {
    }

    double ImplicitSurface2::signedDistance(const Vector2D& otherPoint) const {
        double sd = signedDistanceLocal(transform.toLocal(otherPoint));
        return (isNormalFlipped) ? -sd : sd;
    }

    double ImplicitSurface2::closestDistanceLocal(
        const Vector2D& otherPoint) const {
        return std::fabs(signedDistanceLocal(otherPoint));
    }

    bool ImplicitSurface2::isInsideLocal(const Vector2D& otherPoint) const {
        return isInsideSdf(signedDistanceLocal(otherPoint));
    }

}  // namespace jet

#endif  // INCLUDE_JET_IMPLICIT_SURFACE2_H_