#ifndef INCLUDE_JET_IMPLICIT_SURFACE3_H_
#define INCLUDE_JET_IMPLICIT_SURFACE3_H_

#include "surface3.h"
#include "level_set_utils.h"

namespace jet {

    // Abstract base class for 3-D implicit surface.
    class ImplicitSurface3 : public Surface3 {
    public:
        // Default constructor.
        ImplicitSurface3(
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        // Copy constructor.
        ImplicitSurface3(const ImplicitSurface3& other);

        // Default destructor.
        virtual ~ImplicitSurface3();

        // Returns signed distance from the given point otherPoint.
        double signedDistance(const Vector3D& otherPoint) const;

    protected:
        // Returns signed distance from the given point otherPoint in local
        // space.
        virtual double signedDistanceLocal(const Vector3D& otherPoint) const = 0;

    private:
        double closestDistanceLocal(const Vector3D& otherPoint) const override;

        bool isInsideLocal(const Vector3D& otherPoint) const override;
    };

    // Shared pointer type for the ImplicitSurface3.
    typedef std::shared_ptr<ImplicitSurface3> ImplicitSurface3Ptr;

    ImplicitSurface3::ImplicitSurface3(
        const Transform3& transform_, bool isNormalFlipped_)
        : Surface3(transform_, isNormalFlipped_) {
    }

    ImplicitSurface3::ImplicitSurface3(const ImplicitSurface3& other) :
        Surface3(other) {
    }

    ImplicitSurface3::~ImplicitSurface3() {
    }

    double ImplicitSurface3::signedDistance(const Vector3D& otherPoint) const {
        double sd = signedDistanceLocal(transform.toLocal(otherPoint));
        return (isNormalFlipped) ? -sd : sd;
    }

    double ImplicitSurface3::closestDistanceLocal(
        const Vector3D& otherPoint) const {
        return std::fabs(signedDistanceLocal(otherPoint));
    }

    bool ImplicitSurface3::isInsideLocal(const Vector3D& otherPoint) const {
        return isInsideSdf(signedDistanceLocal(otherPoint));
    }

}  // namespace jet

#endif  // INCLUDE_JET_IMPLICIT_SURFACE3_H_