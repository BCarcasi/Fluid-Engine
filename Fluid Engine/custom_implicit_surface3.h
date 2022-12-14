#ifndef INCLUDE_JET_CUSTOM_IMPLICIT_SURFACE3_H_
#define INCLUDE_JET_CUSTOM_IMPLICIT_SURFACE3_H_

#include "implicit_surface3.h"
#include "scalar_field3.h"

#include "constants.h"
#include "level_set_utils.h"

namespace jet {

    // Custom 3-D implicit surface using arbitrary function.
    class CustomImplicitSurface3 final : public ImplicitSurface3 {
    public:
        class Builder;

        //
        // Constructs an implicit surface using the given signed-distance function.
        //
        // func Custom SDF function object.
        // domain Bounding box of the SDF if exists.
        // resolution Finite differencing resolution for derivatives.
        // rayMarchingResolution Ray marching resolution for ray tests.
        // maxNumOfIterations Number of iterations for closest point search.
        // transform Local-to-world transform.
        // isNormalFlipped True if normal is flipped.
        //
        CustomImplicitSurface3(const std::function<double(const Vector3D&)>& func,
            const BoundingBox3D& domain = BoundingBox3D(),
            double resolution = 1e-3,
            double rayMarchingResolution = 1e-6,
            unsigned int maxNumOfIterations = 5,
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        // Destructor.
        virtual ~CustomImplicitSurface3();

        // Returns builder for CustomImplicitSurface3.
        static Builder builder();

    private:
        std::function<double(const Vector3D&)> _func;
        BoundingBox3D _domain;
        double _resolution = 1e-3;
        double _rayMarchingResolution = 1e-6;
        unsigned int _maxNumOfIterations = 5;

        Vector3D closestPointLocal(const Vector3D& otherPoint) const override;

        bool intersectsLocal(const Ray3D& ray) const override;

        BoundingBox3D boundingBoxLocal() const override;

        Vector3D closestNormalLocal(const Vector3D& otherPoint) const override;

        double signedDistanceLocal(const Vector3D& otherPoint) const override;

        SurfaceRayIntersection3 closestIntersectionLocal(
            const Ray3D& ray) const override;

        Vector3D gradientLocal(const Vector3D& x) const;
    };

    // Shared pointer type for the CustomImplicitSurface3.
    typedef std::shared_ptr<CustomImplicitSurface3> CustomImplicitSurface3Ptr;

    //
    // Front-end to create CustomImplicitSurface3 objects step by step.
    //
    class CustomImplicitSurface3::Builder final
        : public SurfaceBuilderBase3<CustomImplicitSurface3::Builder> {
    public:
        // Returns builder with custom signed-distance function
        Builder& withSignedDistanceFunction(
            const std::function<double(const Vector3D&)>& func);

        // Returns builder with domain.
        Builder& withDomain(const BoundingBox3D& domain);

        // Returns builder with finite differencing resolution.
        Builder& withResolution(double resolution);

        // Returns builder with ray marching resolution which determines the ray
        // intersection quality.
        Builder& withRayMarchingResolution(double rayMarchingResolution);

        // Returns builder with number of iterations for closest point/normal
        // searches.
        Builder& withMaxNumberOfIterations(unsigned int numIter);

        // Builds CustomImplicitSurface3.
        CustomImplicitSurface3 build() const;

        // Builds shared pointer of CustomImplicitSurface3 instance.
        CustomImplicitSurface3Ptr makeShared() const;

    private:
        std::function<double(const Vector3D&)> _func;
        BoundingBox3D _domain;
        double _resolution = 1e-3;
        double _rayMarchingResolution = 1e-6;
        unsigned int _maxNumOfIterations = 5;
    };

    CustomImplicitSurface3::CustomImplicitSurface3(
        const std::function<double(const Vector3D&)>& func,
        const BoundingBox3D& domain, double resolution,
        double rayMarchingResolution, unsigned int maxNumOfIterations,
        const Transform3& transform, bool isNormalFlipped)
        : ImplicitSurface3(transform, isNormalFlipped),
        _func(func),
        _domain(domain),
        _resolution(resolution),
        _rayMarchingResolution(rayMarchingResolution),
        _maxNumOfIterations(maxNumOfIterations) {}

    CustomImplicitSurface3::~CustomImplicitSurface3() {}

    Vector3D CustomImplicitSurface3::closestPointLocal(
        const Vector3D& otherPoint) const {
        Vector3D pt = clamp(otherPoint, _domain.lowerCorner, _domain.upperCorner);
        for (unsigned int iter = 0; iter < _maxNumOfIterations; ++iter) {
            double sdf = signedDistanceLocal(pt);
            if (std::fabs(sdf) < kEpsilonD) {
                break;
            }
            Vector3D g = gradientLocal(pt);
            pt = pt - sdf * g;
        }
        return pt;
    }

    bool CustomImplicitSurface3::intersectsLocal(const Ray3D& ray) const {
        BoundingBoxRayIntersection3D intersection =
            _domain.closestIntersection(ray);

        if (intersection.isIntersecting) {
            double tStart, tEnd;
            if (intersection.tFar == kMaxD) {
                tStart = 0.0;
                tEnd = intersection.tNear;
            }
            else {
                tStart = intersection.tNear;
                tEnd = intersection.tFar;
            }

            double t = tStart;
            Vector3D pt = ray.pointAt(t);
            double prevPhi = _func(pt);
            while (t <= tEnd) {
                pt = ray.pointAt(t);
                const double newPhi = _func(pt);
                const double newPhiAbs = std::fabs(newPhi);

                if (newPhi * prevPhi < 0.0) {
                    return true;
                }

                t += std::max(newPhiAbs, _rayMarchingResolution);
                prevPhi = newPhi;
            }
        }

        return false;
    }

    BoundingBox3D CustomImplicitSurface3::boundingBoxLocal() const {
        return _domain;
    }

    double CustomImplicitSurface3::signedDistanceLocal(
        const Vector3D& otherPoint) const {
        if (_func) {
            return _func(otherPoint);
        }
        else {
            return kMaxD;
        }
    }

    Vector3D CustomImplicitSurface3::closestNormalLocal(
        const Vector3D& otherPoint) const {
        Vector3D pt = closestPointLocal(otherPoint);
        Vector3D g = gradientLocal(pt);
        if (g.lengthSquared() > 0.0) {
            return g.normalized();
        }
        else {
            return g;
        }
    }

    SurfaceRayIntersection3 CustomImplicitSurface3::closestIntersectionLocal(
        const Ray3D& ray) const {
        SurfaceRayIntersection3 result;

        BoundingBoxRayIntersection3D intersection =
            _domain.closestIntersection(ray);

        if (intersection.isIntersecting) {
            double tStart, tEnd;
            if (intersection.tFar == kMaxD) {
                tStart = 0.0;
                tEnd = intersection.tNear;
            }
            else {
                tStart = intersection.tNear;
                tEnd = intersection.tFar;
            }

            double t = tStart;
            double tPrev = t;
            Vector3D pt = ray.pointAt(t);
            double prevPhi = _func(pt);

            while (t <= tEnd) {
                pt = ray.pointAt(t);
                const double newPhi = _func(pt);
                const double newPhiAbs = std::fabs(newPhi);

                if (newPhi * prevPhi < 0.0) {
                    const double frac = prevPhi / (prevPhi - newPhi);
                    const double tSub = tPrev + _rayMarchingResolution * frac;

                    result.isIntersecting = true;
                    result.distance = tSub;
                    result.point = ray.pointAt(tSub);
                    result.normal = gradientLocal(result.point);
                    if (result.normal.length() > 0.0) {
                        result.normal.normalize();
                    }

                    return result;
                }

                tPrev = t;
                t += std::max(newPhiAbs, _rayMarchingResolution);
                prevPhi = newPhi;
            }
        }

        return result;
    }

    Vector3D CustomImplicitSurface3::gradientLocal(const Vector3D& x) const {
        double left = _func(x - Vector3D(0.5 * _resolution, 0.0, 0.0));
        double right = _func(x + Vector3D(0.5 * _resolution, 0.0, 0.0));
        double bottom = _func(x - Vector3D(0.0, 0.5 * _resolution, 0.0));
        double top = _func(x + Vector3D(0.0, 0.5 * _resolution, 0.0));
        double back = _func(x - Vector3D(0.0, 0.0, 0.5 * _resolution));
        double front = _func(x + Vector3D(0.0, 0.0, 0.5 * _resolution));

        return Vector3D((right - left) / _resolution, (top - bottom) / _resolution,
            (front - back) / _resolution);
    }

    CustomImplicitSurface3::Builder CustomImplicitSurface3::builder() {
        return Builder();
    }

    CustomImplicitSurface3::Builder&
        CustomImplicitSurface3::Builder::withSignedDistanceFunction(
            const std::function<double(const Vector3D&)>& func) {
        _func = func;
        return *this;
    }

    CustomImplicitSurface3::Builder& CustomImplicitSurface3::Builder::withDomain(
        const BoundingBox3D& domain) {
        _domain = domain;
        return *this;
    }

    CustomImplicitSurface3::Builder&
        CustomImplicitSurface3::Builder::withResolution(double resolution) {
        _resolution = resolution;
        return *this;
    }

    CustomImplicitSurface3::Builder&
        CustomImplicitSurface3::Builder::withRayMarchingResolution(double resolution) {
        _rayMarchingResolution = resolution;
        return *this;
    }

    CustomImplicitSurface3::Builder&
        CustomImplicitSurface3::Builder::withMaxNumberOfIterations(
            unsigned int numIter) {
        _maxNumOfIterations = numIter;
        return *this;
    }

    CustomImplicitSurface3 CustomImplicitSurface3::Builder::build() const {
        return CustomImplicitSurface3(_func, _domain, _resolution,
            _rayMarchingResolution, _maxNumOfIterations,
            _transform, _isNormalFlipped);
    }

    CustomImplicitSurface3Ptr CustomImplicitSurface3::Builder::makeShared() const {
        return std::shared_ptr<CustomImplicitSurface3>(
            new CustomImplicitSurface3(_func, _domain, _resolution,
                _rayMarchingResolution, _maxNumOfIterations,
                _transform, _isNormalFlipped),
            [](CustomImplicitSurface3* obj) { delete obj; });
    }
}  // namespace jet

#endif  // INCLUDE_JET_CUSTOM_IMPLICIT_SURFACE3_H_