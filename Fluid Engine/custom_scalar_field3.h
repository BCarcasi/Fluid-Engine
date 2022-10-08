#ifndef INCLUDE_JET_CUSTOM_SCALAR_FIELD3_H_
#define INCLUDE_JET_CUSTOM_SCALAR_FIELD3_H_

#include "scalar_field3.h"

namespace jet {

    // 3-D scalar field with custom field function.
    class CustomScalarField3 final : public ScalarField3 {
    public:
        class Builder;

        //
        // Constructs a field with given function.
        //
        // This constructor creates a field with user-provided function object.
        // To compute derivatives, such as gradient and Laplacian, finite
        // differencing is used. Thus, the differencing resolution also can be
        // provided as the last parameter.
        //
        CustomScalarField3(
            const std::function<double(const Vector3D&)>& customFunction,
            double derivativeResolution = 1e-3);

        //
        // Constructs a field with given field and gradient function.
        //
        // This constructor creates a field with user-provided field and gradient
        // function objects. To compute Laplacian, finite differencing is used.
        // Thus, the differencing resolution also can be provided as the last
        // parameter.
        //
        CustomScalarField3(
            const std::function<double(const Vector3D&)>& customFunction,
            const std::function<Vector3D(const Vector3D&)>& customGradientFunction,
            double derivativeResolution = 1e-3);

        // Constructs a field with given field, gradient, and Laplacian function.
        CustomScalarField3(
            const std::function<double(const Vector3D&)>& customFunction,
            const std::function<Vector3D(const Vector3D&)>& customGradientFunction,
            const std::function<double(const Vector3D&)>& customLaplacianFunction);

        // Returns the sampled value at given position x.
        double sample(const Vector3D& x) const override;

        // Returns the sampler function.
        std::function<double(const Vector3D&)> sampler() const override;

        // Returns the gradient vector at given position x.
        Vector3D gradient(const Vector3D& x) const override;

        // Returns the Laplacian at given position x.
        double laplacian(const Vector3D& x) const override;

        // Returns builder fox CustomScalarField3.
        static Builder builder();

    private:
        std::function<double(const Vector3D&)> _customFunction;
        std::function<Vector3D(const Vector3D&)> _customGradientFunction;
        std::function<double(const Vector3D&)> _customLaplacianFunction;
        double _resolution = 1e-3;
    };

    // Shared pointer type for the CustomScalarField3.
    typedef std::shared_ptr<CustomScalarField3> CustomScalarField3Ptr;


    //
    // Front-end to create CustomScalarField3 objects step by step.
    //
    class CustomScalarField3::Builder final {
    public:
        // Returns builder with field function.
        Builder& withFunction(
            const std::function<double(const Vector3D&)>& func);

        // Returns builder with divergence function.
        Builder& withGradientFunction(
            const std::function<Vector3D(const Vector3D&)>& func);

        // Returns builder with curl function.
        Builder& withLaplacianFunction(
            const std::function<double(const Vector3D&)>& func);

        // Returns builder with derivative resolution.
        Builder& withDerivativeResolution(double resolution);

        // Builds CustomScalarField3.
        CustomScalarField3 build() const;

        // Builds shared pointer of CustomScalarField3 instance.
        CustomScalarField3Ptr makeShared() const;

    private:
        double _resolution = 1e-3;
        std::function<double(const Vector3D&)> _customFunction;
        std::function<Vector3D(const Vector3D&)> _customGradientFunction;
        std::function<double(const Vector3D&)> _customLaplacianFunction;
    };

    CustomScalarField3::CustomScalarField3(
        const std::function<double(const Vector3D&)>& customFunction,
        double derivativeResolution) :
        _customFunction(customFunction),
        _resolution(derivativeResolution) {
    }

    CustomScalarField3::CustomScalarField3(
        const std::function<double(const Vector3D&)>& customFunction,
        const std::function<Vector3D(const Vector3D&)>& customGradientFunction,
        double derivativeResolution) :
        _customFunction(customFunction),
        _customGradientFunction(customGradientFunction),
        _resolution(derivativeResolution) {
    }

    CustomScalarField3::CustomScalarField3(
        const std::function<double(const Vector3D&)>& customFunction,
        const std::function<Vector3D(const Vector3D&)>& customGradientFunction,
        const std::function<double(const Vector3D&)>& customLaplacianFunction) :
        _customFunction(customFunction),
        _customGradientFunction(customGradientFunction),
        _customLaplacianFunction(customLaplacianFunction),
        _resolution(1e-3) {
    }

    double CustomScalarField3::sample(const Vector3D& x) const {
        return _customFunction(x);
    }

    std::function<double(const Vector3D&)> CustomScalarField3::sampler() const {
        return _customFunction;
    }

    Vector3D CustomScalarField3::gradient(const Vector3D& x) const {
        if (_customGradientFunction) {
            return _customGradientFunction(x);
        }
        else {
            double left
                = _customFunction(x - Vector3D(0.5 * _resolution, 0.0, 0.0));
            double right
                = _customFunction(x + Vector3D(0.5 * _resolution, 0.0, 0.0));
            double bottom
                = _customFunction(x - Vector3D(0.0, 0.5 * _resolution, 0.0));
            double top
                = _customFunction(x + Vector3D(0.0, 0.5 * _resolution, 0.0));
            double back
                = _customFunction(x - Vector3D(0.0, 0.0, 0.5 * _resolution));
            double front
                = _customFunction(x + Vector3D(0.0, 0.0, 0.5 * _resolution));

            return Vector3D(
                (right - left) / _resolution,
                (top - bottom) / _resolution,
                (front - back) / _resolution);
        }
    }

    double CustomScalarField3::laplacian(const Vector3D& x) const {
        if (_customLaplacianFunction) {
            return _customLaplacianFunction(x);
        }
        else {
            double center = _customFunction(x);
            double left
                = _customFunction(x - Vector3D(0.5 * _resolution, 0.0, 0.0));
            double right
                = _customFunction(x + Vector3D(0.5 * _resolution, 0.0, 0.0));
            double bottom
                = _customFunction(x - Vector3D(0.0, 0.5 * _resolution, 0.0));
            double top
                = _customFunction(x + Vector3D(0.0, 0.5 * _resolution, 0.0));
            double back
                = _customFunction(x - Vector3D(0.0, 0.0, 0.5 * _resolution));
            double front
                = _customFunction(x + Vector3D(0.0, 0.0, 0.5 * _resolution));

            return (left + right + bottom + top + back + front - 6.0 * center)
                / (_resolution * _resolution);
        }
    }

    CustomScalarField3::Builder CustomScalarField3::builder() {
        return Builder();
    }


    CustomScalarField3::Builder&
        CustomScalarField3::Builder::withFunction(
            const std::function<double(const Vector3D&)>& func) {
        _customFunction = func;
        return *this;
    }

    CustomScalarField3::Builder&
        CustomScalarField3::Builder::withGradientFunction(
            const std::function<Vector3D(const Vector3D&)>& func) {
        _customGradientFunction = func;
        return *this;
    }

    CustomScalarField3::Builder&
        CustomScalarField3::Builder::withLaplacianFunction(
            const std::function<double(const Vector3D&)>& func) {
        _customLaplacianFunction = func;
        return *this;
    }

    CustomScalarField3::Builder&
        CustomScalarField3::Builder::withDerivativeResolution(double resolution) {
        _resolution = resolution;
        return *this;
    }

    CustomScalarField3 CustomScalarField3::Builder::build() const {
        if (_customLaplacianFunction) {
            return CustomScalarField3(
                _customFunction,
                _customGradientFunction,
                _customLaplacianFunction);
        }
        else {
            return CustomScalarField3(
                _customFunction,
                _customGradientFunction,
                _resolution);
        }
    }

    CustomScalarField3Ptr CustomScalarField3::Builder::makeShared() const {
        if (_customLaplacianFunction) {
            return std::shared_ptr<CustomScalarField3>(
                new CustomScalarField3(
                    _customFunction,
                    _customGradientFunction,
                    _customLaplacianFunction),
                [](CustomScalarField3* obj) {
                    delete obj;
                });
        }
        else {
            return std::shared_ptr<CustomScalarField3>(
                new CustomScalarField3(
                    _customFunction,
                    _customGradientFunction,
                    _resolution),
                [](CustomScalarField3* obj) {
                    delete obj;
                });
        }
    }

}  // namespace jet

#endif  // INCLUDE_JET_CUSTOM_SCALAR_FIELD3_H_