#pragma once

#include "field3.h"
#include "vector3.h"
#include <functional>
#include <memory>

namespace jet {

    // Abstract base class for 3-D scalar field.
    class ScalarField3 : public Field3 {
    public:
        // Default constructor.
        ScalarField3();

        // Default destructor.
        virtual ~ScalarField3();

        // Returns sampled value at given position x.
        virtual double sample(const Vector3D& x) const = 0;

        // Returns gradient vector at given position x.
        virtual Vector3D gradient(const Vector3D& x) const;

        // Returns Laplacian at given position x.
        virtual double laplacian(const Vector3D& x) const;

        // Returns sampler function object.
        virtual std::function<double(const Vector3D&)> sampler() const;
    };

    // Shared pointer for the ScalarField3 type.
    typedef std::shared_ptr<ScalarField3> ScalarField3Ptr;

    ScalarField3::ScalarField3() {
    }

    ScalarField3::~ScalarField3() {
    }

    Vector3D ScalarField3::gradient(const Vector3D&) const {
        return Vector3D();
    }

    double ScalarField3::laplacian(const Vector3D&) const {
        return 0.0;
    }

    std::function<double(const Vector3D&)> ScalarField3::sampler() const {
        const ScalarField3* self = this;
        return [self](const Vector3D& x) -> double {
            return self->sample(x);
        };
    }

}  // namespace jet
