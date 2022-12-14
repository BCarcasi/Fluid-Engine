#ifndef INCLUDE_JET_VECTOR_FIELD3_H_
#define INCLUDE_JET_VECTOR_FIELD3_H_

#include "field3.h"
#include "vector3.h"
#include <functional>
#include <memory>

namespace jet {

    // Abstract base class for 3-D vector field.
    class VectorField3 : public Field3 {
    public:
        // Default constructor.
        VectorField3();

        // Default destructor.
        virtual ~VectorField3();

        // Returns sampled value at given position x.
        virtual Vector3D sample(const Vector3D& x) const = 0;

        // Returns divergence at given position x.
        virtual double divergence(const Vector3D& x) const;

        // Returns curl at given position x.
        virtual Vector3D curl(const Vector3D& x) const;

        // Returns sampler function object.
        virtual std::function<Vector3D(const Vector3D&)> sampler() const;
    };

    // Shared pointer for the VectorField3 type.
    typedef std::shared_ptr<VectorField3> VectorField3Ptr;

    VectorField3::VectorField3() {
    }

    VectorField3::~VectorField3() {
    }

    double VectorField3::divergence(const Vector3D&) const {
        return 0.0;
    }

    Vector3D VectorField3::curl(const Vector3D&) const {
        return Vector3D();
    }

    std::function<Vector3D(const Vector3D&)> VectorField3::sampler() const {
        const VectorField3* self = this;
        return [self](const Vector3D& x) -> Vector3D {
            return self->sample(x);
        };
    }

}  // namespace jet

#endif  // INCLUDE_JET_VECTOR_FIELD3_H_