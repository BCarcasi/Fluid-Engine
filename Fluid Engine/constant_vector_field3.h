#pragma once

#ifndef INCLUDE_JET_CONSTANT_VECTOR_FIELD3_H_
#define INCLUDE_JET_CONSTANT_VECTOR_FIELD3_H_

#include "vector_field3.h"
#include <memory>
#include "private_helpers.h"

namespace jet {

    // 3-D constant vector field.
    class ConstantVectorField3 final : public VectorField3 {
    public:
        class Builder;

        // Constructs a constant vector field with given value.
        explicit ConstantVectorField3(const Vector3D& value);

        // Returns the sampled value at given position x.
        Vector3D sample(const Vector3D& x) const override;

        // Returns the sampler function.
        std::function<Vector3D(const Vector3D&)> sampler() const override;

        // Returns builder fox ConstantVectorField3.
        static Builder builder();

    private:
        Vector3D _value;
    };

    // Shared pointer for the ConstantVectorField3 type.
    typedef std::shared_ptr<ConstantVectorField3> ConstantVectorField3Ptr;


    //
    // Front-end to create ConstantVectorField3 objects step by step.
    //
    class ConstantVectorField3::Builder final {
    public:
        // Returns builder with value.
        Builder& withValue(const Vector3D& value);

        // Builds ConstantVectorField3.
        ConstantVectorField3 build() const;

        // Builds shared pointer of ConstantVectorField3 instance.
        ConstantVectorField3Ptr makeShared() const;

    private:
        Vector3D _value{ 0, 0, 0 };
    };


    ConstantVectorField3::ConstantVectorField3(const Vector3D& value) :
        _value(value) {
    }

    Vector3D ConstantVectorField3::sample(const Vector3D& x) const {
        UNUSED_VARIABLE(x);

        return _value;
    }

    std::function<Vector3D(const Vector3D&)> ConstantVectorField3::sampler() const {
        return [this](const Vector3D&) -> Vector3D {
            return _value;
        };
    }

    ConstantVectorField3::Builder ConstantVectorField3::builder() {
        return Builder();
    }


    ConstantVectorField3::Builder&
        ConstantVectorField3::Builder::withValue(const Vector3D& value) {
        _value = value;
        return *this;
    }

    ConstantVectorField3 ConstantVectorField3::Builder::build() const {
        return ConstantVectorField3(_value);
    }

    ConstantVectorField3Ptr ConstantVectorField3::Builder::makeShared() const {
        return std::shared_ptr<ConstantVectorField3>(
            new ConstantVectorField3(_value),
            [](ConstantVectorField3* obj) {
                delete obj;
            });
    }
}  // namespace jet

#endif  // INCLUDE_JET_CONSTANT_VECTOR_FIELD3_H_