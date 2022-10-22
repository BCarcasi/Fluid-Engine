#ifndef INCLUDE_JET_VECTOR_FIELD2_H_
#define INCLUDE_JET_VECTOR_FIELD2_H_

#include "field2.h"
#include "vector2.h"
#include <functional>
#include <memory>
#include "pch.h"

namespace jet {

    //! Abstract base class for 2-D vector field.
    class VectorField2 : public Field2 {
    public:
        //! Default constructor.
        VectorField2();

        //! Default destructor.
        virtual ~VectorField2();

        //! Returns sampled value at given position \p x.
        virtual Vector2D sample(const Vector2D& x) const = 0;

        //! Returns divergence at given position \p x.
        virtual double divergence(const Vector2D& x) const;

        //! Returns curl at given position \p x.
        virtual double curl(const Vector2D& x) const;

        //! Returns sampler function object.
        virtual std::function<Vector2D(const Vector2D&)> sampler() const;
    };

    //! Shared pointer for the VectorField2 type.
    typedef std::shared_ptr<VectorField2> VectorField2Ptr;

    VectorField2::VectorField2() {
    }

    VectorField2::~VectorField2() {
    }

    double VectorField2::divergence(const Vector2D&) const {
        return 0.0;
    }

    double VectorField2::curl(const Vector2D&) const {
        return 0.0;
    }

    std::function<Vector2D(const Vector2D&)> VectorField2::sampler() const {
        const VectorField2* self = this;
        return [self](const Vector2D& x) -> Vector2D {
            return self->sample(x);
        };
    }

}  // namespace jet

#endif  // INCLUDE_JET_VECTOR_FIELD2_H_