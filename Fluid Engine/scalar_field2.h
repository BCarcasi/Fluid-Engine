#ifndef INCLUDE_JET_SCALAR_FIELD2_H_
#define INCLUDE_JET_SCALAR_FIELD2_H_

#include "field2.h"
#include "vector2.h"
#include <functional>
#include <memory>

#include "pch.h"

namespace jet {

    //! Abstract base class for 2-D scalar field.
    class ScalarField2 : public Field2 {
    public:
        //! Default constructor.
        ScalarField2();

        //! Default destructor.
        virtual ~ScalarField2();

        //! Returns sampled value at given position \p x.
        virtual double sample(const Vector2D& x) const = 0;

        //! Returns gradient vector at given position \p x.
        virtual Vector2D gradient(const Vector2D& x) const;

        //! Returns Laplacian at given position \p x.
        virtual double laplacian(const Vector2D& x) const;

        //! Returns sampler function object.
        virtual std::function<double(const Vector2D&)> sampler() const;
    };

    //! Shared pointer for the ScalarField2 type.
    typedef std::shared_ptr<ScalarField2> ScalarField2Ptr;

    ScalarField2::ScalarField2() {
    }

    ScalarField2::~ScalarField2() {
    }

    Vector2D ScalarField2::gradient(const Vector2D&) const {
        return Vector2D();
    }

    double ScalarField2::laplacian(const Vector2D&) const {
        return 0.0;
    }

    std::function<double(const Vector2D&)> ScalarField2::sampler() const {
        const ScalarField2* self = this;
        return [self](const Vector2D& x) -> double {
            return self->sample(x);
        };
    }

}  // namespace jet

#endif  // INCLUDE_JET_SCALAR_FIELD2_H_