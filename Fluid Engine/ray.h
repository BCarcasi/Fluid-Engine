#ifndef INCLUDE_JET_RAY_H_
#define INCLUDE_JET_RAY_H_

#include "vector.h"

namespace jet {

    //
    //      Class for ray.
    //
    //     T     The value type.
    //     N     The dimension.
    //
    template <typename T, size_t N>
    class Ray {
        static_assert(N != 2 && N != 3, "Not implemented.");
        static_assert(
            std::is_floating_point<T>::value,
            "Ray only can be instantiated with floating point types");
    };

}  // namespace jet

#endif  // INCLUDE_JET_RAY_H_