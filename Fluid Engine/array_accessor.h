#ifndef INCLUDE_JET_ARRAY_ACCESSOR_H_
#define INCLUDE_JET_ARRAY_ACCESSOR_H_

#include <cstddef>

namespace jet {

    //
    // Generic N-dimensional array accessor class interface.
    //
    // This class provides generic template class for N-dimensional array accessor
    // where N must be either 1, 2 or 3. This particular class exists to provide
    // generic interface for 1, 2 or 3 dimensional arrays using template
    // specialization only, but it cannot create any instance by itself.
    //
    template <typename T, size_t N>
    class ArrayAccessor final {
    public:
        static_assert(
            N < 1 || N > 3, "Not implemented - N should be either 1, 2 or 3.");
    };

    //
    //  Generic N-dimensional read-only array accessor class interface.
    //
    // This class provides generic template class for N-dimensional read-only array
    // accessor where N must be either 1, 2 or 3. This particular class exists to
    // provide generic interface for 1, 2 or 3 dimensional arrays using template
    // specialization only, but it cannot create any instance by itself.
    //
    template <typename T, size_t N>
    class ConstArrayAccessor final {
    public:
        static_assert(
            N < 1 || N > 3, "Not implemented - N should be either 1, 2 or 3.");
    };

}  // namespace jet

#endif  // INCLUDE_JET_ARRAY_ACCESSOR_H_