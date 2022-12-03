#ifndef INCLUDE_JET_CPP_UTILS_H_
#define INCLUDE_JET_CPP_UTILS_H_

#include <algorithm>

namespace jet {

    template <class ForwardIt, class T, class Compare = std::less<T>>
    ForwardIt binaryFind(ForwardIt first, ForwardIt last, const T& value,
        Compare comp = {});

    template <class ForwardIt, class T, class Compare>
    ForwardIt binaryFind(ForwardIt first, ForwardIt last, const T& value,
        Compare comp) {
        // Note: BOTH type T and the type after ForwardIt is dereferenced
        // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
        // This is stricter than lower_bound requirement (see above)

        first = std::lower_bound(first, last, value, comp);
        return first != last && !comp(value, *first) ? first : last;
    }
}  // namespace jet


#endif  // INCLUDE_JET_CPP_UTILS_H_