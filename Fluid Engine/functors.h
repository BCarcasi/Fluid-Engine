#ifndef INCLUDE_JET_FUNCTORS_H_
#define INCLUDE_JET_FUNCTORS_H_

#include <functional>

namespace jet {

    //! Type casting operator.
    template <typename T, typename U>
    struct TypeCast {
        constexpr U operator()(const T& a) const;
    };

    //! Reverse minus operator.
    template <typename T>
    struct RMinus {
        constexpr T operator()(const T& a, const T& b) const;
    };

    //! Reverse divides operator.
    template <typename T>
    struct RDivides {
        constexpr T operator()(const T& a, const T& b) const;
    };

    template <typename T, typename U>
    constexpr U TypeCast<T, U>::operator()(const T& a) const {
        return static_cast<U>(a);
    }

    template <typename T>
    constexpr T RMinus<T>::operator()(const T& a, const T& b) const {
        return b - a;
    }

    template <typename T>
    constexpr T RDivides<T>::operator()(const T& a, const T& b) const {
        return b / a;
    }
}


#endif  // INCLUDE_JET_FUNCTORS_H_