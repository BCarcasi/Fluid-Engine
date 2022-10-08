#ifndef INCLUDE_JET_MATH_UTILS_H_
#define INCLUDE_JET_MATH_UTILS_H_

#include "macros.h"
#include <cstddef>
#include <limits>
#include "constants.h"
#include <algorithm>
#include <cmath>


namespace jet {

    //
    //      Returns true if x and y are similar.
    //
    // \param[in]  x     The first value.
    // \param[in]  y     The second value.
    // \param[in]  eps   The tolerance.
    //
    //     T     Value type.
    //
    // \return     True if similar.
    //
    template <typename T>
    inline bool similar(T x, T y, T eps = std::numeric_limits<T>::epsilon());

    //
    //      Returns the sign of the value.
    //
    // \param[in]  x     Input value.
    //
    //     T     Value type.
    //
    // \return     The sign.
    //
    template <typename T>
    inline T sign(T x);

    //
    //      Returns the minimum value among three inputs.
    //
    // \param[in]  x     The first value.
    // \param[in]  y     The second value.
    // \param[in]  z     The three value.
    //
    //     T     Value type.
    //
    // \return     The minimum value.
    //
    template <typename T>
    inline T min3(T x, T y, T z);

    //
    //      Returns the maximum value among three inputs.
    //
    // \param[in]  x     The first value.
    // \param[in]  y     The second value.
    // \param[in]  z     The three value.
    //
    //     T     Value type.
    //
    // \return     The maximum value.
    //
    template <typename T>
    inline T max3(T x, T y, T z);

    // Returns minimum among n-elements.
    template <typename T>
    inline T minn(const T* x, size_t n);

    // Returns maximum among n-elements.
    template <typename T>
    inline T maxn(const T* x, size_t n);

    //
    //      Returns the absolute minimum value among the two inputs.
    //
    // \param[in]  x     The first value.
    // \param[in]  y     The second value.
    //
    //     T     Value type.
    //
    // \return     The absolute minimum.
    //
    template <typename T>
    inline T absmin(T x, T y);

    //
    //      Returns the absolute maximum value among the two inputs.
    //
    // \param[in]  x     The first value.
    // \param[in]  y     The second value.
    //
    //     T     Value type.
    //
    // \return     The absolute maximum.
    //
    template <typename T>
    inline T absmax(T x, T y);

    // Returns absolute minimum among n-elements.
    template <typename T>
    inline T absminn(const T* x, size_t n);

    // Returns absolute maximum among n-elements.
    template <typename T>
    inline T absmaxn(const T* x, size_t n);

    template <typename T>
    inline size_t argmin2(T x, T y);

    template <typename T>
    inline size_t argmax2(T x, T y);

    template <typename T>
    inline size_t argmin3(T x, T y, T z);

    template <typename T>
    inline size_t argmax3(T x, T y, T z);

    //
    //      Returns the square of x.
    //
    // \param[in]  x     The input.
    //
    //     T     Value type.
    //
    // \return     The squared value.
    //
    template <typename T>
    inline T square(T x);

    //
    //      Returns the cubic of x.
    //
    // \param[in]  x     The input.
    //
    //     T     Value type.
    //
    // \return     The cubic of x.
    //
    template <typename T>
    inline T cubic(T x);

    //
    //      Returns the clamped value.
    //
    // \param[in]  val   The value.
    // \param[in]  low   The low value.
    // \param[in]  high  The high value.
    //
    //     T     Value type.
    //
    // \return     The clamped value.
    //
    template <typename T>
    inline T clamp(T val, T low, T high);

    //
    //      Converts degrees to radians.
    //
    // \param[in]  angleInDegrees The angle in degrees.
    //
    //     T              Value type.
    //
    // \return     Angle in radians.
    //
    template <typename T>
    inline T degreesToRadians(T angleInDegrees);

    //
    //      Converts radians to degrees.
    //
    // \param[in]  angleInDegrees The angle in radians.
    //
    //     T              Value type.
    //
    // \return     Angle in degrees.
    //
    template <typename T>
    inline T radiansToDegrees(T angleInRadians);

    //
    //      Gets the barycentric coordinate.
    //
    // \param[in]  x     The input value.
    // \param[in]  iLow  The lowest index.
    // \param[in]  iHigh The highest index.
    //      i     The output index.
    //      t     The offset from i.
    //
    //     T     Value type.
    //
    template <class T>
    inline void getBarycentric(T x, ssize_t iLow, ssize_t iHigh, ssize_t* i, T* t);

    //
    //      Computes linear interpolation.
    //
    // \param[in]  f0    The first value.
    // \param[in]  f1    The second value.
    // \param[in]  t     Relative offset [0, 1] from the first value.
    //
    //     S     Input value type.
    //     T     Offset type.
    //
    // \return     The interpolated value.
    //
    template <typename S, typename T>
    inline S lerp(const S& f0, const S& f1, T t);

    //      Computes bilinear interpolation.
    template <typename S, typename T>
    inline S bilerp(const S& f00, const S& f10, const S& f01, const S& f11, T tx,
        T ty);

    //      Computes trilinear interpolation.
    template <typename S, typename T>
    inline S trilerp(const S& f000, const S& f100, const S& f010, const S& f110,
        const S& f001, const S& f101, const S& f011, const S& f111,
        T tx, T ty, T tz);

    //      Computes Catmull-Rom interpolation.
    template <typename S, typename T>
    inline S catmullRom(const S& f0, const S& f1, const S& f2, const S& f3, T t);

    //      Computes monotonic Catmull-Rom interpolation.
    template <typename T>
    inline T monotonicCatmullRom(const T& f0, const T& f1, const T& f2, const T& f3,
        T t);

    template <typename T>
    inline bool similar(T x, T y, T eps) {
        return (std::abs(x - y) <= eps);
    }

    template <typename T>
    inline T sign(T x) {
        if (x >= 0) {
            return 1;
        }
        else {
            return -1;
        }
    }

    template <typename T>
    inline T min3(T x, T y, T z) {
        return std::min(std::min(x, y), z);
    }

    template <typename T>
    inline T max3(T x, T y, T z) {
        return std::max(std::max(x, y), z);
    }

    template <typename T>
    inline T minn(const T* x, size_t n) {
        T m = x[0];
        for (size_t i = 1; i < n; i++) {
            m = std::min(m, x[i]);
        }
        return m;
    }

    template <typename T>
    inline T maxn(const T* x, size_t n) {
        T m = x[0];
        for (size_t i = 1; i < n; i++) {
            m = std::max(m, x[i]);
        }
        return m;
    }

    template <typename T>
    inline T absmin(T x, T y) {
        return (x * x < y* y) ? x : y;
    }

    template <typename T>
    inline T absmax(T x, T y) {
        return (x * x > y * y) ? x : y;
    }

    template <typename T>
    inline T absminn(const T* x, size_t n) {
        T m = x[0];
        for (size_t i = 1; i < n; i++) {
            m = absmin(m, x[i]);
        }
        return m;
    }

    template <typename T>
    inline T absmaxn(const T* x, size_t n) {
        T m = x[0];
        for (size_t i = 1; i < n; i++) {
            m = absmax(m, x[i]);
        }
        return m;
    }

    template <typename T>
    inline size_t argmin2(T x, T y) {
        return (x < y) ? 0 : 1;
    }

    template <typename T>
    inline size_t argmax2(T x, T y) {
        return (x > y) ? 0 : 1;
    }

    template <typename T>
    inline size_t argmin3(T x, T y, T z) {
        if (x < y) {
            return (x < z) ? 0 : 2;
        }
        else {
            return (y < z) ? 1 : 2;
        }
    }

    template <typename T>
    inline size_t argmax3(T x, T y, T z) {
        if (x > y) {
            return (x > z) ? 0 : 2;
        }
        else {
            return (y > z) ? 1 : 2;
        }
    }

    template <typename T>
    inline T square(T x) {
        return x * x;
    }

    template <typename T>
    inline T cubic(T x) {
        return x * x * x;
    }

    template <typename T>
    inline T clamp(T val, T low, T high) {
        if (val < low) {
            return low;
        }
        else if (val > high) {
            return high;
        }
        else {
            return val;
        }
    }

    template <typename T>
    inline T degreesToRadians(T angleInDegrees) {
        return angleInDegrees * pi<T>() / 180;
    }

    template <typename T>
    inline T radiansToDegrees(T angleInRadians) {
        return angleInRadians * 180 / pi<T>();
    }

    template<typename T>
    inline void getBarycentric(
        T x,
        ssize_t iLow,
        ssize_t iHigh,
        ssize_t* i,
        T* f) {
        T s = std::floor(x);
        *i = static_cast<ssize_t>(s);

        ssize_t offset = -iLow;
        iLow += offset;
        iHigh += offset;

        if (iLow == iHigh) {
            *i = iLow;
            *f = 0;
        }
        else if (*i < iLow) {
            *i = iLow;
            *f = 0;
        }
        else if (*i > iHigh - 1) {
            *i = iHigh - 1;
            *f = 1;
        }
        else {
            *f = static_cast<T>(x - s);
        }

        *i -= offset;
    }

    template<typename S, typename T>
    inline S lerp(const S& value0, const S& value1, T f) {
        return (1 - f) * value0 + f * value1;
    }

    template<typename S, typename T>
    inline S bilerp(
        const S& f00,
        const S& f10,
        const S& f01,
        const S& f11,
        T tx, T ty) {
        return lerp(
            lerp(f00, f10, tx),
            lerp(f01, f11, tx),
            ty);
    }

    template<typename S, typename T>
    inline S trilerp(
        const S& f000,
        const S& f100,
        const S& f010,
        const S& f110,
        const S& f001,
        const S& f101,
        const S& f011,
        const S& f111,
        T tx,
        T ty,
        T fz) {
        return lerp(
            bilerp(f000, f100, f010, f110, tx, ty),
            bilerp(f001, f101, f011, f111, tx, ty),
            fz);
    }

    template <typename S, typename T>
    inline S catmullRom(
        const S& f0,
        const S& f1,
        const S& f2,
        const S& f3,
        T f) {
        S d1 = (f2 - f0) / 2;
        S d2 = (f3 - f1) / 2;
        S D1 = f2 - f1;

        S a3 = d1 + d2 - 2 * D1;
        S a2 = 3 * D1 - 2 * d1 - d2;
        S a1 = d1;
        S a0 = f1;

        return a3 * cubic(f) + a2 * square(f) + a1 * f + a0;
    }

    template <typename T>
    inline T monotonicCatmullRom(
        const T& f0,
        const T& f1,
        const T& f2,
        const T& f3,
        T f) {
        T d1 = (f2 - f0) / 2;
        T d2 = (f3 - f1) / 2;
        T D1 = f2 - f1;

        if (std::fabs(D1) < kEpsilonD) {
            d1 = d2 = 0;
        }

        if (sign(D1) != sign(d1)) {
            d1 = 0;
        }

        if (sign(D1) != sign(d2)) {
            d2 = 0;
        }

        T a3 = d1 + d2 - 2 * D1;
        T a2 = 3 * D1 - 2 * d1 - d2;
        T a1 = d1;
        T a0 = f1;

        return a3 * cubic(f) + a2 * square(f) + a1 * f + a0;
    }


}  // namespace jet



#endif  // INCLUDE_JET_MATH_UTILS_H_