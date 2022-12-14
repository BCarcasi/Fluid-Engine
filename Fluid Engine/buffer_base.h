#ifndef FLATBUFFERS_BASE_H_
#define FLATBUFFERS_BASE_H_

#include <assert.h>

#ifndef ARDUINO
#include <cstdint>
#endif
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string>
#ifndef ARDUINO
#include <utility>
#else
#include <utility.h>
#endif
#include <type_traits>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <memory>

#ifdef _STLPORT_VERSION
#define FLATBUFFERS_CPP98_STL
#endif
#ifndef FLATBUFFERS_CPP98_STL
#include <functional>
#endif

/// @cond FLATBUFFERS_INTERNAL
#if __cplusplus <= 199711L && \
    (!defined(_MSC_VER) || _MSC_VER < 1600) && \
    (!defined(__GNUC__) || \
      (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__ < 40400))
#error A C++11 compatible compiler with support for the auto typing is \
         required for FlatBuffers.
#error __cplusplus _MSC_VER __GNUC__  __GNUC_MINOR__  __GNUC_PATCHLEVEL__
#endif

#if !defined(__clang__) && \
    defined(__GNUC__) && \
    (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__ < 40600)
  // Backwards compatability for g++ 4.4, and 4.5 which don't have the nullptr
  // and constexpr keywords. Note the __clang__ check is needed, because clang
  // presents itself as an older GNUC compiler.
#ifndef nullptr_t
const class nullptr_t {
public:
    template<class T> inline operator T* () const { return 0; }
private:
    void operator&() const;
} nullptr = {};
#endif
#ifndef constexpr
#define constexpr const
#endif
#endif

// The wire format uses a little endian encoding (since that's efficient for
// the common platforms).
#if defined(__s390x__)
#define FLATBUFFERS_LITTLEENDIAN 0
#endif // __s390x__
#if !defined(FLATBUFFERS_LITTLEENDIAN)
#if defined(__GNUC__) || defined(__clang__)
#ifdef __BIG_ENDIAN__
#define FLATBUFFERS_LITTLEENDIAN 0
#else
#define FLATBUFFERS_LITTLEENDIAN 1
#endif // __BIG_ENDIAN__
#elif defined(_MSC_VER)
#if defined(_M_PPC)
#define FLATBUFFERS_LITTLEENDIAN 0
#else
#define FLATBUFFERS_LITTLEENDIAN 1
#endif
#else
#error Unable to determine endianness, define FLATBUFFERS_LITTLEENDIAN.
#endif
#endif // !defined(FLATBUFFERS_LITTLEENDIAN)

#define FLATBUFFERS_VERSION_MAJOR 1
#define FLATBUFFERS_VERSION_MINOR 7
#define FLATBUFFERS_VERSION_REVISION 0
#define FLATBUFFERS_STRING_EXPAND(X) #X
#define FLATBUFFERS_STRING(X) FLATBUFFERS_STRING_EXPAND(X)

#if (!defined(_MSC_VER) || _MSC_VER > 1600) && \
    (!defined(__GNUC__) || (__GNUC__ * 100 + __GNUC_MINOR__ >= 407))
#define FLATBUFFERS_FINAL_CLASS final
#define FLATBUFFERS_OVERRIDE override
#else
#define FLATBUFFERS_FINAL_CLASS
#define FLATBUFFERS_OVERRIDE
#endif

#if (!defined(_MSC_VER) || _MSC_VER >= 1900) && \
    (!defined(__GNUC__) || (__GNUC__ * 100 + __GNUC_MINOR__ >= 406))
#define FLATBUFFERS_CONSTEXPR constexpr
#else
#define FLATBUFFERS_CONSTEXPR
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) && __GNUC__ * 10 + __GNUC_MINOR__ >= 46 || \
    defined(_MSC_FULL_VER) && _MSC_FULL_VER >= 190023026
#define FLATBUFFERS_NOEXCEPT noexcept
#else
#define FLATBUFFERS_NOEXCEPT
#endif

// NOTE: the FLATBUFFERS_DELETE_FUNC macro may change the access mode to
// private, so be sure to put it at the end or reset access mode explicitly.
#if (!defined(_MSC_VER) || _MSC_FULL_VER >= 180020827) && \
    (!defined(__GNUC__) || (__GNUC__ * 100 + __GNUC_MINOR__ >= 404))
#define FLATBUFFERS_DELETE_FUNC(func) func = delete;
#else
#define FLATBUFFERS_DELETE_FUNC(func) private: func;
#endif

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable: 4127) // C4127: conditional expression is constant
#endif

/// @endcond

/// @file
namespace flatbuffers {

    /// @cond FLATBUFFERS_INTERNAL
    // Our default offset / size type, 32bit on purpose on 64bit systems.
    // Also, using a consistent offset type maintains compatibility of serialized
    // offset values between 32bit and 64bit systems.
    typedef uint32_t uoffset_t;

    // Signed offsets for references that can go in both directions.
    typedef int32_t soffset_t;

    // Offset/index used in v-tables, can be changed to uint8_t in
    // format forks to save a bit of space if desired.
    typedef uint16_t voffset_t;

    typedef uintmax_t largest_scalar_t;

    // In 32bits, this evaluates to 2GB - 1
#define FLATBUFFERS_MAX_BUFFER_SIZE ((1ULL << (sizeof(soffset_t) * 8 - 1)) - 1)

// We support aligning the contents of buffers up to this size.
#define FLATBUFFERS_MAX_ALIGNMENT 16

    template<typename T> T EndianScalar(T t) {
#if FLATBUFFERS_LITTLEENDIAN
        return t;
#else
        return EndianSwap(t);
#endif
    }

    template<typename T> T ReadScalar(const void* p) {
        return EndianScalar(*reinterpret_cast<const T*>(p));
    }

    template<typename T> void WriteScalar(void* p, T t) {
        *reinterpret_cast<T*>(p) = EndianScalar(t);
    }

    // Computes how many bytes you'd have to pad to be able to write an
    // "scalar_size" scalar if the buffer had grown to "buf_size" (downwards in
    // memory).
    inline size_t PaddingBytes(size_t buf_size, size_t scalar_size) {
        return ((~buf_size) + 1) & (scalar_size - 1);
    }

}
#endif  // FLATBUFFERS_BASE_H_