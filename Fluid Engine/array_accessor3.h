#ifndef INCLUDE_JET_ARRAY_ACCESSOR3_H_
#define INCLUDE_JET_ARRAY_ACCESSOR3_H_

#include "array_accessor.h"
#include "point3.h"
#include "size3.h"
#include <utility>
#include "macros.h"
#include "parallel.h"


namespace jet {

    //
    // 3-D array accessor class.
    //
    // This class represents 3-D array accessor. Array accessor provides array-like
    // data read/write functions, but does not handle memory management. Thus, it
    // is more like a random access iterator, but with multi-dimension support.
    // Similar to Array<T, 3>, this class interprets a linear array as a 3-D array
    // using i-major indexing.
    //
    // Array<T, 3>
    //
    //  T - Array value type.
    //
    template <typename T>
    class ArrayAccessor<T, 3> final {
    public:
        // Constructs empty 3-D array accessor.
        ArrayAccessor();

        // Constructs an array accessor that wraps given array.
        //size Size of the 3-D array.
        //data Raw array pointer.
        ArrayAccessor(const Size3& size, T* const data);

        // Constructs an array accessor that wraps given array.
        //width Width of the 3-D array.
        //height Height of the 3-D array.
        //depth Depth of the 3-D array.
        //data Raw array pointer.
        ArrayAccessor(
            size_t width, size_t height, size_t depth, T* const data);

        // Copy constructor.
        ArrayAccessor(const ArrayAccessor& other);

        // Replaces the content with given other array accessor.
        void set(const ArrayAccessor& other);

        // Resets the array.
        void reset(const Size3& size, T* const data);

        // Resets the array.
        void reset(size_t width, size_t height, size_t depth, T* const data);

        // Returns the reference to the i-th element.
        T& at(size_t i);

        // Returns the const reference to the i-th element.
        const T& at(size_t i) const;

        // Returns the reference to the element at (pt.x, pt.y, pt.z).
        T& at(const Point3UI& pt);

        // Returns the const reference to the element at (pt.x, pt.y, pt.z).
        const T& at(const Point3UI& pt) const;

        // Returns the reference to the element at (i, j, k).
        T& at(size_t i, size_t j, size_t k);

        // Returns the const reference to the element at (i, j, k).
        const T& at(size_t i, size_t j, size_t k) const;

        // Returns the begin iterator of the array.
        T* const begin() const;

        // Returns the end iterator of the array.
        T* const end() const;

        // Returns the begin iterator of the array.
        T* begin();

        // Returns the end iterator of the array.
        T* end();

        // Returns the size of the array.
        Size3 size() const;

        // Returns the width of the array.
        size_t width() const;

        // Returns the height of the array.
        size_t height() const;

        // Returns the depth of the array.
        size_t depth() const;

        // Returns the raw pointer to the array data.
        T* const data() const;

        // Swaps the content of with other array accessor.
        void swap(ArrayAccessor& other);

        //
        // Iterates the array and invoke given func for each index.
        //
        // This function iterates the array elements and invoke the callback
        // function func. The callback function takes array's element as its
        // input. 
        //
        template <typename Callback>
        void forEach(Callback func) const;

        //
        // Iterates the array and invoke given func for each index.
        //
        // This function iterates the array elements and invoke the callback
        // function func. The callback function takes three parameters which are
        // the (i, j, j) indices of the array. 
        //
        template <typename Callback>
        void forEachIndex(Callback func) const;

        //
        // Iterates the array and invoke given func for each index in
        //     parallel.
        //
        // This function iterates the array elements and invoke the callback
        // function func. The callback function takes array's element as its
        // input. The order of execution will be non-deterministic since it runs in
        // parallel. 
        //
        template <typename Callback>
        void parallelForEach(Callback func);

        //
        // Iterates the array and invoke given func for each index in
        //     parallel using multi-threading.
        //
        // This function iterates the array elements and invoke the callback
        // function func in parallel using multi-threading. The callback
        // function takes two parameters which are the (i, j, k) indices of the
        // array. The order of execution will be non-deterministic since it runs in
        // parallel.
        //
        template <typename Callback>
        void parallelForEachIndex(Callback func) const;

        // Returns the linear index of the given 3-D coordinate (pt.x, pt.y, pt.z).
        size_t index(const Point3UI& pt) const;

        // Returns the linear index of the given =3-D coordinate (i, j, k).
        size_t index(size_t i, size_t j, size_t k) const;

        // Returns the reference to the i-th element.
        T& operator[](size_t i);

        // Returns the const reference to the i-th element.
        const T& operator[](size_t i) const;

        // Returns the reference to the element at (pt.x, pt.y, pt.z).
        T& operator()(const Point3UI& pt);

        // Returns the const reference to the element at (pt.x, pt.y, pt.z).
        const T& operator()(const Point3UI& pt) const;

        // Returns the reference to the element at (i, j, k).
        T& operator()(size_t i, size_t j, size_t k);

        // Returns the const reference to the element at (i, j, k).
        const T& operator()(size_t i, size_t j, size_t k) const;

        // Copies given array other to this array.
        ArrayAccessor& operator=(const ArrayAccessor& other);

        // Casts type to ConstArrayAccessor.
        operator ConstArrayAccessor<T, 3>() const;

    private:
        Size3 _size;
        T* _data;
    };

    // Type alias for 3-D array accessor.
    template <typename T> using ArrayAccessor3 = ArrayAccessor<T, 3>;


    //
    // 3-D read-only array accessor class.
    //
    // This class represents 3-D read-only array accessor. Array accessor provides
    // array-like data read/write functions, but does not handle memory management.
    // Thus, it is more like a random access iterator, but with multi-dimension
    // support.Similar to Array<T, 3>, this class interprets a linear array as a
    // 3-D array using i-major indexing.
    //
    template <typename T>
    class ConstArrayAccessor<T, 3> {
    public:
        // Constructs empty 3-D read-only array accessor.
        ConstArrayAccessor();

        // Constructs a read-only array accessor that wraps given array.
        //size Size of the 3-D array.
        //data Raw array pointer.
        ConstArrayAccessor(const Size3& size, const T* const data);

        // Constructs a read-only array accessor that wraps given array.
        //width Width of the 3-D array.
        //height Height of the 3-D array.
        //depth Depth of the 3-D array.
        //data Raw array pointer.
        ConstArrayAccessor(
            size_t width, size_t height, size_t depth, const T* const data);

        // Constructs a read-only array accessor from read/write accessor.
        explicit ConstArrayAccessor(const ArrayAccessor<T, 3>& other);

        // Copy constructor.
        ConstArrayAccessor(const ConstArrayAccessor& other);

        // Returns the const reference to the i-th element.
        const T& at(size_t i) const;

        // Returns the const reference to the element at (pt.x, pt.y, pt.z).
        const T& at(const Point3UI& pt) const;

        // Returns the const reference to the element at (i, j, k).
        const T& at(size_t i, size_t j, size_t k) const;

        // Returns the begin iterator of the array.
        const T* const begin() const;

        // Returns the end iterator of the array.
        const T* const end() const;

        // Returns the size of the array.
        Size3 size() const;

        // Returns the width of the array.
        size_t width() const;

        // Returns the height of the array.
        size_t height() const;

        // Returns the depth of the array.
        size_t depth() const;

        // Returns the raw pointer to the array data.
        const T* const data() const;

        //
        // Iterates the array and invoke given func for each index.
        //
        // This function iterates the array elements and invoke the callback
        // function func. The callback function takes array's element as its
        // input.
        //
        template <typename Callback>
        void forEach(Callback func) const;

        //
        // Iterates the array and invoke given func for each index.
        //
        // This function iterates the array elements and invoke the callback
        // function func. The callback function takes three parameters which are
        // the (i, j, j) indices of the array.
        //
        template <typename Callback>
        void forEachIndex(Callback func) const;

        //
        // Iterates the array and invoke given func for each index in
        //     parallel using multi-threading.
        //
        // This function iterates the array elements and invoke the callback
        // function func in parallel using multi-threading. The callback
        // function takes two parameters which are the (i, j, k) indices of the
        // array. The order of execution will be non-deterministic since it runs in
        // parallel.
        //
        template <typename Callback>
        void parallelForEachIndex(Callback func) const;

        // Returns the linear index of the given 3-D coordinate (pt.x, pt.y, pt.z).
        size_t index(const Point3UI& pt) const;

        // Returns the linear index of the given =3-D coordinate (i, j, k).
        size_t index(size_t i, size_t j, size_t k) const;

        // Returns the const reference to the i-th element.
        const T& operator[](size_t i) const;

        // Returns the const reference to the element at (pt.x, pt.y, pt.z).
        const T& operator()(const Point3UI& pt) const;

        // Returns the reference to the element at (i, j, k).
        const T& operator()(size_t i, size_t j, size_t k) const;

    private:
        Size3 _size;
        const T* _data;
    };

    // Type alias for 3-D const array accessor.
    template <typename T> using ConstArrayAccessor3 = ConstArrayAccessor<T, 3>;

    template <typename T>
    ArrayAccessor<T, 3>::ArrayAccessor() : _data(nullptr) {
    }

    template <typename T>
    ArrayAccessor<T, 3>::ArrayAccessor(const Size3& size, T* const data) {
        reset(size, data);
    }

    template <typename T>
    ArrayAccessor<T, 3>::ArrayAccessor(
        size_t width, size_t height, size_t depth, T* const data) {
        reset(width, height, depth, data);
    }

    template <typename T>
    ArrayAccessor<T, 3>::ArrayAccessor(const ArrayAccessor& other) {
        set(other);
    }

    template <typename T>
    void ArrayAccessor<T, 3>::set(const ArrayAccessor& other) {
        reset(other._size, other._data);
    }

    template <typename T>
    void ArrayAccessor<T, 3>::reset(const Size3& size, T* const data) {
        _size = size;
        _data = data;
    }

    template <typename T>
    void ArrayAccessor<T, 3>::reset(
        size_t width, size_t height, size_t depth, T* const data) {
        reset(Size3(width, height, depth), data);
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::at(size_t i) {
        JET_ASSERT(i < _size.x* _size.y* _size.z);
        return _data[i];
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::at(size_t i) const {
        JET_ASSERT(i < _size.x* _size.y* _size.z);
        return _data[i];
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::at(const Point3UI& pt) {
        return at(pt.x, pt.y, pt.z);
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::at(const Point3UI& pt) const {
        return at(pt.x, pt.y, pt.z);
    }

    template <typename T>
    T* const ArrayAccessor<T, 3>::begin() const {
        return _data;
    }

    template <typename T>
    T* const ArrayAccessor<T, 3>::end() const {
        return _data + _size.x * _size.y * _size.z;
    }

    template <typename T>
    T* ArrayAccessor<T, 3>::begin() {
        return _data;
    }

    template <typename T>
    T* ArrayAccessor<T, 3>::end() {
        return _data + _size.x * _size.y * _size.z;
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::operator()(const Point3UI& pt) {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y&& pt.z < _size.z);
        return _data[pt.x + _size.x * (pt.y + _size.y * pt.z)];
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::operator()(const Point3UI& pt) const {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y&& pt.z < _size.z);
        return _data[pt.x + _size.x * (pt.y + _size.y * pt.z)];
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::at(size_t i, size_t j, size_t k) {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::at(size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    Size3 ArrayAccessor<T, 3>::size() const {
        return _size;
    }

    template <typename T>
    size_t ArrayAccessor<T, 3>::width() const {
        return _size.x;
    }

    template <typename T>
    size_t ArrayAccessor<T, 3>::height() const {
        return _size.y;
    }

    template <typename T>
    size_t ArrayAccessor<T, 3>::depth() const {
        return _size.z;
    }

    template <typename T>
    T* const ArrayAccessor<T, 3>::data() const {
        return _data;
    }

    template <typename T>
    void ArrayAccessor<T, 3>::swap(ArrayAccessor& other) {
        std::swap(other._data, _data);
        std::swap(other._size, _size);
    }

    template <typename T>
    template <typename Callback>
    void ArrayAccessor<T, 3>::forEach(Callback func) const {
        for (size_t k = 0; k < _size.z; ++k) {
            for (size_t j = 0; j < _size.y; ++j) {
                for (size_t i = 0; i < _size.x; ++i) {
                    func(at(i, j, k));
                }
            }
        }
    }

    template <typename T>
    template <typename Callback>
    void ArrayAccessor<T, 3>::forEachIndex(Callback func) const {
        for (size_t k = 0; k < _size.z; ++k) {
            for (size_t j = 0; j < _size.y; ++j) {
                for (size_t i = 0; i < _size.x; ++i) {
                    func(i, j, k);
                }
            }
        }
    }

    template <typename T>
    template <typename Callback>
    void ArrayAccessor<T, 3>::parallelForEach(Callback func) {
        parallelFor(kZeroSize, _size.x, kZeroSize, _size.y, kZeroSize, _size.z,
            [&](size_t i, size_t j, size_t k) {
                func(at(i, j, k));
            });
    }

    template <typename T>
    template <typename Callback>
    void ArrayAccessor<T, 3>::parallelForEachIndex(Callback func) const {
        parallelFor(
            kZeroSize, _size.x, kZeroSize, _size.y, kZeroSize, _size.z, func);
    }

    template <typename T>
    size_t ArrayAccessor<T, 3>::index(const Point3UI& pt) const {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y&& pt.z < _size.z);
        return pt.x + _size.x * (pt.y + _size.y * pt.z);
    }

    template <typename T>
    size_t ArrayAccessor<T, 3>::index(size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return i + _size.x * (j + _size.y * k);
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::operator[](size_t i) {
        return _data[i];
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::operator[](size_t i) const {
        return _data[i];
    }

    template <typename T>
    T& ArrayAccessor<T, 3>::operator()(size_t i, size_t j, size_t k) {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    const T& ArrayAccessor<T, 3>::operator()(size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    ArrayAccessor<T, 3>& ArrayAccessor<T, 3>::operator=(
        const ArrayAccessor& other) {
        set(other);
        return *this;
    }

    template <typename T>
    ArrayAccessor<T, 3>::operator ConstArrayAccessor<T, 3>() const {
        return ConstArrayAccessor<T, 3>(*this);
    }


    template <typename T>
    ConstArrayAccessor<T, 3>::ConstArrayAccessor() : _data(nullptr) {
    }

    template <typename T>
    ConstArrayAccessor<T, 3>::ConstArrayAccessor(
        const Size3& size, const T* const data) {
        _size = size;
        _data = data;
    }

    template <typename T>
    ConstArrayAccessor<T, 3>::ConstArrayAccessor(
        size_t width, size_t height, size_t depth, const T* const data) {
        _size = Size3(width, height, depth);
        _data = data;
    }

    template <typename T>
    ConstArrayAccessor<T, 3>::ConstArrayAccessor(const ArrayAccessor<T, 3>& other) {
        _size = other.size();
        _data = other.data();
    }

    template <typename T>
    ConstArrayAccessor<T, 3>::ConstArrayAccessor(const ConstArrayAccessor& other) {
        _size = other._size;
        _data = other._data;
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::at(size_t i) const {
        JET_ASSERT(i < _size.x* _size.y* _size.z);
        return _data[i];
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::at(const Point3UI& pt) const {
        return at(pt.x, pt.y, pt.z);
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::at(size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    const T* const ConstArrayAccessor<T, 3>::begin() const {
        return _data;
    }

    template <typename T>
    const T* const ConstArrayAccessor<T, 3>::end() const {
        return _data + _size.x * _size.y * _size.z;
    }

    template <typename T>
    Size3 ConstArrayAccessor<T, 3>::size() const {
        return _size;
    }

    template <typename T>
    size_t ConstArrayAccessor<T, 3>::width() const {
        return _size.x;
    }

    template <typename T>
    size_t ConstArrayAccessor<T, 3>::height() const {
        return _size.y;
    }

    template <typename T>
    size_t ConstArrayAccessor<T, 3>::depth() const {
        return _size.z;
    }

    template <typename T>
    const T* const ConstArrayAccessor<T, 3>::data() const {
        return _data;
    }

    template <typename T>
    template <typename Callback>
    void ConstArrayAccessor<T, 3>::forEach(Callback func) const {
        for (size_t k = 0; k < _size.z; ++k) {
            for (size_t j = 0; j < _size.y; ++j) {
                for (size_t i = 0; i < _size.x; ++i) {
                    func(at(i, j, k));
                }
            }
        }
    }

    template <typename T>
    template <typename Callback>
    void ConstArrayAccessor<T, 3>::forEachIndex(Callback func) const {
        for (size_t k = 0; k < _size.z; ++k) {
            for (size_t j = 0; j < _size.y; ++j) {
                for (size_t i = 0; i < _size.x; ++i) {
                    func(i, j, k);
                }
            }
        }
    }

    template <typename T>
    template <typename Callback>
    void ConstArrayAccessor<T, 3>::parallelForEachIndex(Callback func) const {
        parallelFor(
            kZeroSize, _size.x, kZeroSize, _size.y, kZeroSize, _size.z, func);
    }

    template <typename T>
    size_t ConstArrayAccessor<T, 3>::index(const Point3UI& pt) const {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y&& pt.z < _size.z);
        return pt.x + _size.x * (pt.y + _size.y * pt.z);
    }

    template <typename T>
    size_t ConstArrayAccessor<T, 3>::index(size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return i + _size.x * (j + _size.y * k);
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::operator[](size_t i) const {
        return _data[i];
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::operator()(
        size_t i, size_t j, size_t k) const {
        JET_ASSERT(i < _size.x&& j < _size.y&& k < _size.z);
        return _data[i + _size.x * (j + _size.y * k)];
    }

    template <typename T>
    const T& ConstArrayAccessor<T, 3>::operator()(const Point3UI& pt) const {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y&& pt.z < _size.z);
        return _data[pt.x + _size.x * (pt.y + _size.y * pt.z)];
    }

}  // namespace jet


#endif  // INCLUDE_JET_ARRAY_ACCESSOR3_H_