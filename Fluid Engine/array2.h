#ifndef INCLUDE_JET_ARRAY2_H_
#define INCLUDE_JET_ARRAY2_H_

#include "array.h"
#include "array_accessor2.h"
#include "size2.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <utility>  // just make cpplint happy..
#include <vector>

#include "macros.h"
#include "parallel.h"
#include <algorithm>

namespace jet {

    //
    // 2-D array class.
    //
    // This class represents 2-D array data structure. Internally, the 2-D data is
    // mapped to a linear array such that (i, j) element is actually stroed at
    // (i + (width * j))th element of the linear array.
    // 
    // T - Type to store in the array.
    //
    template <typename T>
    class Array<T, 2> final {
    public:
        typedef std::vector<T> ContainerType;
        typedef typename ContainerType::iterator Iterator;
        typedef typename ContainerType::const_iterator ConstIterator;

        // Constructs zero-sized 2-D array.
        Array();

        // Constructs 2-D array with given size and fill it with initVal.
        // size Initial size of the array.
        // initVal Initial value of each array element.
        explicit Array(const Size2& size, const T& initVal = T());

        // Constructs 2-D array with size width x height and fill it with
        // initVal.
        // width Initial width of the array.
        // height Initial height of the array.
        // initVal Initial value of each array element.
        Array(size_t width, size_t height, const T& initVal = T());

        //
        // Constructs 2-D array with given initializer list lst.
        //
        Array(const std::initializer_list<std::initializer_list<T>>& lst);

        // Copy constructor.
        Array(const Array& other);

        // Move constructor.
        Array(Array&& other);

        // Sets entire array with given value.
        void set(const T& value);

        // Copies given array other to this array.
        void set(const Array& other);

        //
        // Copies given initializer list lst to this array.
        //
        void set(const std::initializer_list<std::initializer_list<T>>& lst);

        // Clears the array and resizes to zero.
        void clear();

        // Resizes the array with size and fill the new element with initVal.
        void resize(const Size2& size, const T& initVal = T());

        // Resizes the array with size width x height and fill the new
        // element with initVal.
        void resize(size_t width, size_t height, const T& initVal = T());

        //
        // Returns the reference to the i-th element.
        //
        // This function returns the reference to the i-th element of the array
        // where i is the index of linearly mapped elements such that
        // i = x + (width * y) (x and y are the 2-D coordinates of the element).
        //
        T& at(size_t i);

        //
        // Returns the const reference to the i-th element.
        //
        // This function returns the const reference to the i-th element of the
        // array where i is the index of linearly mapped elements such that
        // i = x + (width * y) (x and y are the 2-D coordinates of the element).
        //
        const T& at(size_t i) const;

        // Returns the reference to the element at (pt.x, pt.y).
        T& at(const Point2UI& pt);

        // Returns the const reference to the element at (pt.x, pt.y).
        const T& at(const Point2UI& pt) const;

        // Returns the reference to the element at (i, j).
        T& at(size_t i, size_t j);

        // Returns the const reference to the element at (i, j).
        const T& at(size_t i, size_t j) const;

        // Returns the size of the array.
        Size2 size() const;

        // Returns the width of the array.
        size_t width() const;

        // Returns the height of the array.
        size_t height() const;

        // Returns the raw pointer to the array data.
        T* data();

        // Returns the const raw pointer to the array data.
        const T* const data() const;

        // Returns the begin iterator of the array.
        Iterator begin();

        // Returns the begin const iterator of the array.
        ConstIterator begin() const;

        // Returns the end iterator of the array.
        Iterator end();

        // Returns the end const iterator of the array.
        ConstIterator end() const;

        // Returns the array accessor.
        ArrayAccessor2<T> accessor();

        // Returns the const array accessor.
        ConstArrayAccessor2<T> constAccessor() const;

        // Swaps the content of the array with other array.
        void swap(Array& other);

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
        // function func. The callback function takes two parameters which are
        // the (i, j) indices of the array.
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
        // function takes two parameters which are the (i, j) indices of the array.
        // The order of execution will be non-deterministic since it runs in
        // parallel.
        //
        template <typename Callback>
        void parallelForEachIndex(Callback func) const;

        //
        // Returns the reference to the i-th element.
        //
        // This function returns the reference to the i-th element of the array
        // where i is the index of linearly mapped elements such that
        // i = x + (width * y) (x and y are the 2-D coordinates of the element).
        //
        T& operator[](size_t i);

        //
        // Returns the const reference to the i-th element.
        //
        // This function returns the const reference to the i-th element of the
        // array where i is the index of linearly mapped elements such that
        // i = x + (width * y) (x and y are the 2-D coordinates of the element).
        //
        const T& operator[](size_t i) const;

        // Returns the reference to the element at (pt.x, pt.y).
        T& operator()(const Point2UI& pt);

        // Returns the const reference to the element at (pt.x, pt.y).
        const T& operator()(const Point2UI& pt) const;

        // Returns the reference to the element at (i, j).
        T& operator()(size_t i, size_t j);

        // Returns the const reference to the element at (i, j).
        const T& operator()(size_t i, size_t j) const;

        // Sets entire array with given value.
        Array& operator=(const T& other);

        // Copies given array other to this array.
        Array& operator=(const Array& other);

        // Move assignment.
        Array& operator=(Array&& other);

        //
        // Copies given initializer list lst to this array.
        //
        // lst Initializer list that should be copy to the new array.
        //
        Array& operator=(
            const std::initializer_list<std::initializer_list<T>>& lst);

        // Casts to array accessor.
        operator ArrayAccessor2<T>();

        // Casts to const array accessor.
        operator ConstArrayAccessor2<T>() const;

    private:
        Size2 _size;
        std::vector<T> _data;
    };

    // Type alias for 2-D array.
    template <typename T>
    using Array2 = Array<T, 2>;

    template <typename T>
    Array<T, 2>::Array() {}

    template <typename T>
    Array<T, 2>::Array(const Size2& size, const T& initVal) {
        resize(size, initVal);
    }

    template <typename T>
    Array<T, 2>::Array(size_t width, size_t height, const T& initVal) {
        resize(width, height, initVal);
    }

    template <typename T>
    Array<T, 2>::Array(const std::initializer_list<std::initializer_list<T>>& lst) {
        set(lst);
    }

    template <typename T>
    Array<T, 2>::Array(const Array& other) {
        set(other);
    }

    template <typename T>
    Array<T, 2>::Array(Array&& other) {
        (*this) = std::move(other);
    }

    template <typename T>
    void Array<T, 2>::set(const T& value) {
        for (auto& v : _data) {
            v = value;
        }
    }

    template <typename T>
    void Array<T, 2>::set(const Array& other) {
        _data.resize(other._data.size());
        std::copy(other._data.begin(), other._data.end(), _data.begin());
        _size = other._size;
    }

    template <typename T>
    void Array<T, 2>::set(
        const std::initializer_list<std::initializer_list<T>>& lst) {
        size_t height = lst.size();
        size_t width = (height > 0) ? lst.begin()->size() : 0;
        resize(Size2(width, height));
        auto rowIter = lst.begin();
        for (size_t j = 0; j < height; ++j) {
            JET_ASSERT(width == rowIter->size());
            auto colIter = rowIter->begin();
            for (size_t i = 0; i < width; ++i) {
                (*this)(i, j) = *colIter;
                ++colIter;
            }
            ++rowIter;
        }
    }

    template <typename T>
    void Array<T, 2>::clear() {
        _data.clear();
        _size = Size2(0, 0);
    }

    template <typename T>
    void Array<T, 2>::resize(const Size2& size, const T& initVal) {
        Array grid;
        grid._data.resize(size.x * size.y, initVal);
        grid._size = size;
        size_t iMin = std::min(size.x, _size.x);
        size_t jMin = std::min(size.y, _size.y);
        for (size_t j = 0; j < jMin; ++j) {
            for (size_t i = 0; i < iMin; ++i) {
                grid(i, j) = at(i, j);
            }
        }

        swap(grid);
    }

    template <typename T>
    void Array<T, 2>::resize(size_t width, size_t height, const T& initVal) {
        resize(Size2(width, height), initVal);
    }

    template <typename T>
    T& Array<T, 2>::at(size_t i) {
        JET_ASSERT(i < _size.x* _size.y);
        return _data[i];
    }

    template <typename T>
    const T& Array<T, 2>::at(size_t i) const {
        JET_ASSERT(i < _size.x* _size.y);
        return _data[i];
    }

    template <typename T>
    T& Array<T, 2>::at(const Point2UI& pt) {
        return at(pt.x, pt.y);
    }

    template <typename T>
    const T& Array<T, 2>::at(const Point2UI& pt) const {
        return at(pt.x, pt.y);
    }

    template <typename T>
    T& Array<T, 2>::at(size_t i, size_t j) {
        JET_ASSERT(i < _size.x&& j < _size.y);
        return _data[i + _size.x * j];
    }

    template <typename T>
    const T& Array<T, 2>::at(size_t i, size_t j) const {
        JET_ASSERT(i < _size.x&& j < _size.y);
        return _data[i + _size.x * j];
    }

    template <typename T>
    Size2 Array<T, 2>::size() const {
        return _size;
    }

    template <typename T>
    size_t Array<T, 2>::width() const {
        return _size.x;
    }

    template <typename T>
    size_t Array<T, 2>::height() const {
        return _size.y;
    }

    template <typename T>
    T* Array<T, 2>::data() {
        return _data.data();
    }

    template <typename T>
    const T* const Array<T, 2>::data() const {
        return _data.data();
    }

    template <typename T>
    typename Array<T, 2>::ContainerType::iterator Array<T, 2>::begin() {
        return _data.begin();
    }

    template <typename T>
    typename Array<T, 2>::ContainerType::const_iterator Array<T, 2>::begin() const {
        return _data.cbegin();
    }

    template <typename T>
    typename Array<T, 2>::ContainerType::iterator Array<T, 2>::end() {
        return _data.end();
    }

    template <typename T>
    typename Array<T, 2>::ContainerType::const_iterator Array<T, 2>::end() const {
        return _data.cend();
    }

    template <typename T>
    ArrayAccessor2<T> Array<T, 2>::accessor() {
        return ArrayAccessor2<T>(size(), data());
    }

    template <typename T>
    ConstArrayAccessor2<T> Array<T, 2>::constAccessor() const {
        return ConstArrayAccessor2<T>(size(), data());
    }

    template <typename T>
    void Array<T, 2>::swap(Array& other) {
        std::swap(other._data, _data);
        std::swap(other._size, _size);
    }

    template <typename T>
    template <typename Callback>
    void Array<T, 2>::forEach(Callback func) const {
        constAccessor().forEach(func);
    }

    template <typename T>
    template <typename Callback>
    void Array<T, 2>::forEachIndex(Callback func) const {
        constAccessor().forEachIndex(func);
    }

    template <typename T>
    template <typename Callback>
    void Array<T, 2>::parallelForEach(Callback func) {
        accessor().parallelForEach(func);
    }

    template <typename T>
    template <typename Callback>
    void Array<T, 2>::parallelForEachIndex(Callback func) const {
        constAccessor().parallelForEachIndex(func);
    }

    template <typename T>
    T& Array<T, 2>::operator[](size_t i) {
        return _data[i];
    }

    template <typename T>
    const T& Array<T, 2>::operator[](size_t i) const {
        return _data[i];
    }

    template <typename T>
    T& Array<T, 2>::operator()(size_t i, size_t j) {
        JET_ASSERT(i < _size.x&& j < _size.y);
        return _data[i + _size.x * j];
    }

    template <typename T>
    const T& Array<T, 2>::operator()(size_t i, size_t j) const {
        JET_ASSERT(i < _size.x&& j < _size.y);
        return _data[i + _size.x * j];
    }

    template <typename T>
    T& Array<T, 2>::operator()(const Point2UI& pt) {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y);
        return _data[pt.x + _size.x * pt.y];
    }

    template <typename T>
    const T& Array<T, 2>::operator()(const Point2UI& pt) const {
        JET_ASSERT(pt.x < _size.x&& pt.y < _size.y);
        return _data[pt.x + _size.x * pt.y];
    }

    template <typename T>
    Array<T, 2>& Array<T, 2>::operator=(const T& value) {
        set(value);
        return *this;
    }

    template <typename T>
    Array<T, 2>& Array<T, 2>::operator=(const Array& other) {
        set(other);
        return *this;
    }

    template <typename T>
    Array<T, 2>& Array<T, 2>::operator=(Array&& other) {
        _data = std::move(other._data);
        _size = other._size;
        other._size = Size2();
        return *this;
    }

    template <typename T>
    Array<T, 2>& Array<T, 2>::operator=(
        const std::initializer_list<std::initializer_list<T>>& lst) {
        set(lst);
        return *this;
    }

    template <typename T>
    Array<T, 2>::operator ArrayAccessor2<T>() {
        return accessor();
    }

    template <typename T>
    Array<T, 2>::operator ConstArrayAccessor2<T>() const {
        return constAccessor();
    }

}  // namespace jet


#endif  // INCLUDE_JET_ARRAY2_H_