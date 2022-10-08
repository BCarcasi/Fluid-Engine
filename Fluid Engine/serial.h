#ifndef INCLUDE_JET_SERIAL_H_
#define INCLUDE_JET_SERIAL_H_


#include "macros.h"
#include <algorithm>
#include <functional>
#include <vector>

namespace jet {

    //
    //      Fills from begin to end with value.
    //
    // This function fills a container specified by begin and end iterators with
    // single thread. The order of the filling is deterministic.
    //
    //  begin          The begin iterator of a container.
    //  end            The end iterator of a container.
    //  value          The value to fill a container.
    //
    //     RandomIterator Random iterator type.
    //     T              Value type of a container.
    //
    template <typename RandomIterator, typename T>
    void serialFill(
        const RandomIterator& begin, const RandomIterator& end, const T& value);

    //
    //      Makes a for-loop from beginIndex to endIndex.
    //
    // This function makes a for-loop specified by begin and end indices with
    // single thread. The order of the visit is deterministic.
    //
    //  beginIndex The begin index.
    //  endIndex   The end index.
    //  function   The function to call for each index.
    //
    //     IndexType  Index type.
    //     Function   Function type.
    //
    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndex, IndexType endIndex, const Function& function);

    //
    //      Makes a 2D nested for-loop.
    //
    // This function makes a 2D nested for-loop specified by begin and end indices
    // for each dimension. X will be the inner-most loop while Y is the outer-most.
    // The order of the visit is deterministic.
    //
    //  beginIndexX The begin index in X dimension.
    //  endIndexX   The end index in X dimension.
    //  beginIndexY The begin index in Y dimension.
    //  endIndexY   The end index in Y dimension.
    //  function    The function to call for each index (i, j).
    //
    //     IndexType  Index type.
    //     Function   Function type.
    //
    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndexX,
        IndexType endIndexX,
        IndexType beginIndexY,
        IndexType endIndexY,
        const Function& function);

    //
    //      Makes a 3D nested for-loop.
    //
    // This function makes a 3D nested for-loop specified by begin and end indices
    // for each dimension. X will be the inner-most loop while Z is the outer-most.
    // The order of the visit is deterministic.
    //
    //  beginIndexX The begin index in X dimension.
    //  endIndexX   The end index in X dimension.
    //  beginIndexY The begin index in Y dimension.
    //  endIndexY   The end index in Y dimension.
    //  beginIndexZ The begin index in Z dimension.
    //  endIndexZ   The end index in Z dimension.
    //  function    The function to call for each index (i, j, k).
    //
    //     IndexType   Index type.
    //     Function    Function type.
    //
    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndexX,
        IndexType endIndexX,
        IndexType beginIndexY,
        IndexType endIndexY,
        IndexType beginIndexZ,
        IndexType endIndexZ,
        const Function& function);

    //
    //      Sorts a container.
    //
    // This function sorts a container specified by begin and end iterators.
    //
    //  begin          The begin random access iterator.
    //  end            The end random access iterator.
    //
    //     RandomIterator Iterator type.
    //
    template<typename RandomIterator>
    void serialSort(RandomIterator begin, RandomIterator end);

    //
    //      Sorts a container with a custom compare function.
    //
    // This function sorts a container specified by begin and end iterators. It
    // takes extra compare function which returns true if the first argument is
    // less than the second argument.
    //
    //  begin           The begin random access iterator.
    //  end             The end random access iterator.
    //  compare         The compare function.
    //
    //     RandomIterator  Iterator type.
    //     CompareFunction Compare function type.
    //
    template<typename RandomIterator, typename SortingFunction>
    void serialSort(
        RandomIterator begin,
        RandomIterator end,
        const SortingFunction& sortingFunction);

    template <typename RandomIterator, typename T>
    void serialFill(
        const RandomIterator& begin,
        const RandomIterator& end,
        const T& value) {
        size_t size = static_cast<size_t>(end - begin);
        serialFor(size_t(0), size, [begin, value](size_t i) {
            begin[i] = value;
            });
    }

    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndex,
        IndexType endIndex,
        const Function& function) {
        for (IndexType i = beginIndex; i < endIndex; ++i) {
            function(i);
        }
    }

    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndexX,
        IndexType endIndexX,
        IndexType beginIndexY,
        IndexType endIndexY,
        const Function& function) {
        for (IndexType j = beginIndexY; j < endIndexY; ++j) {
            for (IndexType i = beginIndexX; i < endIndexX; ++i) {
                function(i, j);
            }
        }
    }

    template <typename IndexType, typename Function>
    void serialFor(
        IndexType beginIndexX,
        IndexType endIndexX,
        IndexType beginIndexY,
        IndexType endIndexY,
        IndexType beginIndexZ,
        IndexType endIndexZ,
        const Function& function) {
        for (IndexType k = beginIndexZ; k < endIndexZ; ++k) {
            for (IndexType j = beginIndexY; j < endIndexY; ++j) {
                for (IndexType i = beginIndexX; i < endIndexX; ++i) {
                    function(i, j, k);
                }
            }
        }
    }

    template<typename RandomIterator, typename SortingFunction>
    void serialSort(
        RandomIterator begin,
        RandomIterator end,
        const SortingFunction& sortingFunction) {
        std::sort(begin, end, sortingFunction);
    }

    template<typename RandomIterator>
    void serialSort(RandomIterator begin, RandomIterator end) {
        serialSort(begin, end, std::less<typename RandomIterator::value_type>());
    }

}  // namespace jet


#endif  // INCLUDE_JET_SERIAL_H_