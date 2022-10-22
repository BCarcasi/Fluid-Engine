#ifndef INCLUDE_JET_POINT_PARALLEL_HASH_GRID_SEARCHER3_H_
#define INCLUDE_JET_POINT_PARALLEL_HASH_GRID_SEARCHER3_H_

#include "point_neighbor_searcher3.h"
#include "point3.h"
#include "size3.h"
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"
#include "fbs_helpers.h"
#include "point_parallel_hash_grid_searcher3_generated.h"

#include "constants.h"
#include "parallel.h"

#include <algorithm>
#include <vector>

namespace jet {

    //!
    //! \brief Parallel version of hash grid-based 3-D point searcher.
    //!
    //! This class implements parallel version of 3-D point searcher by using hash
    //! grid for its internal acceleration data structure. Each point is recorded to
    //! its corresponding bucket where the hashing function is 3-D grid mapping.
    //!
    class PointParallelHashGridSearcher3 final : public PointNeighborSearcher3 {
    public:
        JET_NEIGHBOR_SEARCHER3_TYPE_NAME(PointParallelHashGridSearcher3)

            class Builder;

        //!
        //! \brief      Constructs hash grid with given resolution and grid spacing.
        //!
        //! This constructor takes hash grid resolution and its grid spacing as
        //! its input parameters. The grid spacing must be 2x or greater than
        //! search radius.
        //!
        //! \param[in]  resolution  The resolution.
        //! \param[in]  gridSpacing The grid spacing.
        //!
        PointParallelHashGridSearcher3(
            const Size3& resolution, double gridSpacing);

        //!
        //! \brief      Constructs hash grid with given resolution and grid spacing.
        //!
        //! This constructor takes hash grid resolution and its grid spacing as
        //! its input parameters. The grid spacing must be 2x or greater than
        //! search radius.
        //!
        //! \param[in]  resolutionX The resolution x.
        //! \param[in]  resolutionY The resolution y.
        //! \param[in]  resolutionZ The resolution z.
        //! \param[in]  gridSpacing The grid spacing.
        //!
        PointParallelHashGridSearcher3(
            size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ,
            double gridSpacing);

        //! Copy constructor
        PointParallelHashGridSearcher3(const PointParallelHashGridSearcher3& other);

        //!
        //! \brief Builds internal acceleration structure for given points list.
        //!
        //! This function builds the hash grid for given points in parallel.
        //!
        //! \param[in]  points The points to be added.
        //!
        void build(const ConstArrayAccessor1<Vector3D>& points) override;

        //!
        //! Invokes the callback function for each nearby point around the origin
        //! within given radius.
        //!
        //! \param[in]  origin   The origin position.
        //! \param[in]  radius   The search radius.
        //! \param[in]  callback The callback function.
        //!
        void forEachNearbyPoint(
            const Vector3D& origin,
            double radius,
            const ForEachNearbyPointFunc& callback) const override;

        //!
        //! Returns true if there are any nearby points for given origin within
        //! radius.
        //!
        //! \param[in]  origin The origin.
        //! \param[in]  radius The radius.
        //!
        //! \return     True if has nearby point, false otherwise.
        //!
        bool hasNearbyPoint(
            const Vector3D& origin, double radius) const override;

        //!
        //! \brief      Returns the hash key list.
        //!
        //! The hash key list maps sorted point index i to its hash key value.
        //! The sorting order is based on the key value itself.
        //!
        //! \return     The hash key list.
        //!
        const std::vector<size_t>& keys() const;

        //!
        //! \brief      Returns the start index table.
        //!
        //! The start index table maps the hash grid bucket index to starting index
        //! of the sorted point list. Assume the hash key list looks like:
        //!
        //! \code
        //! [5|8|8|10|10|10]
        //! \endcode
        //!
        //! Then startIndexTable and endIndexTable should be like:
        //!
        //! \code
        //! [.....|0|...|1|..|3|..]
        //! [.....|1|...|3|..|6|..]
        //!       ^5    ^8   ^10
        //! \endcode
        //!
        //! So that endIndexTable[i] - startIndexTable[i] is the number points
        //! in i-th table bucket.
        //!
        //! \return     The start index table.
        //!
        const std::vector<size_t>& startIndexTable() const;

        //!
        //! \brief      Returns the end index table.
        //!
        //! The end index table maps the hash grid bucket index to starting index
        //! of the sorted point list. Assume the hash key list looks like:
        //!
        //! \code
        //! [5|8|8|10|10|10]
        //! \endcode
        //!
        //! Then startIndexTable and endIndexTable should be like:
        //!
        //! \code
        //! [.....|0|...|1|..|3|..]
        //! [.....|1|...|3|..|6|..]
        //!       ^5    ^8   ^10
        //! \endcode
        //!
        //! So that endIndexTable[i] - startIndexTable[i] is the number points
        //! in i-th table bucket.
        //!
        //! \return     The end index table.
        //!
        const std::vector<size_t>& endIndexTable() const;

        //!
        //! \brief      Returns the sorted indices of the points.
        //!
        //! When the hash grid is built, it sorts the points in hash key order. But
        //! rather than sorting the original points, this class keeps the shuffled
        //! indices of the points. The list this function returns maps sorted index
        //! i to original index j.
        //!
        //! \return     The sorted indices of the points.
        //!
        const std::vector<size_t>& sortedIndices() const;

        //!
        //! Returns the hash value for given 3-D bucket index.
        //!
        //! \param[in]  bucketIndex The bucket index.
        //!
        //! \return     The hash key from bucket index.
        //!
        size_t getHashKeyFromBucketIndex(const Point3I& bucketIndex) const;

        //!
        //! Gets the bucket index from a point.
        //!
        //! \param[in]  position The position of the point.
        //!
        //! \return     The bucket index.
        //!
        Point3I getBucketIndex(const Vector3D& position) const;

        //!
        //! \brief      Creates a new instance of the object with same properties
        //!             than original.
        //!
        //! \return     Copy of this object.
        //!
        PointNeighborSearcher3Ptr clone() const override;

        //! Assignment operator.
        PointParallelHashGridSearcher3& operator=(
            const PointParallelHashGridSearcher3& other);

        //! Copy from the other instance.
        void set(const PointParallelHashGridSearcher3& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointParallelHashGridSearcher3.
        static Builder builder();

    private:
        double _gridSpacing = 1.0;
        Point3I _resolution = Point3I(1, 1, 1);
        std::vector<Vector3D> _points;
        std::vector<size_t> _keys;
        std::vector<size_t> _startIndexTable;
        std::vector<size_t> _endIndexTable;
        std::vector<size_t> _sortedIndices;

        size_t getHashKeyFromPosition(const Vector3D& position) const;

        void getNearbyKeys(const Vector3D& position, size_t* bucketIndices) const;
    };

    //! Shared pointer for the PointParallelHashGridSearcher3 type.
    typedef std::shared_ptr<PointParallelHashGridSearcher3>
        PointParallelHashGridSearcher3Ptr;

    //!
    //! \brief Front-end to create PointParallelHashGridSearcher3 objects step by
    //!        step.
    //!
    class PointParallelHashGridSearcher3::Builder final
        : public PointNeighborSearcherBuilder3 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size3& resolution);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(double gridSpacing);

        //! Builds PointParallelHashGridSearcher3 instance.
        PointParallelHashGridSearcher3 build() const;

        //! Builds shared pointer of PointParallelHashGridSearcher3 instance.
        PointParallelHashGridSearcher3Ptr makeShared() const;

        //! Returns shared pointer of PointNeighborSearcher3 type.
        PointNeighborSearcher3Ptr buildPointNeighborSearcher() const override;

    private:
        Size3 _resolution{ 64, 64, 64 };
        double _gridSpacing = 1.0;
    };

    PointParallelHashGridSearcher3::PointParallelHashGridSearcher3(
        const Size3& resolution,
        double gridSpacing) :
        PointParallelHashGridSearcher3(
            resolution.x, resolution.y, resolution.z, gridSpacing) {
    }

    PointParallelHashGridSearcher3::PointParallelHashGridSearcher3(
        size_t resolutionX,
        size_t resolutionY,
        size_t resolutionZ,
        double gridSpacing) :
        _gridSpacing(gridSpacing) {
        _resolution.x = std::max(static_cast<ssize_t>(resolutionX), kOneSSize);
        _resolution.y = std::max(static_cast<ssize_t>(resolutionY), kOneSSize);
        _resolution.z = std::max(static_cast<ssize_t>(resolutionZ), kOneSSize);

        _startIndexTable.resize(
            _resolution.x * _resolution.y * _resolution.z, kMaxSize);
        _endIndexTable.resize(
            _resolution.x * _resolution.y * _resolution.z, kMaxSize);
    }

    PointParallelHashGridSearcher3::PointParallelHashGridSearcher3(
        const PointParallelHashGridSearcher3& other) {
        set(other);
    }

    void PointParallelHashGridSearcher3::build(
        const ConstArrayAccessor1<Vector3D>& points) {
        _points.clear();
        _keys.clear();
        _startIndexTable.clear();
        _endIndexTable.clear();
        _sortedIndices.clear();

        // Allocate memory chuncks
        size_t numberOfPoints = points.size();
        std::vector<size_t> tempKeys(numberOfPoints);
        _startIndexTable.resize(_resolution.x * _resolution.y * _resolution.z);
        _endIndexTable.resize(_resolution.x * _resolution.y * _resolution.z);
        parallelFill(_startIndexTable.begin(), _startIndexTable.end(), kMaxSize);
        parallelFill(_endIndexTable.begin(), _endIndexTable.end(), kMaxSize);
        _keys.resize(numberOfPoints);
        _sortedIndices.resize(numberOfPoints);
        _points.resize(numberOfPoints);

        if (numberOfPoints == 0) {
            return;
        }

        // Initialize indices array and generate hash key for each point
        parallelFor(
            kZeroSize,
            numberOfPoints,
            [&](size_t i) {
                _sortedIndices[i] = i;
                _points[i] = points[i];
                tempKeys[i] = getHashKeyFromPosition(points[i]);
            });

        // Sort indices based on hash key
        parallelSort(
            _sortedIndices.begin(),
            _sortedIndices.end(),
            [&tempKeys](size_t indexA, size_t indexB) {
                return tempKeys[indexA] < tempKeys[indexB];
            });

        // Re-order point and key arrays
        parallelFor(
            kZeroSize,
            numberOfPoints,
            [&](size_t i) {
                _points[i] = points[_sortedIndices[i]];
                _keys[i] = tempKeys[_sortedIndices[i]];
            });

        // Now _points and _keys are sorted by points' hash key values.
        // Let's fill in start/end index table with _keys.

        // Assume that _keys array looks like:
        // [5|8|8|10|10|10]
        // Then _startIndexTable and _endIndexTable should be like:
        // [.....|0|...|1|..|3|..]
        // [.....|1|...|3|..|6|..]
        //       ^5    ^8   ^10
        // So that _endIndexTable[i] - _startIndexTable[i] is the number points
        // in i-th table bucket.

        _startIndexTable[_keys[0]] = 0;
        _endIndexTable[_keys[numberOfPoints - 1]] = numberOfPoints;

        parallelFor(
            (size_t)1,
            numberOfPoints,
            [&](size_t i) {
                if (_keys[i] > _keys[i - 1]) {
                    _startIndexTable[_keys[i]] = i;
                    _endIndexTable[_keys[i - 1]] = i;
                }
            });

        size_t sumNumberOfPointsPerBucket = 0;
        size_t maxNumberOfPointsPerBucket = 0;
        size_t numberOfNonEmptyBucket = 0;
        for (size_t i = 0; i < _startIndexTable.size(); ++i) {
            if (_startIndexTable[i] != kMaxSize) {
                size_t numberOfPointsInBucket
                    = _endIndexTable[i] - _startIndexTable[i];
                sumNumberOfPointsPerBucket += numberOfPointsInBucket;
                maxNumberOfPointsPerBucket =
                    std::max(maxNumberOfPointsPerBucket, numberOfPointsInBucket);
                ++numberOfNonEmptyBucket;
            }
        }

        JET_INFO << "Average number of points per non-empty bucket: "
            << static_cast<float>(sumNumberOfPointsPerBucket)
            / static_cast<float>(numberOfNonEmptyBucket);
        JET_INFO << "Max number of points per bucket: "
            << maxNumberOfPointsPerBucket;
    }

    void PointParallelHashGridSearcher3::forEachNearbyPoint(
        const Vector3D& origin,
        double radius,
        const ForEachNearbyPointFunc& callback) const {
        size_t nearbyKeys[8];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 8; i++) {
            size_t nearbyKey = nearbyKeys[i];
            size_t start = _startIndexTable[nearbyKey];
            size_t end = _endIndexTable[nearbyKey];

            // Empty bucket -- continue to next bucket
            if (start == kMaxSize) {
                continue;
            }

            for (size_t j = start; j < end; ++j) {
                Vector3D direction = _points[j] - origin;
                double distanceSquared = direction.lengthSquared();
                if (distanceSquared <= queryRadiusSquared) {
                    double distance = 0.0;
                    if (distanceSquared > 0) {
                        distance = std::sqrt(distanceSquared);
                        direction /= distance;
                    }

                    callback(_sortedIndices[j], _points[j]);
                }
            }
        }
    }

    bool PointParallelHashGridSearcher3::hasNearbyPoint(
        const Vector3D& origin,
        double radius) const {
        size_t nearbyKeys[8];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 8; i++) {
            size_t nearbyKey = nearbyKeys[i];
            size_t start = _startIndexTable[nearbyKey];
            size_t end = _endIndexTable[nearbyKey];

            // Empty bucket -- continue to next bucket
            if (start == kMaxSize) {
                continue;
            }

            for (size_t j = start; j < end; ++j) {
                Vector3D direction = _points[j] - origin;
                double distanceSquared = direction.lengthSquared();
                if (distanceSquared <= queryRadiusSquared) {
                    return true;
                }
            }
        }

        return false;
    }

    const std::vector<size_t>& PointParallelHashGridSearcher3::keys() const {
        return _keys;
    }

    const std::vector<size_t>&
        PointParallelHashGridSearcher3::startIndexTable() const {
        return _startIndexTable;
    }

    const std::vector<size_t>&
        PointParallelHashGridSearcher3::endIndexTable() const {
        return _endIndexTable;
    }

    const std::vector<size_t>&
        PointParallelHashGridSearcher3::sortedIndices() const {
        return _sortedIndices;
    }

    Point3I PointParallelHashGridSearcher3::getBucketIndex(
        const Vector3D& position) const {
        Point3I bucketIndex;
        bucketIndex.x = static_cast<ssize_t>(
            std::floor(position.x / _gridSpacing));
        bucketIndex.y = static_cast<ssize_t>(
            std::floor(position.y / _gridSpacing));
        bucketIndex.z = static_cast<ssize_t>(
            std::floor(position.z / _gridSpacing));
        return bucketIndex;
    }

    size_t PointParallelHashGridSearcher3::getHashKeyFromPosition(
        const Vector3D& position) const {
        Point3I bucketIndex = getBucketIndex(position);

        return getHashKeyFromBucketIndex(bucketIndex);
    }

    size_t PointParallelHashGridSearcher3::getHashKeyFromBucketIndex(
        const Point3I& bucketIndex) const {
        Point3I wrappedIndex = bucketIndex;
        wrappedIndex.x = bucketIndex.x % _resolution.x;
        wrappedIndex.y = bucketIndex.y % _resolution.y;
        wrappedIndex.z = bucketIndex.z % _resolution.z;
        if (wrappedIndex.x < 0) {
            wrappedIndex.x += _resolution.x;
        }
        if (wrappedIndex.y < 0) {
            wrappedIndex.y += _resolution.y;
        }
        if (wrappedIndex.z < 0) {
            wrappedIndex.z += _resolution.z;
        }
        return static_cast<size_t>(
            (wrappedIndex.z * _resolution.y + wrappedIndex.y) * _resolution.x
            + wrappedIndex.x);
    }

    void PointParallelHashGridSearcher3::getNearbyKeys(
        const Vector3D& position,
        size_t* nearbyKeys) const {
        Point3I originIndex = getBucketIndex(position), nearbyBucketIndices[8];

        for (int i = 0; i < 8; i++) {
            nearbyBucketIndices[i] = originIndex;
        }

        if ((originIndex.x + 0.5f) * _gridSpacing <= position.x) {
            nearbyBucketIndices[4].x += 1;
            nearbyBucketIndices[5].x += 1;
            nearbyBucketIndices[6].x += 1;
            nearbyBucketIndices[7].x += 1;
        }
        else {
            nearbyBucketIndices[4].x -= 1;
            nearbyBucketIndices[5].x -= 1;
            nearbyBucketIndices[6].x -= 1;
            nearbyBucketIndices[7].x -= 1;
        }

        if ((originIndex.y + 0.5f) * _gridSpacing <= position.y) {
            nearbyBucketIndices[2].y += 1;
            nearbyBucketIndices[3].y += 1;
            nearbyBucketIndices[6].y += 1;
            nearbyBucketIndices[7].y += 1;
        }
        else {
            nearbyBucketIndices[2].y -= 1;
            nearbyBucketIndices[3].y -= 1;
            nearbyBucketIndices[6].y -= 1;
            nearbyBucketIndices[7].y -= 1;
        }

        if ((originIndex.z + 0.5f) * _gridSpacing <= position.z) {
            nearbyBucketIndices[1].z += 1;
            nearbyBucketIndices[3].z += 1;
            nearbyBucketIndices[5].z += 1;
            nearbyBucketIndices[7].z += 1;
        }
        else {
            nearbyBucketIndices[1].z -= 1;
            nearbyBucketIndices[3].z -= 1;
            nearbyBucketIndices[5].z -= 1;
            nearbyBucketIndices[7].z -= 1;
        }

        for (int i = 0; i < 8; i++) {
            nearbyKeys[i] = getHashKeyFromBucketIndex(nearbyBucketIndices[i]);
        }
    }

    PointNeighborSearcher3Ptr PointParallelHashGridSearcher3::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointParallelHashGridSearcher3);
    }

    PointParallelHashGridSearcher3&
        PointParallelHashGridSearcher3::operator=(
            const PointParallelHashGridSearcher3& other) {
        set(other);
        return *this;
    }

    void PointParallelHashGridSearcher3::set(
        const PointParallelHashGridSearcher3& other) {
        _gridSpacing = other._gridSpacing;
        _resolution = other._resolution;
        _points = other._points;
        _keys = other._keys;
        _startIndexTable = other._startIndexTable;
        _endIndexTable = other._endIndexTable;
        _sortedIndices = other._sortedIndices;
    }

    void PointParallelHashGridSearcher3::serialize(
        std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        // Copy simple data
        auto fbsResolution
            = fbs::Size3(_resolution.x, _resolution.y, _resolution.z);

        // Copy points
        std::vector<fbs::Vector3D> points;
        for (const auto& pt : _points) {
            points.push_back(jetToFbs(pt));
        }

        auto fbsPoints
            = builder.CreateVectorOfStructs(points.data(), points.size());

        // Copy key/tables
        std::vector<uint64_t> keys(_keys.begin(), _keys.end());
        std::vector<uint64_t> startIndexTable(
            _startIndexTable.begin(), _startIndexTable.end());
        std::vector<uint64_t> endIndexTable(
            _endIndexTable.begin(), _endIndexTable.end());
        std::vector<uint64_t> sortedIndices(
            _sortedIndices.begin(), _sortedIndices.end());

        auto fbsKeys = builder.CreateVector(keys.data(), keys.size());
        auto fbsStartIndexTable
            = builder.CreateVector(startIndexTable.data(), startIndexTable.size());
        auto fbsEndIndexTable
            = builder.CreateVector(endIndexTable.data(), endIndexTable.size());
        auto fbsSortedIndices
            = builder.CreateVector(sortedIndices.data(), sortedIndices.size());

        // Copy the searcher
        auto fbsSearcher = fbs::CreatePointParallelHashGridSearcher3(
            builder,
            _gridSpacing,
            &fbsResolution,
            fbsPoints,
            fbsKeys,
            fbsStartIndexTable,
            fbsEndIndexTable,
            fbsSortedIndices);

        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointParallelHashGridSearcher3::deserialize(
        const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointParallelHashGridSearcher3(buffer.data());

        // Copy simple data
        auto res = fbsToJet(*fbsSearcher->resolution());
        _resolution.set({ res.x, res.y, res.z });
        _gridSpacing = fbsSearcher->gridSpacing();

        // Copy points
        auto fbsPoints = fbsSearcher->points();
        _points.resize(fbsPoints->size());
        for (uint32_t i = 0; i < fbsPoints->size(); ++i) {
            _points[i] = fbsToJet(*fbsPoints->Get(i));
        }

        // Copy key/tables
        auto fbsKeys = fbsSearcher->keys();
        _keys.resize(fbsKeys->size());
        for (uint32_t i = 0; i < fbsKeys->size(); ++i) {
            _keys[i] = static_cast<size_t>(fbsKeys->Get(i));
        }

        auto fbsStartIndexTable = fbsSearcher->startIndexTable();
        _startIndexTable.resize(fbsStartIndexTable->size());
        for (uint32_t i = 0; i < fbsStartIndexTable->size(); ++i) {
            _startIndexTable[i] = static_cast<size_t>(fbsStartIndexTable->Get(i));
        }

        auto fbsEndIndexTable = fbsSearcher->endIndexTable();
        _endIndexTable.resize(fbsEndIndexTable->size());
        for (uint32_t i = 0; i < fbsEndIndexTable->size(); ++i) {
            _endIndexTable[i] = static_cast<size_t>(fbsEndIndexTable->Get(i));
        }

        auto fbsSortedIndices = fbsSearcher->sortedIndices();
        _sortedIndices.resize(fbsSortedIndices->size());
        for (uint32_t i = 0; i < fbsSortedIndices->size(); ++i) {
            _sortedIndices[i] = static_cast<size_t>(fbsSortedIndices->Get(i));
        }
    }

    PointParallelHashGridSearcher3::Builder
        PointParallelHashGridSearcher3::builder() {
        return Builder();
    }


    PointParallelHashGridSearcher3::Builder&
        PointParallelHashGridSearcher3::Builder::withResolution(
            const Size3& resolution) {
        _resolution = resolution;
        return *this;
    }

    PointParallelHashGridSearcher3::Builder&
        PointParallelHashGridSearcher3::Builder::withGridSpacing(
            double gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    PointParallelHashGridSearcher3
        PointParallelHashGridSearcher3::Builder::build() const {
        return PointParallelHashGridSearcher3(_resolution, _gridSpacing);
    }

    PointParallelHashGridSearcher3Ptr
        PointParallelHashGridSearcher3::Builder::makeShared() const {
        return std::shared_ptr<PointParallelHashGridSearcher3>(
            new PointParallelHashGridSearcher3(_resolution, _gridSpacing),
            [](PointParallelHashGridSearcher3* obj) {
                delete obj;
            });
    }

    PointNeighborSearcher3Ptr
        PointParallelHashGridSearcher3::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }

}  // namespace jet

#endif  // INCLUDE_JET_POINT_PARALLEL_HASH_GRID_SEARCHER3_H_