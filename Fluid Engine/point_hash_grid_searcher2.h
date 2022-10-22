#ifndef INCLUDE_JET_POINT_HASH_GRID_SEARCHER2_H_
#define INCLUDE_JET_POINT_HASH_GRID_SEARCHER2_H_

#include "point_neighbor_searcher2.h"
#include "point2.h"
#include "size2.h"

#include <vector>

#include "pch.h"

#include "fbs_helpers.h"
#include "point_hash_grid_searcher2_generated.h"

#include "array1.h"

#include <algorithm>
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief Hash grid-based 2-D point searcher.
    //!
    //! This class implements 2-D point searcher by using hash grid for its internal
    //! acceleration data structure. Each point is recorded to its corresponding
    //! bucket where the hashing function is 2-D grid mapping.
    //!
    class PointHashGridSearcher2 final : public PointNeighborSearcher2 {
    public:
        JET_NEIGHBOR_SEARCHER2_TYPE_NAME(PointHashGridSearcher2)

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
        PointHashGridSearcher2(const Size2& resolution, double gridSpacing);

        //!
        //! \brief      Constructs hash grid with given resolution and grid spacing.
        //!
        //! This constructor takes hash grid resolution and its grid spacing as
        //! its input parameters. The grid spacing must be 2x or greater than
        //! search radius.
        //!
        //! \param[in]  resolutionX The resolution x.
        //! \param[in]  resolutionY The resolution y.
        //! \param[in]  gridSpacing The grid spacing.
        //!
        PointHashGridSearcher2(
            size_t resolutionX,
            size_t resolutionY,
            double gridSpacing);

        //! Copy constructor.
        PointHashGridSearcher2(const PointHashGridSearcher2& other);

        //! Builds internal acceleration structure for given points list.
        void build(const ConstArrayAccessor1<Vector2D>& points) override;

        //!
        //! Invokes the callback function for each nearby point around the origin
        //! within given radius.
        //!
        //! \param[in]  origin   The origin position.
        //! \param[in]  radius   The search radius.
        //! \param[in]  callback The callback function.
        //!
        void forEachNearbyPoint(
            const Vector2D& origin,
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
            const Vector2D& origin, double radius) const override;

        //!
        //! \brief      Adds a single point to the hash grid.
        //!
        //! This function adds a single point to the hash grid for future queries.
        //! It can be used for a hash grid that is already built by calling function
        //! PointHashGridSearcher2::build.
        //!
        //! \param[in]  point The point to be added.
        //!
        void add(const Vector2D& point);

        //!
        //! \brief      Returns the internal bucket.
        //!
        //! A bucket is a list of point indices that has same hash value. This
        //! function returns the (immutable) internal bucket structure.
        //!
        //! \return     List of buckets.
        //!
        const std::vector<std::vector<size_t>>& buckets() const;

        //!
        //! Returns the hash value for given 2-D bucket index.
        //!
        //! \param[in]  bucketIndex The bucket index.
        //!
        //! \return     The hash key from bucket index.
        //!
        size_t getHashKeyFromBucketIndex(const Point2I& bucketIndex) const;

        //!
        //! Gets the bucket index from a point.
        //!
        //! \param[in]  position The position of the point.
        //!
        //! \return     The bucket index.
        //!
        Point2I getBucketIndex(const Vector2D& position) const;

        //!
        //! \brief      Creates a new instance of the object with same properties
        //!             than original.
        //!
        //! \return     Copy of this object.
        //!
        PointNeighborSearcher2Ptr clone() const override;

        //! Assignment operator.
        PointHashGridSearcher2& operator=(const PointHashGridSearcher2& other);

        //! Copy from the other instance.
        void set(const PointHashGridSearcher2& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointHashGridSearcher2.
        static Builder builder();

    private:
        double _gridSpacing = 1.0;
        Point2I _resolution = Point2I(1, 1);
        std::vector<Vector2D> _points;
        std::vector<std::vector<size_t>> _buckets;

        size_t getHashKeyFromPosition(const Vector2D& position) const;

        void getNearbyKeys(const Vector2D& position, size_t* bucketIndices) const;
    };

    //! Shared pointer for the PointHashGridSearcher2 type.
    typedef std::shared_ptr<PointHashGridSearcher2> PointHashGridSearcher2Ptr;

    //!
    //! \brief Front-end to create PointHashGridSearcher2 objects step by step.
    //!
    class PointHashGridSearcher2::Builder final
        : public PointNeighborSearcherBuilder2 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size2& resolution);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(double gridSpacing);

        //! Builds PointHashGridSearcher2 instance.
        PointHashGridSearcher2 build() const;

        //! Builds shared pointer of PointHashGridSearcher2 instance.
        PointHashGridSearcher2Ptr makeShared() const;

        //! Returns shared pointer of PointNeighborSearcher3 type.
        PointNeighborSearcher2Ptr buildPointNeighborSearcher() const override;

    private:
        Size2 _resolution{ 64, 64 };
        double _gridSpacing = 1.0;
    };

    PointHashGridSearcher2::PointHashGridSearcher2(
        const Size2& resolution,
        double gridSpacing) :
        PointHashGridSearcher2(resolution.x, resolution.y, gridSpacing) {
    }

    PointHashGridSearcher2::PointHashGridSearcher2(
        size_t resolutionX,
        size_t resolutionY,
        double gridSpacing) :
        _gridSpacing(gridSpacing) {
        _resolution.x = std::max(static_cast<ssize_t>(resolutionX), kOneSSize);
        _resolution.y = std::max(static_cast<ssize_t>(resolutionY), kOneSSize);
    }

    PointHashGridSearcher2::PointHashGridSearcher2(
        const PointHashGridSearcher2& other) {
        set(other);
    }

    void PointHashGridSearcher2::build(
        const ConstArrayAccessor1<Vector2D>& points) {
        _buckets.clear();
        _points.clear();

        // Allocate memory chuncks
        _buckets.resize(_resolution.x * _resolution.y);
        _points.resize(points.size());

        if (points.size() == 0) {
            return;
        }

        // Put points into buckets
        for (size_t i = 0; i < points.size(); ++i) {
            _points[i] = points[i];
            size_t key = getHashKeyFromPosition(points[i]);
            _buckets[key].push_back(i);
        }
    }

    void PointHashGridSearcher2::forEachNearbyPoint(
        const Vector2D& origin,
        double radius,
        const ForEachNearbyPointFunc& callback) const {
        if (_buckets.empty()) {
            return;
        }

        size_t nearbyKeys[4];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 4; i++) {
            const auto& bucket = _buckets[nearbyKeys[i]];
            size_t numberOfPointsInBucket = bucket.size();

            for (size_t j = 0; j < numberOfPointsInBucket; ++j) {
                size_t pointIndex = bucket[j];
                double rSquared = (_points[pointIndex] - origin).lengthSquared();
                if (rSquared <= queryRadiusSquared) {
                    callback(pointIndex, _points[pointIndex]);
                }
            }
        }
    }

    bool PointHashGridSearcher2::hasNearbyPoint(
        const Vector2D& origin,
        double radius) const {
        if (_buckets.empty()) {
            return false;
        }

        size_t nearbyKeys[4];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 4; i++) {
            const auto& bucket = _buckets[nearbyKeys[i]];
            size_t numberOfPointsInBucket = bucket.size();

            for (size_t j = 0; j < numberOfPointsInBucket; ++j) {
                size_t pointIndex = bucket[j];
                double rSquared = (_points[pointIndex] - origin).lengthSquared();
                if (rSquared <= queryRadiusSquared) {
                    return true;
                }
            }
        }

        return false;
    }

    void PointHashGridSearcher2::add(const Vector2D& point) {
        if (_buckets.empty()) {
            Array1<Vector2D> arr = { point };
            build(arr);
        }
        else {
            size_t i = _points.size();
            _points.push_back(point);
            size_t key = getHashKeyFromPosition(point);
            _buckets[key].push_back(i);
        }
    }

    const std::vector<std::vector<size_t>>&
        PointHashGridSearcher2::buckets() const {
        return _buckets;
    }

    Point2I PointHashGridSearcher2::getBucketIndex(const Vector2D& position) const {
        Point2I bucketIndex;
        bucketIndex.x = static_cast<ssize_t>(
            std::floor(position.x / _gridSpacing));
        bucketIndex.y = static_cast<ssize_t>(
            std::floor(position.y / _gridSpacing));
        return bucketIndex;
    }

    size_t PointHashGridSearcher2::getHashKeyFromPosition(
        const Vector2D& position) const {
        Point2I bucketIndex = getBucketIndex(position);

        return getHashKeyFromBucketIndex(bucketIndex);
    }

    size_t PointHashGridSearcher2::getHashKeyFromBucketIndex(
        const Point2I& bucketIndex) const {
        Point2I wrappedIndex = bucketIndex;
        wrappedIndex.x = bucketIndex.x % _resolution.x;
        wrappedIndex.y = bucketIndex.y % _resolution.y;
        if (wrappedIndex.x < 0) {
            wrappedIndex.x += _resolution.x;
        }
        if (wrappedIndex.y < 0) {
            wrappedIndex.y += _resolution.y;
        }
        return static_cast<size_t>(wrappedIndex.y * _resolution.x + wrappedIndex.x);
    }

    void PointHashGridSearcher2::getNearbyKeys(
        const Vector2D& position,
        size_t* nearbyKeys) const {
        Point2I originIndex = getBucketIndex(position), nearbyBucketIndices[4];

        for (int i = 0; i < 4; i++) {
            nearbyBucketIndices[i] = originIndex;
        }

        if ((originIndex.x + 0.5f) * _gridSpacing <= position.x) {
            nearbyBucketIndices[2].x += 1;
            nearbyBucketIndices[3].x += 1;
        }
        else {
            nearbyBucketIndices[2].x -= 1;
            nearbyBucketIndices[3].x -= 1;
        }

        if ((originIndex.y + 0.5f) * _gridSpacing <= position.y) {
            nearbyBucketIndices[1].y += 1;
            nearbyBucketIndices[3].y += 1;
        }
        else {
            nearbyBucketIndices[1].y -= 1;
            nearbyBucketIndices[3].y -= 1;
        }

        for (int i = 0; i < 4; i++) {
            nearbyKeys[i] = getHashKeyFromBucketIndex(nearbyBucketIndices[i]);
        }
    }

    PointNeighborSearcher2Ptr PointHashGridSearcher2::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointHashGridSearcher2);
    }

    PointHashGridSearcher2&
        PointHashGridSearcher2::operator=(const PointHashGridSearcher2& other) {
        set(other);
        return *this;
    }

    void PointHashGridSearcher2::set(const PointHashGridSearcher2& other) {
        _gridSpacing = other._gridSpacing;
        _resolution = other._resolution;
        _points = other._points;
        _buckets = other._buckets;
    }

    void PointHashGridSearcher2::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        // Copy simple data
        auto fbsResolution = fbs::Size2(_resolution.x, _resolution.y);

        // Copy points
        std::vector<fbs::Vector2D> points;
        for (const auto& pt : _points) {
            points.push_back(jetToFbs(pt));
        }

        auto fbsPoints
            = builder.CreateVectorOfStructs(points.data(), points.size());

        // Copy buckets
        std::vector<flatbuffers::Offset<fbs::PointHashGridSearcherBucket2>> buckets;
        for (const auto& bucket : _buckets) {
            std::vector<uint64_t> bucket64(bucket.begin(), bucket.end());
            flatbuffers::Offset<fbs::PointHashGridSearcherBucket2> fbsBucket
                = fbs::CreatePointHashGridSearcherBucket2(
                    builder,
                    builder.CreateVector(bucket64.data(), bucket64.size()));
            buckets.push_back(fbsBucket);
        }

        auto fbsBuckets = builder.CreateVector(buckets);

        // Copy the searcher
        auto fbsSearcher = fbs::CreatePointHashGridSearcher2(
            builder, _gridSpacing, &fbsResolution, fbsPoints, fbsBuckets);

        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointHashGridSearcher2::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointHashGridSearcher2(buffer.data());

        // Copy simple data
        auto res = fbsToJet(*fbsSearcher->resolution());
        _resolution.set({ res.x, res.y });
        _gridSpacing = fbsSearcher->gridSpacing();

        // Copy points
        auto fbsPoints = fbsSearcher->points();
        _points.resize(fbsPoints->size());
        for (uint32_t i = 0; i < fbsPoints->size(); ++i) {
            _points[i] = fbsToJet(*fbsPoints->Get(i));
        }

        // Copy buckets
        auto fbsBuckets = fbsSearcher->buckets();
        _buckets.resize(fbsBuckets->size());
        for (uint32_t i = 0; i < fbsBuckets->size(); ++i) {
            auto fbsBucket = fbsBuckets->Get(i);
            _buckets[i].resize(fbsBucket->data()->size());
            std::transform(
                fbsBucket->data()->begin(),
                fbsBucket->data()->end(),
                _buckets[i].begin(),
                [](uint64_t val) {
                    return static_cast<size_t>(val);
                });
        }
    }

    PointHashGridSearcher2::Builder PointHashGridSearcher2::builder() {
        return Builder();
    }


    PointHashGridSearcher2::Builder&
        PointHashGridSearcher2::Builder::withResolution(const Size2& resolution) {
        _resolution = resolution;
        return *this;
    }

    PointHashGridSearcher2::Builder&
        PointHashGridSearcher2::Builder::withGridSpacing(double gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    PointHashGridSearcher2
        PointHashGridSearcher2::Builder::build() const {
        return PointHashGridSearcher2(_resolution, _gridSpacing);
    }

    PointHashGridSearcher2Ptr
        PointHashGridSearcher2::Builder::makeShared() const {
        return std::shared_ptr<PointHashGridSearcher2>(
            new PointHashGridSearcher2(_resolution, _gridSpacing),
            [](PointHashGridSearcher2* obj) {
                delete obj;
            });
    }

    PointNeighborSearcher2Ptr
        PointHashGridSearcher2::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }
}  // namespace jet

#endif  // INCLUDE_JET_POINT_HASH_GRID_SEARCHER2_H_