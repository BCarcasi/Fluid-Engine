#ifndef INCLUDE_JET_POINT_HASH_GRID_SEARCHER3_H_
#define INCLUDE_JET_POINT_HASH_GRID_SEARCHER3_H_

#include "point_neighbor_searcher3.h"
#include "point3.h"
#include "size3.h"

#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "fbs_helpers.h"
#include "point_hash_grid_searcher3_generated.h"

#include "array1.h"

#include <algorithm>
#include <vector>
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief Hash grid-based 3-D point searcher.
    //!
    //! This class implements 3-D point searcher by using hash grid for its internal
    //! acceleration data structure. Each point is recorded to its corresponding
    //! bucket where the hashing function is 3-D grid mapping.
    //!
    class PointHashGridSearcher3 final : public PointNeighborSearcher3 {
    public:
        JET_NEIGHBOR_SEARCHER3_TYPE_NAME(PointHashGridSearcher3)

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
        PointHashGridSearcher3(const Size3& resolution, double gridSpacing);

        //!
        //! \brief      Constructs hash grid with given resolution and grid spacing.
        //!
        //! This constructor takes hash grid resolution and its grid spacing as
        //! its input parameters. The grid spacing must be 2x or greater than
        //! search radius.
        //!
        //! \param[in]  resolutionX The resolution x.
        //! \param[in]  resolutionY The resolution y.
        //! \param[in]  resolutionY The resolution z.
        //! \param[in]  gridSpacing The grid spacing.
        //!
        PointHashGridSearcher3(
            size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ,
            double gridSpacing);

        //! Copy constructor.
        PointHashGridSearcher3(const PointHashGridSearcher3& other);

        //! Builds internal acceleration structure for given points list.
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
        //! \brief      Adds a single point to the hash grid.
        //!
        //! This function adds a single point to the hash grid for future queries.
        //! It can be used for a hash grid that is already built by calling function
        //! PointHashGridSearcher3::build.
        //!
        //! \param[in]  point The point to be added.
        //!
        void add(const Vector3D& point);

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
        PointHashGridSearcher3& operator=(const PointHashGridSearcher3& other);

        //! Copy from the other instance.
        void set(const PointHashGridSearcher3& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointHashGridSearcher3.
        static Builder builder();

    private:
        double _gridSpacing = 1.0;
        Point3I _resolution = Point3I(1, 1, 1);
        std::vector<Vector3D> _points;
        std::vector<std::vector<size_t>> _buckets;

        size_t getHashKeyFromPosition(const Vector3D& position) const;

        void getNearbyKeys(const Vector3D& position, size_t* bucketIndices) const;
    };

    //! Shared pointer for the PointHashGridSearcher3 type.
    typedef std::shared_ptr<PointHashGridSearcher3> PointHashGridSearcher3Ptr;

    //!
    //! \brief Front-end to create PointHashGridSearcher3 objects step by step.
    //!
    class PointHashGridSearcher3::Builder final
        : public PointNeighborSearcherBuilder3 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size3& resolution);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(double gridSpacing);

        //! Builds PointHashGridSearcher3 instance.
        PointHashGridSearcher3 build() const;

        //! Builds shared pointer of PointHashGridSearcher3 instance.
        PointHashGridSearcher3Ptr makeShared() const;

        //! Returns shared pointer of PointHashGridSearcher3 type.
        PointNeighborSearcher3Ptr buildPointNeighborSearcher() const override;

    private:
        Size3 _resolution{ 64, 64, 64 };
        double _gridSpacing = 1.0;
    };

    PointHashGridSearcher3::PointHashGridSearcher3(
        const Size3& resolution,
        double gridSpacing) :
        PointHashGridSearcher3(
            resolution.x,
            resolution.y,
            resolution.z,
            gridSpacing) {
    }

    PointHashGridSearcher3::PointHashGridSearcher3(
        size_t resolutionX,
        size_t resolutionY,
        size_t resolutionZ,
        double gridSpacing) :
        _gridSpacing(gridSpacing) {
        _resolution.x = std::max(static_cast<ssize_t>(resolutionX), kOneSSize);
        _resolution.y = std::max(static_cast<ssize_t>(resolutionY), kOneSSize);
        _resolution.z = std::max(static_cast<ssize_t>(resolutionZ), kOneSSize);
    }

    PointHashGridSearcher3::PointHashGridSearcher3(
        const PointHashGridSearcher3& other) {
        set(other);
    }

    void PointHashGridSearcher3::build(
        const ConstArrayAccessor1<Vector3D>& points) {
        _buckets.clear();
        _points.clear();

        // Allocate memory chuncks
        _buckets.resize(_resolution.x * _resolution.y * _resolution.z);
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

    void PointHashGridSearcher3::forEachNearbyPoint(
        const Vector3D& origin,
        double radius,
        const std::function<void(size_t, const Vector3D&)>& callback) const {
        if (_buckets.empty()) {
            return;
        }

        size_t nearbyKeys[8];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 8; i++) {
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

    bool PointHashGridSearcher3::hasNearbyPoint(
        const Vector3D& origin,
        double radius) const {
        if (_buckets.empty()) {
            return false;
        }

        size_t nearbyKeys[8];
        getNearbyKeys(origin, nearbyKeys);

        const double queryRadiusSquared = radius * radius;

        for (int i = 0; i < 8; i++) {
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

    void PointHashGridSearcher3::add(const Vector3D& point) {
        if (_buckets.empty()) {
            Array1<Vector3D> arr = { point };
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
        PointHashGridSearcher3::buckets() const {
        return _buckets;
    }

    Point3I PointHashGridSearcher3::getBucketIndex(const Vector3D& position) const {
        Point3I bucketIndex;
        bucketIndex.x = static_cast<ssize_t>(
            std::floor(position.x / _gridSpacing));
        bucketIndex.y = static_cast<ssize_t>(
            std::floor(position.y / _gridSpacing));
        bucketIndex.z = static_cast<ssize_t>(
            std::floor(position.z / _gridSpacing));
        return bucketIndex;
    }

    size_t PointHashGridSearcher3::getHashKeyFromPosition(
        const Vector3D& position) const {
        Point3I bucketIndex = getBucketIndex(position);
        return getHashKeyFromBucketIndex(bucketIndex);
    }

    size_t PointHashGridSearcher3::getHashKeyFromBucketIndex(
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

    void PointHashGridSearcher3::getNearbyKeys(
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

    PointNeighborSearcher3Ptr PointHashGridSearcher3::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointHashGridSearcher3);
    }

    PointHashGridSearcher3&
        PointHashGridSearcher3::operator=(const PointHashGridSearcher3& other) {
        set(other);
        return *this;
    }

    void PointHashGridSearcher3::set(const PointHashGridSearcher3& other) {
        _gridSpacing = other._gridSpacing;
        _resolution = other._resolution;
        _points = other._points;
        _buckets = other._buckets;
    }

    void PointHashGridSearcher3::serialize(std::vector<uint8_t>* buffer) const {
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

        // Copy buckets
        std::vector<flatbuffers::Offset<fbs::PointHashGridSearcherBucket3>> buckets;
        for (const auto& bucket : _buckets) {
            std::vector<uint64_t> bucket64(bucket.begin(), bucket.end());
            flatbuffers::Offset<fbs::PointHashGridSearcherBucket3> fbsBucket
                = fbs::CreatePointHashGridSearcherBucket3(
                    builder,
                    builder.CreateVector(bucket64.data(), bucket64.size()));
            buckets.push_back(fbsBucket);
        }

        auto fbsBuckets = builder.CreateVector(buckets);

        // Copy the searcher
        auto fbsSearcher = fbs::CreatePointHashGridSearcher3(
            builder, _gridSpacing, &fbsResolution, fbsPoints, fbsBuckets);

        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointHashGridSearcher3::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointHashGridSearcher3(buffer.data());

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

    PointHashGridSearcher3::Builder PointHashGridSearcher3::builder() {
        return Builder();
    }


    PointHashGridSearcher3::Builder&
        PointHashGridSearcher3::Builder::withResolution(const Size3& resolution) {
        _resolution = resolution;
        return *this;
    }

    PointHashGridSearcher3::Builder&
        PointHashGridSearcher3::Builder::withGridSpacing(double gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    PointHashGridSearcher3
        PointHashGridSearcher3::Builder::build() const {
        return PointHashGridSearcher3(_resolution, _gridSpacing);
    }

    PointHashGridSearcher3Ptr
        PointHashGridSearcher3::Builder::makeShared() const {
        return std::shared_ptr<PointHashGridSearcher3>(
            new PointHashGridSearcher3(_resolution, _gridSpacing),
            [](PointHashGridSearcher3* obj) {
                delete obj;
            });
    }

    PointNeighborSearcher3Ptr
        PointHashGridSearcher3::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }

}  // namespace jet

#endif  // INCLUDE_JET_POINT_HASH_GRID_SEARCHER3_H_