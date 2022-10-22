#ifndef INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER3_H_
#define INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER3_H_

#include "point_neighbor_searcher3.h"
#include <vector>

#include "pch.h"
#include "fbs_helpers.h"
#include "point_simple_list_searcher3_generated.h"


#include <algorithm>
#include <vector>
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief Simple ad-hoc 3-D point searcher.
    //!
    //! This class implements 3-D point searcher simply by looking up every point in
    //! the list. Thus, this class is not ideal for searches involing large number
    //! of points, but only for small set of items.
    //!
    class PointSimpleListSearcher3 final : public PointNeighborSearcher3 {
    public:
        JET_NEIGHBOR_SEARCHER3_TYPE_NAME(PointSimpleListSearcher3)

            class Builder;

        //! Default constructor.
        PointSimpleListSearcher3();

        //! Copy constructor.
        PointSimpleListSearcher3(const PointSimpleListSearcher3& other);

        //!
        //! \brief      Builds internal structure for given points list.
        //!
        //! For this class, this function simply copies the given point list to the
        //! internal list.
        //!
        //! \param[in]  points The points to search.
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
        //! \brief      Creates a new instance of the object with same properties
        //!             than original.
        //!
        //! \return     Copy of this object.
        //!
        PointNeighborSearcher3Ptr clone() const override;

        //! Assignment operator.
        PointSimpleListSearcher3& operator=(const PointSimpleListSearcher3& other);

        //! Copy from the other instance.
        void set(const PointSimpleListSearcher3& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointSimpleListSearcher3.
        static Builder builder();

    private:
        std::vector<Vector3D> _points;
    };

    //! Shared pointer for the PointSimpleListSearcher3 type.
    typedef std::shared_ptr<PointSimpleListSearcher3> PointSimpleListSearcher3Ptr;

    //!
    //! \brief Front-end to create PointSimpleListSearcher3 objects step by step.
    //!
    class PointSimpleListSearcher3::Builder final
        : public PointNeighborSearcherBuilder3 {
    public:
        //! Builds PointSimpleListSearcher3 instance.
        PointSimpleListSearcher3 build() const;

        //! Builds shared pointer of PointSimpleListSearcher3 instance.
        PointSimpleListSearcher3Ptr makeShared() const;

        //! Returns shared pointer of PointNeighborSearcher3 type.
        PointNeighborSearcher3Ptr buildPointNeighborSearcher() const override;
    };

    PointSimpleListSearcher3::PointSimpleListSearcher3() {
    }

    PointSimpleListSearcher3::PointSimpleListSearcher3(
        const PointSimpleListSearcher3& other) {
        set(other);
    }

    void PointSimpleListSearcher3::build(
        const ConstArrayAccessor1<Vector3D>& points) {
        _points.resize(points.size());
        std::copy(points.data(), points.data() + points.size(), _points.begin());
    }

    void PointSimpleListSearcher3::forEachNearbyPoint(
        const Vector3D& origin,
        double radius,
        const std::function<void(size_t, const Vector3D&)>& callback) const {
        double radiusSquared = radius * radius;
        for (size_t i = 0; i < _points.size(); ++i) {
            Vector3D r = _points[i] - origin;
            double distanceSquared = r.dot(r);
            if (distanceSquared <= radiusSquared) {
                callback(i, _points[i]);
            }
        }
    }

    bool PointSimpleListSearcher3::hasNearbyPoint(
        const Vector3D& origin,
        double radius) const {
        double radiusSquared = radius * radius;
        for (size_t i = 0; i < _points.size(); ++i) {
            Vector3D r = _points[i] - origin;
            double distanceSquared = r.dot(r);
            if (distanceSquared <= radiusSquared) {
                return true;
            }
        }

        return false;
    }

    PointNeighborSearcher3Ptr PointSimpleListSearcher3::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointSimpleListSearcher3);
    }

    PointSimpleListSearcher3&
        PointSimpleListSearcher3::operator=(const PointSimpleListSearcher3& other) {
        set(other);
        return *this;
    }

    void PointSimpleListSearcher3::set(const PointSimpleListSearcher3& other) {
        _points = other._points;
    }

    void PointSimpleListSearcher3::serialize(
        std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        // Copy points
        std::vector<fbs::Vector3D> points;
        for (const auto& pt : _points) {
            points.push_back(jetToFbs(pt));
        }

        auto fbsPoints
            = builder.CreateVectorOfStructs(points.data(), points.size());

        // Copy the searcher
        auto fbsSearcher = fbs::CreatePointSimpleListSearcher3(builder, fbsPoints);

        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointSimpleListSearcher3::deserialize(
        const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointSimpleListSearcher3(buffer.data());

        // Copy points
        auto fbsPoints = fbsSearcher->points();
        _points.resize(fbsPoints->size());
        for (uint32_t i = 0; i < fbsPoints->size(); ++i) {
            _points[i] = fbsToJet(*fbsPoints->Get(i));
        }
    }

    PointSimpleListSearcher3
        PointSimpleListSearcher3::Builder::build() const {
        return PointSimpleListSearcher3();
    }

    PointSimpleListSearcher3Ptr
        PointSimpleListSearcher3::Builder::makeShared() const {
        return std::shared_ptr<PointSimpleListSearcher3>(
            new PointSimpleListSearcher3(),
            [](PointSimpleListSearcher3* obj) {
                delete obj;
            });
    }

    PointNeighborSearcher3Ptr
        PointSimpleListSearcher3::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }

}  // namespace jet

#endif  // INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER3_H_