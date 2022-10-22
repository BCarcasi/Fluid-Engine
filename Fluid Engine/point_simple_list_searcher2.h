#ifndef INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER2_H_
#define INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER2_H_

#include "point_neighbor_searcher2.h"
#include <vector>

#include "pch.h"
#include "fbs_helpers.h"
#include "point_simple_list_searcher2_generated.h"


#include <algorithm>
#include "private_helpers.h"


namespace jet {

    //!
    //! \brief Simple ad-hoc 2-D point searcher.
    //!
    //! This class implements 2-D point searcher simply by looking up every point in
    //! the list. Thus, this class is not ideal for searches involing large number
    //! of points, but only for small set of items.
    //!
    class PointSimpleListSearcher2 final : public PointNeighborSearcher2 {
    public:
        JET_NEIGHBOR_SEARCHER2_TYPE_NAME(PointSimpleListSearcher2)

            class Builder;

        //! Default constructor.
        PointSimpleListSearcher2();

        //! Copy constructor.
        PointSimpleListSearcher2(const PointSimpleListSearcher2& other);

        //!
        //! \brief      Builds internal structure for given points list.
        //!
        //! For this class, this function simply copies the given point list to the
        //! internal list.
        //!
        //! \param[in]  points The points to search.
        //!
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
        //! \brief      Creates a new instance of the object with same properties
        //!             than original.
        //!
        //! \return     Copy of this object.
        //!
        PointNeighborSearcher2Ptr clone() const override;

        //! Assignment operator.
        PointSimpleListSearcher2& operator=(const PointSimpleListSearcher2& other);

        //! Copy from the other instance.
        void set(const PointSimpleListSearcher2& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointSimpleListSearcher2.
        static Builder builder();

    private:
        std::vector<Vector2D> _points;
    };

    //! Shared pointer for the PointSimpleListSearcher2 type.
    typedef std::shared_ptr<PointSimpleListSearcher2> PointSimpleListSearcher2Ptr;

    //!
    //! \brief Front-end to create PointSimpleListSearcher2 objects step by step.
    //!
    class PointSimpleListSearcher2::Builder final
        : public PointNeighborSearcherBuilder2 {
    public:
        //! Builds PointSimpleListSearcher2 instance.
        PointSimpleListSearcher2 build() const;

        //! Builds shared pointer of PointSimpleListSearcher2 instance.
        PointSimpleListSearcher2Ptr makeShared() const;

        //! Returns shared pointer of PointNeighborSearcher3 type.
        PointNeighborSearcher2Ptr buildPointNeighborSearcher() const override;
    };

    PointSimpleListSearcher2::PointSimpleListSearcher2() {
    }

    PointSimpleListSearcher2::PointSimpleListSearcher2(
        const PointSimpleListSearcher2& other) {
        set(other);
    }

    void PointSimpleListSearcher2::build(
        const ConstArrayAccessor1<Vector2D>& points) {
        _points.resize(points.size());
        std::copy(points.data(), points.data() + points.size(), _points.begin());
    }

    void PointSimpleListSearcher2::forEachNearbyPoint(
        const Vector2D& origin,
        double radius,
        const ForEachNearbyPointFunc& callback) const {
        double radiusSquared = radius * radius;
        for (size_t i = 0; i < _points.size(); ++i) {
            Vector2D r = _points[i] - origin;
            double distanceSquared = r.dot(r);
            if (distanceSquared <= radiusSquared) {
                callback(i, _points[i]);
            }
        }
    }

    bool PointSimpleListSearcher2::hasNearbyPoint(
        const Vector2D& origin,
        double radius) const {
        double radiusSquared = radius * radius;
        for (size_t i = 0; i < _points.size(); ++i) {
            Vector2D r = _points[i] - origin;
            double distanceSquared = r.dot(r);
            if (distanceSquared <= radiusSquared) {
                return true;
            }
        }

        return false;
    }

    PointNeighborSearcher2Ptr PointSimpleListSearcher2::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointSimpleListSearcher2);
    }

    PointSimpleListSearcher2&
        PointSimpleListSearcher2::operator=(const PointSimpleListSearcher2& other) {
        set(other);
        return *this;
    }

    void PointSimpleListSearcher2::set(const PointSimpleListSearcher2& other) {
        _points = other._points;
    }

    void PointSimpleListSearcher2::serialize(
        std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        // Copy points
        std::vector<fbs::Vector2D> points;
        for (const auto& pt : _points) {
            points.push_back(jetToFbs(pt));
        }

        auto fbsPoints
            = builder.CreateVectorOfStructs(points.data(), points.size());

        // Copy the searcher
        auto fbsSearcher = fbs::CreatePointSimpleListSearcher2(builder, fbsPoints);

        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointSimpleListSearcher2::deserialize(
        const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointSimpleListSearcher2(buffer.data());

        // Copy points
        auto fbsPoints = fbsSearcher->points();
        _points.resize(fbsPoints->size());
        for (uint32_t i = 0; i < fbsPoints->size(); ++i) {
            _points[i] = fbsToJet(*fbsPoints->Get(i));
        }
    }

    PointSimpleListSearcher2
        PointSimpleListSearcher2::Builder::build() const {
        return PointSimpleListSearcher2();
    }

    PointSimpleListSearcher2Ptr
        PointSimpleListSearcher2::Builder::makeShared() const {
        return std::shared_ptr<PointSimpleListSearcher2>(
            new PointSimpleListSearcher2(),
            [](PointSimpleListSearcher2* obj) {
                delete obj;
            });
    }

    PointNeighborSearcher2Ptr
        PointSimpleListSearcher2::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }

}  // namespace jet

#endif  // INCLUDE_JET_POINT_SIMPLE_LIST_SEARCHER2_H_