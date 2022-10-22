#ifndef INCLUDE_JET_POINT_KDTREE_SEARCHER3_H
#define INCLUDE_JET_POINT_KDTREE_SEARCHER3_H

#include "kdtree.h"
#include "point3.h"
#include "point_neighbor_searcher3.h"
#include "size3.h"

#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "fbs_helpers.h"
#include "point_kdtree_searcher3_generated.h"

#include "bounding_box3.h"

#include <numeric>
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief KdTree-based 3-D point searcher.
    //!
    //! This class implements 3-D point searcher by using KdTree for its internal
    //! acceleration data structure.
    //!
    class PointKdTreeSearcher3 final : public PointNeighborSearcher3 {
    public:
        JET_NEIGHBOR_SEARCHER3_TYPE_NAME(PointKdTreeSearcher3)

            class Builder;

        //! Constructs an empty kD-tree instance.
        PointKdTreeSearcher3();

        //! Copy constructor.
        PointKdTreeSearcher3(const PointKdTreeSearcher3& other);

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
            const Vector3D& origin, double radius,
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
        bool hasNearbyPoint(const Vector3D& origin, double radius) const override;

        //!
        //! \brief      Creates a new instance of the object with same properties
        //!             than original.
        //!
        //! \return     Copy of this object.
        //!
        PointNeighborSearcher3Ptr clone() const override;

        //! Assignment operator.
        PointKdTreeSearcher3& operator=(const PointKdTreeSearcher3& other);

        //! Copy from the other instance.
        void set(const PointKdTreeSearcher3& other);

        //! Serializes the neighbor searcher into the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the neighbor searcher from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Returns builder fox PointKdTreeSearcher3.
        static Builder builder();

    private:
        KdTree<double, 3> _tree;
    };

    //! Shared pointer for the PointKdTreeSearcher3 type.
    typedef std::shared_ptr<PointKdTreeSearcher3> PointKdTreeSearcher3Ptr;

    //!
    //! \brief Front-end to create PointKdTreeSearcher3 objects step by step.
    //!
    class PointKdTreeSearcher3::Builder final
        : public PointNeighborSearcherBuilder3 {
    public:
        //! Builds PointKdTreeSearcher3 instance.
        PointKdTreeSearcher3 build() const;

        //! Builds shared pointer of PointKdTreeSearcher3 instance.
        PointKdTreeSearcher3Ptr makeShared() const;

        //! Returns shared pointer of PointNeighborSearcher3 type.
        PointNeighborSearcher3Ptr buildPointNeighborSearcher() const override;
    };

    PointKdTreeSearcher3::PointKdTreeSearcher3() {}

    PointKdTreeSearcher3::PointKdTreeSearcher3(const PointKdTreeSearcher3& other) {
        set(other);
    }

    void PointKdTreeSearcher3::build(const ConstArrayAccessor1<Vector3D>& points) {
        _tree.build(points);
    }

    void PointKdTreeSearcher3::forEachNearbyPoint(
        const Vector3D& origin, double radius,
        const ForEachNearbyPointFunc& callback) const {
        _tree.forEachNearbyPoint(origin, radius, callback);
    }

    bool PointKdTreeSearcher3::hasNearbyPoint(const Vector3D& origin,
        double radius) const {
        return _tree.hasNearbyPoint(origin, radius);
    }

    PointNeighborSearcher3Ptr PointKdTreeSearcher3::clone() const {
        return CLONE_W_CUSTOM_DELETER(PointKdTreeSearcher3);
    }

    PointKdTreeSearcher3& PointKdTreeSearcher3::operator=(
        const PointKdTreeSearcher3& other) {
        set(other);
        return *this;
    }

    void PointKdTreeSearcher3::set(const PointKdTreeSearcher3& other) {
        _tree = other._tree;
    }

    void PointKdTreeSearcher3::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        // Copy points
        std::vector<fbs::Vector3D> points;
        for (const auto& iter : _tree) {
            points.push_back(jetToFbs(iter));
        }

        auto fbsPoints =
            builder.CreateVectorOfStructs(points.data(), points.size());

        // Copy nodes
        std::vector<fbs::PointKdTreeSearcherNode3> nodes;
        for (auto iter = _tree.beginNode(); iter != _tree.endNode(); ++iter) {
            nodes.emplace_back(iter->flags, iter->child, iter->item);
        }

        auto fbsNodes = builder.CreateVectorOfStructs(nodes);

        // Copy the searcher
        auto fbsSearcher =
            fbs::CreatePointKdTreeSearcher3(builder, fbsPoints, fbsNodes);

        // Finish
        builder.Finish(fbsSearcher);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void PointKdTreeSearcher3::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsSearcher = fbs::GetPointKdTreeSearcher3(buffer.data());

        auto fbsPoints = fbsSearcher->points();
        auto fbsNodes = fbsSearcher->nodes();

        _tree.reserve(fbsPoints->size(), fbsNodes->size());

        // Copy points
        auto pointsIter = _tree.begin();
        for (uint32_t i = 0; i < fbsPoints->size(); ++i) {
            pointsIter[i] = fbsToJet(*fbsPoints->Get(i));
        }

        // Copy nodes
        auto nodesIter = _tree.beginNode();
        for (uint32_t i = 0; i < fbsNodes->size(); ++i) {
            const auto fbsNode = fbsNodes->Get(i);
            nodesIter[i].flags = fbsNode->flags();
            nodesIter[i].child = fbsNode->child();
            nodesIter[i].item = fbsNode->item();
            nodesIter[i].point = pointsIter[fbsNode->item()];
        }
    }

    PointKdTreeSearcher3::Builder PointKdTreeSearcher3::builder() {
        return Builder{};
    }

    //

    PointKdTreeSearcher3 PointKdTreeSearcher3::Builder::build() const {
        return PointKdTreeSearcher3{};
    }

    PointKdTreeSearcher3Ptr PointKdTreeSearcher3::Builder::makeShared() const {
        return std::shared_ptr<PointKdTreeSearcher3>(
            new PointKdTreeSearcher3,
            [](PointKdTreeSearcher3* obj) { delete obj; });
    }

    PointNeighborSearcher3Ptr
        PointKdTreeSearcher3::Builder::buildPointNeighborSearcher() const {
        return makeShared();
    }

}  // namespace jet

#endif  // INCLUDE_JET_POINT_KDTREE_SEARCHER3_H