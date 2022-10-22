#ifndef INCLUDE_JET_KDTREE_H
#define INCLUDE_JET_KDTREE_H

#include "bounding_box2.h"
#include "bounding_box3.h"
#include "vector2.h"
#include "vector3.h"
#include "vector4.h"

#include <vector>


#include <numeric>

namespace jet {

    //! Generic k-d tree structure.
    template <typename T, size_t K>
    class KdTree final {
    public:
        typedef Vector<T, K> Point;
        typedef BoundingBox<T, K> BBox;

        //! Simple K-d tree node.
        struct Node {
            //! Split axis if flags < K, leaf indicator if flags == K.
            size_t flags = 0;

            //! \brief Right child index.
            //! Note that left child index is this node index + 1.
            size_t child = kMaxSize;

            //! Item index.
            size_t item = kMaxSize;

            //! Point stored in the node.
            Point point;

            //! Default contructor.
            Node();

            //! Initializes leaf node.
            void initLeaf(size_t it, const Point& pt);

            //! Initializes internal node.
            void initInternal(size_t axis, size_t it, size_t c, const Point& pt);

            //! Returns true if leaf.
            bool isLeaf() const;
        };

        typedef std::vector<Point> ContainerType;
        typedef typename ContainerType::iterator Iterator;
        typedef typename ContainerType::const_iterator ConstIterator;

        typedef std::vector<Node> NodeContainerType;
        typedef typename NodeContainerType::iterator NodeIterator;
        typedef typename NodeContainerType::const_iterator ConstNodeIterator;

        //! Constructs an empty kD-tree instance.
        KdTree();

        //! Builds internal acceleration structure for given points list.
        void build(const ConstArrayAccessor1<Point>& points);

        //!
        //! Invokes the callback function for each nearby point around the origin
        //! within given radius.
        //!
        //! \param[in]  origin   The origin position.
        //! \param[in]  radius   The search radius.
        //! \param[in]  callback The callback function.
        //!
        void forEachNearbyPoint(
            const Point& origin, T radius,
            const std::function<void(size_t, const Point&)>& callback) const;

        //!
        //! Returns true if there are any nearby points for given origin within
        //! radius.
        //!
        //! \param[in]  origin The origin.
        //! \param[in]  radius The radius.
        //!
        //! \return     True if has nearby point, false otherwise.
        //!
        bool hasNearbyPoint(const Point& origin, T radius) const;

        //! Returns index of the nearest point.
        size_t nearestPoint(const Point& origin) const;

        //! Returns the mutable begin iterator of the item.
        Iterator begin();

        //! Returns the mutable end iterator of the item.
        Iterator end();

        //! Returns the immutable begin iterator of the item.
        ConstIterator begin() const;

        //! Returns the immutable end iterator of the item.
        ConstIterator end() const;

        //! Returns the mutable begin iterator of the node.
        NodeIterator beginNode();

        //! Returns the mutable end iterator of the node.
        NodeIterator endNode();

        //! Returns the immutable begin iterator of the node.
        ConstNodeIterator beginNode() const;

        //! Returns the immutable end iterator of the node.
        ConstNodeIterator endNode() const;

        //! Reserves memory space for this tree.
        void reserve(size_t numPoints, size_t numNodes);

    private:
        std::vector<Point> _points;
        std::vector<Node> _nodes;

        size_t build(size_t nodeIndex, size_t* itemIndices, size_t nItems,
            size_t currentDepth);
    };

    template <typename T, size_t K>
    KdTree<T, K>::Node::Node() {
        child = kMaxSize;
    }

    template <typename T, size_t K>
    void KdTree<T, K>::Node::initLeaf(size_t it, const Point& pt) {
        flags = K;
        item = it;
        child = kMaxSize;
        point = pt;
    }

    template <typename T, size_t K>
    void KdTree<T, K>::Node::initInternal(size_t axis, size_t it, size_t c,
        const Point& pt) {
        flags = axis;
        item = it;
        child = c;
        point = pt;
    }

    template <typename T, size_t K>
    bool KdTree<T, K>::Node::isLeaf() const {
        return flags == K;
    }

    //

    template <typename T, size_t K>
    KdTree<T, K>::KdTree() {}

    template <typename T, size_t K>
    void KdTree<T, K>::build(const ConstArrayAccessor1<Point>& points) {
        _points.resize(points.size());
        std::copy(points.begin(), points.end(), _points.begin());

        if (_points.empty()) {
            return;
        }

        _nodes.clear();

        std::vector<size_t> itemIndices(_points.size());
        std::iota(std::begin(itemIndices), std::end(itemIndices), 0);

        build(0, itemIndices.data(), _points.size(), 0);
    }

    template <typename T, size_t K>
    void KdTree<T, K>::forEachNearbyPoint(
        const Point& origin, T radius,
        const std::function<void(size_t, const Point&)>& callback) const {
        const T r2 = radius * radius;

        // prepare to traverse the tree for sphere
        static const int kMaxTreeDepth = 8 * sizeof(size_t);
        const Node* todo[kMaxTreeDepth];
        size_t todoPos = 0;

        // traverse the tree nodes for sphere
        const Node* node = _nodes.data();

        while (node != nullptr) {
            if (node->item != kMaxSize &&
                (node->point - origin).lengthSquared() <= r2) {
                callback(node->item, node->point);
            }

            if (node->isLeaf()) {
                // grab next node to process from todo stack
                if (todoPos > 0) {
                    // Dequeue
                    --todoPos;
                    node = todo[todoPos];
                }
                else {
                    break;
                }
            }
            else {
                // get node children pointers for sphere
                const Node* firstChild = node + 1;
                const Node* secondChild = (Node*)&_nodes[node->child];

                // advance to next child node, possibly enqueue other child
                const size_t axis = node->flags;
                const T plane = node->point[axis];
                if (plane - origin[axis] > radius) {
                    node = firstChild;
                }
                else if (origin[axis] - plane > radius) {
                    node = secondChild;
                }
                else {
                    // enqueue secondChild in todo stack
                    todo[todoPos] = secondChild;
                    ++todoPos;
                    node = firstChild;
                }
            }
        }
    }

    template <typename T, size_t K>
    bool KdTree<T, K>::hasNearbyPoint(const Point& origin, T radius) const {
        const T r2 = radius * radius;

        // prepare to traverse the tree for sphere
        static const int kMaxTreeDepth = 8 * sizeof(size_t);
        const Node* todo[kMaxTreeDepth];
        size_t todoPos = 0;

        // traverse the tree nodes for sphere
        const Node* node = _nodes.data();

        while (node != nullptr) {
            if (node->item != kMaxSize &&
                (node->point - origin).lengthSquared() <= r2) {
                return true;
            }

            if (node->isLeaf()) {
                // grab next node to process from todo stack
                if (todoPos > 0) {
                    // Dequeue
                    --todoPos;
                    node = todo[todoPos];
                }
                else {
                    break;
                }
            }
            else {
                // get node children pointers for sphere
                const Node* firstChild = node + 1;
                const Node* secondChild = (Node*)&_nodes[node->child];

                // advance to next child node, possibly enqueue other child
                const size_t axis = node->flags;
                const T plane = node->point[axis];
                if (origin[axis] < plane && plane - origin[axis] > radius) {
                    node = firstChild;
                }
                else if (origin[axis] > plane && origin[axis] - plane > radius) {
                    node = secondChild;
                }
                else {
                    // enqueue secondChild in todo stack
                    todo[todoPos] = secondChild;
                    ++todoPos;
                    node = firstChild;
                }
            }
        }

        return false;
    }

    template <typename T, size_t K>
    size_t KdTree<T, K>::nearestPoint(const Point& origin) const {
        // prepare to traverse the tree for sphere
        static const int kMaxTreeDepth = 8 * sizeof(size_t);
        const Node* todo[kMaxTreeDepth];
        size_t todoPos = 0;

        // traverse the tree nodes for sphere
        const Node* node = _nodes.data();
        size_t nearest = 0;
        T minDist2 = (node->point - origin).lengthSquared();

        while (node != nullptr) {
            const T newDist2 = (node->point - origin).lengthSquared();
            if (newDist2 <= minDist2) {
                nearest = node->item;
                minDist2 = newDist2;
            }

            if (node->isLeaf()) {
                // grab next node to process from todo stack
                if (todoPos > 0) {
                    // Dequeue
                    --todoPos;
                    node = todo[todoPos];
                }
                else {
                    break;
                }
            }
            else {
                // get node children pointers for sphere
                const Node* firstChild = node + 1;
                const Node* secondChild = (Node*)&_nodes[node->child];

                // advance to next child node, possibly enqueue other child
                const size_t axis = node->flags;
                const T plane = node->point[axis];
                const T minDist = std::sqrt(minDist2);
                if (plane - origin[axis] > minDist) {
                    node = firstChild;
                }
                else if (origin[axis] - plane > minDist) {
                    node = secondChild;
                }
                else {
                    // enqueue secondChild in todo stack
                    todo[todoPos] = secondChild;
                    ++todoPos;
                    node = firstChild;
                }
            }
        }

        return nearest;
    }

    template <typename T, size_t K>
    void KdTree<T, K>::reserve(size_t numPoints, size_t numNodes) {
        _points.resize(numPoints);
        _nodes.resize(numNodes);
    }

    template <typename T, size_t K>
    typename KdTree<T, K>::Iterator KdTree<T, K>::begin() {
        return _points.begin();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::Iterator KdTree<T, K>::end() {
        return _points.end();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::ConstIterator KdTree<T, K>::begin() const {
        return _points.begin();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::ConstIterator KdTree<T, K>::end() const {
        return _points.end();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::NodeIterator KdTree<T, K>::beginNode() {
        return _nodes.begin();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::NodeIterator KdTree<T, K>::endNode() {
        return _nodes.end();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::ConstNodeIterator KdTree<T, K>::beginNode() const {
        return _nodes.begin();
    };

    template <typename T, size_t K>
    typename KdTree<T, K>::ConstNodeIterator KdTree<T, K>::endNode() const {
        return _nodes.end();
    };

    template <typename T, size_t K>
    size_t KdTree<T, K>::build(size_t nodeIndex, size_t* itemIndices, size_t nItems,
        size_t currentDepth) {
        // add a node
        _nodes.emplace_back();

        // initialize leaf node if termination criteria met
        if (nItems == 0) {
            _nodes[nodeIndex].initLeaf(kMaxSize, {});
            return currentDepth + 1;
        }
        if (nItems == 1) {
            _nodes[nodeIndex].initLeaf(itemIndices[0], _points[itemIndices[0]]);
            return currentDepth + 1;
        }

        // choose which axis to split along
        BBox nodeBound;
        for (size_t i = 0; i < nItems; ++i) {
            nodeBound.merge(_points[itemIndices[i]]);
        }
        Point d = nodeBound.upperCorner - nodeBound.lowerCorner;
        size_t axis = static_cast<size_t>(d.dominantAxis());

        // pick mid point
        std::nth_element(itemIndices, itemIndices + nItems / 2,
            itemIndices + nItems, [&](size_t a, size_t b) {
                return _points[a][axis] < _points[b][axis];
            });
        size_t midPoint = nItems / 2;

        // recursively initialize children nodes
        size_t d0 = build(nodeIndex + 1, itemIndices, midPoint, currentDepth + 1);
        _nodes[nodeIndex].initInternal(axis, itemIndices[midPoint], _nodes.size(),
            _points[itemIndices[midPoint]]);
        size_t d1 = build(_nodes[nodeIndex].child, itemIndices + midPoint + 1,
            nItems - midPoint - 1, currentDepth + 1);

        return std::max(d0, d1);
    }

}  // namespace jet

#endif  // INCLUDE_JET_KDTREE_H