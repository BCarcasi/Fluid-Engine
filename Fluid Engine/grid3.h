#ifndef INCLUDE_JET_GRID3_H_
#define INCLUDE_JET_GRID3_H_

#include "bounding_box3.h"
#include "serialization.h"
#include "size3.h"

#include <functional>
#include <memory>
#include <string>
#include <utility>  // just make cpplint happy..
#include <vector>

#include "parallel.h"
#include "serial.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

namespace jet {

    //!
    //! \brief Abstract base class for 3-D cartesian grid structure.
    //!
    //! This class represents 3-D cartesian grid structure. This class is an
    //! abstract base class and does not store any data. The class only stores the
    //! shape of the grid. The grid structure is axis-aligned and can have different
    //! grid spacing per axis.
    //!
    class Grid3 : public Serializable {
    public:
        //! Function type for mapping data index to actual position.
        typedef std::function<Vector3D(size_t, size_t, size_t)> DataPositionFunc;

        //! Constructs an empty grid.
        Grid3();

        //! Default destructor.
        virtual ~Grid3();

        //! Returns the type name of derived grid.
        virtual std::string typeName() const = 0;

        //! Returns the grid resolution.
        const Size3& resolution() const;

        //! Returns the grid origin.
        const Vector3D& origin() const;

        //! Returns the grid spacing.
        const Vector3D& gridSpacing() const;

        //! Returns the bounding box of the grid.
        const BoundingBox3D& boundingBox() const;

        //! Returns the function that maps grid index to the cell-center position.
        DataPositionFunc cellCenterPosition() const;

        //!
        //! \brief Invokes the given function \p func for each grid cell.
        //!
        //! This function invokes the given function object \p func for each grid
        //! cell in serial manner. The input parameters are i, j, and k indices of a
        //! grid cell. The order of execution is i-first, j-next, k-last.
        //!
        void forEachCellIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each grid cell parallelly.
        //!
        //! This function invokes the given function object \p func for each grid
        //! cell in parallel manner. The input parameters are i, j, and k indices of
        //! a grid cell. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachCellIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //! Serializes the grid instance to the output buffer.
        virtual void serialize(std::vector<uint8_t>* buffer) const = 0;

        //! Deserializes the input buffer to the grid instance.
        virtual void deserialize(const std::vector<uint8_t>& buffer) = 0;

        //! Returns true if resolution, grid-spacing and origin are same.
        bool hasSameShape(const Grid3& other) const;

        //! Swaps the data with other grid.
        virtual void swap(Grid3* other) = 0;

    protected:
        //! Sets the size parameters including the resolution, grid spacing, and
        //! origin.
        void setSizeParameters(const Size3& resolution, const Vector3D& gridSpacing,
            const Vector3D& origin);

        //! Swaps the size parameters with given grid \p other.
        void swapGrid(Grid3* other);

        //! Sets the size parameters with given grid \p other.
        void setGrid(const Grid3& other);

        //! Fetches the data into a continuous linear array.
        virtual void getData(std::vector<double>* data) const = 0;

        //! Sets the data from a continuous linear array.
        virtual void setData(const std::vector<double>& data) = 0;

    private:
        Size3 _resolution;
        Vector3D _gridSpacing = Vector3D(1, 1, 1);
        Vector3D _origin;
        BoundingBox3D _boundingBox = BoundingBox3D(Vector3D(), Vector3D());
    };

    typedef std::shared_ptr<Grid3> Grid3Ptr;

#define JET_GRID3_TYPE_NAME(DerivedClassName) \
    std::string typeName() const override { return #DerivedClassName; }

    Grid3::Grid3() {}

    Grid3::~Grid3() {}

    const Size3& Grid3::resolution() const { return _resolution; }

    const Vector3D& Grid3::origin() const { return _origin; }

    const Vector3D& Grid3::gridSpacing() const { return _gridSpacing; }

    const BoundingBox3D& Grid3::boundingBox() const { return _boundingBox; }

    Grid3::DataPositionFunc Grid3::cellCenterPosition() const {
        Vector3D h = _gridSpacing;
        Vector3D o = _origin;
        return [h, o](size_t i, size_t j, size_t k) {
            return o + h * Vector3D(i + 0.5, j + 0.5, k + 0.5);
        };
    }

    void Grid3::forEachCellIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        serialFor(kZeroSize, _resolution.x, kZeroSize, _resolution.y, kZeroSize,
            _resolution.z,
            [&func](size_t i, size_t j, size_t k) { func(i, j, k); });
    }

    void Grid3::parallelForEachCellIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        parallelFor(kZeroSize, _resolution.x, kZeroSize, _resolution.y, kZeroSize,
            _resolution.z,
            [&func](size_t i, size_t j, size_t k) { func(i, j, k); });
    }

    bool Grid3::hasSameShape(const Grid3& other) const {
        return _resolution.x == other._resolution.x &&
            _resolution.y == other._resolution.y &&
            _resolution.z == other._resolution.z &&
            similar(_gridSpacing.x, other._gridSpacing.x) &&
            similar(_gridSpacing.y, other._gridSpacing.y) &&
            similar(_gridSpacing.z, other._gridSpacing.z) &&
            similar(_origin.x, other._origin.x) &&
            similar(_origin.y, other._origin.y) &&
            similar(_origin.z, other._origin.z);
    }

    void Grid3::setSizeParameters(const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin) {
        _resolution = resolution;
        _origin = origin;
        _gridSpacing = gridSpacing;

        Vector3D resolutionD = Vector3D(static_cast<double>(resolution.x),
            static_cast<double>(resolution.y),
            static_cast<double>(resolution.z));

        _boundingBox = BoundingBox3D(origin, origin + gridSpacing * resolutionD);
    }

    void Grid3::swapGrid(Grid3* other) {
        std::swap(_resolution, other->_resolution);
        std::swap(_gridSpacing, other->_gridSpacing);
        std::swap(_origin, other->_origin);
        std::swap(_boundingBox, other->_boundingBox);
    }

    void Grid3::setGrid(const Grid3& other) {
        _resolution = other._resolution;
        _gridSpacing = other._gridSpacing;
        _origin = other._origin;
        _boundingBox = other._boundingBox;
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID3_H_