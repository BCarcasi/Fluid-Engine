#ifndef INCLUDE_JET_GRID2_H_
#define INCLUDE_JET_GRID2_H_

#include "bounding_box2.h"
#include "serialization.h"
#include "size2.h"

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

namespace jet {

    //!
    //! \brief Abstract base class for 2-D cartesian grid structure.
    //!
    //! This class represents 2-D cartesian grid structure. This class is an
    //! abstract base class and does not store any data. The class only stores the
    //! shape of the grid. The grid structure is axis-aligned and can have different
    //! grid spacing per axis.
    //!
    class Grid2 : public Serializable {
    public:
        //! Function type for mapping data index to actual position.
        typedef std::function<Vector2D(size_t, size_t)> DataPositionFunc;

        //! Constructs an empty grid.
        Grid2();

        //! Default destructor.
        virtual ~Grid2();

        //! Returns the type name of derived grid.
        virtual std::string typeName() const = 0;

        //! Returns the grid resolution.
        const Size2& resolution() const;

        //! Returns the grid origin.
        const Vector2D& origin() const;

        //! Returns the grid spacing.
        const Vector2D& gridSpacing() const;

        //! Returns the bounding box of the grid.
        const BoundingBox2D& boundingBox() const;

        //! Returns the function that maps grid index to the cell-center position.
        DataPositionFunc cellCenterPosition() const;

        //!
        //! \brief Invokes the given function \p func for each grid cell.
        //!
        //! This function invokes the given function object \p func for each grid
        //! cell in serial manner. The input parameters are i and j indices of a
        //! grid cell. The order of execution is i-first, j-last.
        //!
        void forEachCellIndex(
            const std::function<void(size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each grid cell parallelly.
        //!
        //! This function invokes the given function object \p func for each grid
        //! cell in parallel manner. The input parameters are i and j indices of a
        //! grid cell. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachCellIndex(
            const std::function<void(size_t, size_t)>& func) const;

        //! Returns true if resolution, grid-spacing and origin are same.
        bool hasSameShape(const Grid2& other) const;

        //! Swaps the data with other grid.
        virtual void swap(Grid2* other) = 0;

    protected:
        //! Sets the size parameters including the resolution, grid spacing, and
        //! origin.
        void setSizeParameters(const Size2& resolution, const Vector2D& gridSpacing,
            const Vector2D& origin);

        //! Swaps the size parameters with given grid \p other.
        void swapGrid(Grid2* other);

        //! Sets the size parameters with given grid \p other.
        void setGrid(const Grid2& other);

        //! Fetches the data into a continuous linear array.
        virtual void getData(std::vector<double>* data) const = 0;

        //! Sets the data from a continuous linear array.
        virtual void setData(const std::vector<double>& data) = 0;

    private:
        Size2 _resolution;
        Vector2D _gridSpacing = Vector2D(1, 1);
        Vector2D _origin;
        BoundingBox2D _boundingBox = BoundingBox2D(Vector2D(), Vector2D());
    };

    typedef std::shared_ptr<Grid2> Grid2Ptr;

#define JET_GRID2_TYPE_NAME(DerivedClassName) \
    std::string typeName() const override { return #DerivedClassName; }

    Grid2::Grid2() {}

    Grid2::~Grid2() {}

    const Size2& Grid2::resolution() const { return _resolution; }

    const Vector2D& Grid2::origin() const { return _origin; }

    const Vector2D& Grid2::gridSpacing() const { return _gridSpacing; }

    const BoundingBox2D& Grid2::boundingBox() const { return _boundingBox; }

    Grid2::DataPositionFunc Grid2::cellCenterPosition() const {
        Vector2D h = _gridSpacing;
        Vector2D o = _origin;
        return [h, o](size_t i, size_t j) {
            return o + h * Vector2D(i + 0.5, j + 0.5);
        };
    }

    void Grid2::forEachCellIndex(
        const std::function<void(size_t, size_t)>& func) const {
        serialFor(kZeroSize, _resolution.x, kZeroSize, _resolution.y,
            [&func](size_t i, size_t j) { func(i, j); });
    }

    void Grid2::parallelForEachCellIndex(
        const std::function<void(size_t, size_t)>& func) const {
        parallelFor(kZeroSize, _resolution.x, kZeroSize, _resolution.y,
            [&func](size_t i, size_t j) { func(i, j); });
    }

    bool Grid2::hasSameShape(const Grid2& other) const {
        return _resolution.x == other._resolution.x &&
            _resolution.y == other._resolution.y &&
            similar(_gridSpacing.x, other._gridSpacing.x) &&
            similar(_gridSpacing.y, other._gridSpacing.y) &&
            similar(_origin.x, other._origin.x) &&
            similar(_origin.y, other._origin.y);
    }

    void Grid2::setSizeParameters(const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& origin) {
        _resolution = resolution;
        _origin = origin;
        _gridSpacing = gridSpacing;

        Vector2D resolutionD = Vector2D(static_cast<double>(resolution.x),
            static_cast<double>(resolution.y));

        _boundingBox = BoundingBox2D(origin, origin + gridSpacing * resolutionD);
    }

    void Grid2::swapGrid(Grid2* other) {
        std::swap(_resolution, other->_resolution);
        std::swap(_gridSpacing, other->_gridSpacing);
        std::swap(_origin, other->_origin);
        std::swap(_boundingBox, other->_boundingBox);
    }

    void Grid2::setGrid(const Grid2& other) {
        _resolution = other._resolution;
        _gridSpacing = other._gridSpacing;
        _origin = other._origin;
        _boundingBox = other._boundingBox;
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID2_H_