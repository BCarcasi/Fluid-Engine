#ifndef INCLUDE_JET_CELL_CENTERED_SCALAR_GRID2_H_
#define INCLUDE_JET_CELL_CENTERED_SCALAR_GRID2_H_

#include "scalar_grid2.h"
#include <utility>  

#include "pch.h"

#include <algorithm>
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief 2-D Cell-centered scalar grid structure.
    //!
    //! This class represents 2-D cell-centered scalar grid which extends
    //! ScalarGrid2. As its name suggests, the class defines the data point at the
    //! center of a grid cell. Thus, the dimension of data points are equal to the
    //! dimension of the cells.
    //!
    class CellCenteredScalarGrid2 final : public ScalarGrid2 {
    public:
        JET_GRID2_TYPE_NAME(CellCenteredScalarGrid2)

            class Builder;

        //! Constructs zero-sized grid.
        CellCenteredScalarGrid2();

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        CellCenteredScalarGrid2(
            size_t resolutionX,
            size_t resolutionY,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double initialValue = 0.0);

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        CellCenteredScalarGrid2(
            const Size2& resolution,
            const Vector2D& gridSpacing = Vector2D(1.0, 1.0),
            const Vector2D& origin = Vector2D(),
            double initialValue = 0.0);

        //! Copy constructor.
        CellCenteredScalarGrid2(const CellCenteredScalarGrid2& other);

        //! Returns the actual data point size.
        Size2 dataSize() const override;

        //! Returns data position for the grid point at (0, 0).
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        Vector2D dataOrigin() const override;

        //!
        //! \brief Swaps the contents with the given \p other grid.
        //!
        //! This function swaps the contents of the grid instance with the given
        //! grid object \p other only if \p other has the same type with this grid.
        //!
        void swap(Grid2* other) override;

        //! Sets the contents with the given \p other grid.
        void set(const CellCenteredScalarGrid2& other);

        //! Sets the contents with the given \p other grid.
        CellCenteredScalarGrid2& operator=(const CellCenteredScalarGrid2& other);

        //! Returns the copy of the grid instance.
        std::shared_ptr<ScalarGrid2> clone() const override;

        //! Returns builder fox CellCenteredScalarGrid2.
        static Builder builder();
    };

    //! Shared pointer for the CellCenteredScalarGrid2 type.
    typedef std::shared_ptr<CellCenteredScalarGrid2> CellCenteredScalarGrid2Ptr;


    //!
    //! \brief Front-end to create CellCenteredScalarGrid2 objects step by step.
    //!
    class CellCenteredScalarGrid2::Builder final : public ScalarGridBuilder2 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size2& resolution);

        //! Returns builder with resolution.
        Builder& withResolution(size_t resolutionX, size_t resolutionY);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(const Vector2D& gridSpacing);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(double gridSpacingX, double gridSpacingY);

        //! Returns builder with grid origin.
        Builder& withOrigin(const Vector2D& gridOrigin);

        //! Returns builder with grid origin.
        Builder& withOrigin(double gridOriginX, double gridOriginY);

        //! Returns builder with initial value.
        Builder& withInitialValue(double initialVal);

        //! Builds CellCenteredScalarGrid2 instance.
        CellCenteredScalarGrid2 build() const;

        //! Builds shared pointer of CellCenteredScalarGrid2 instance.
        CellCenteredScalarGrid2Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of CellCenteredScalarGrid2 instance.
        //!
        //! This is an overriding function that implements ScalarGridBuilder2.
        //!
        ScalarGrid2Ptr build(
            const Size2& resolution,
            const Vector2D& gridSpacing,
            const Vector2D& gridOrigin,
            double initialVal) const override;

    private:
        Size2 _resolution{ 1, 1 };
        Vector2D _gridSpacing{ 1, 1 };
        Vector2D _gridOrigin{ 0, 0 };
        double _initialVal = 0.0;
    };

    CellCenteredScalarGrid2::CellCenteredScalarGrid2() {
    }

    CellCenteredScalarGrid2::CellCenteredScalarGrid2(
        size_t resolutionX,
        size_t resolutionY,
        double gridSpacingX,
        double gridSpacingY,
        double originX,
        double originY,
        double initialValue) {
        resize(
            resolutionX,
            resolutionY,
            gridSpacingX,
            gridSpacingY,
            originX,
            originY,
            initialValue);
    }

    CellCenteredScalarGrid2::CellCenteredScalarGrid2(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& origin,
        double initialValue) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    CellCenteredScalarGrid2::CellCenteredScalarGrid2(
        const CellCenteredScalarGrid2& other) {
        set(other);
    }

    Size2 CellCenteredScalarGrid2::dataSize() const {
        // The size of the data should be the same as the grid resolution.
        return resolution();
    }

    Vector2D CellCenteredScalarGrid2::dataOrigin() const {
        return origin() + 0.5 * gridSpacing();
    }

    std::shared_ptr<ScalarGrid2> CellCenteredScalarGrid2::clone() const {
        return CLONE_W_CUSTOM_DELETER(CellCenteredScalarGrid2);
    }

    void CellCenteredScalarGrid2::swap(Grid2* other) {
        CellCenteredScalarGrid2* sameType
            = dynamic_cast<CellCenteredScalarGrid2*>(other);
        if (sameType != nullptr) {
            swapScalarGrid(sameType);
        }
    }

    void CellCenteredScalarGrid2::set(const CellCenteredScalarGrid2& other) {
        setScalarGrid(other);
    }

    CellCenteredScalarGrid2&
        CellCenteredScalarGrid2::operator=(const CellCenteredScalarGrid2& other) {
        set(other);
        return *this;
    }

    CellCenteredScalarGrid2::Builder CellCenteredScalarGrid2::builder() {
        return Builder();
    }


    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withResolution(const Size2& resolution) {
        _resolution = resolution;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withResolution(
            size_t resolutionX, size_t resolutionY) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withGridSpacing(const Vector2D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withGridSpacing(
            double gridSpacingX, double gridSpacingY) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withOrigin(const Vector2D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withOrigin(
            double gridOriginX, double gridOriginY) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        return *this;
    }

    CellCenteredScalarGrid2::Builder&
        CellCenteredScalarGrid2::Builder::withInitialValue(double initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    CellCenteredScalarGrid2 CellCenteredScalarGrid2::Builder::build() const {
        return CellCenteredScalarGrid2(
            _resolution,
            _gridSpacing,
            _gridOrigin,
            _initialVal);
    }

    ScalarGrid2Ptr CellCenteredScalarGrid2::Builder::build(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& gridOrigin,
        double initialVal) const {
        return std::shared_ptr<CellCenteredScalarGrid2>(
            new CellCenteredScalarGrid2(
                resolution,
                gridSpacing,
                gridOrigin,
                initialVal),
            [](CellCenteredScalarGrid2* obj) {
                delete obj;
            });
    }

    CellCenteredScalarGrid2Ptr
        CellCenteredScalarGrid2::Builder::makeShared() const {
        return std::shared_ptr<CellCenteredScalarGrid2>(
            new CellCenteredScalarGrid2(
                _resolution,
                _gridSpacing,
                _gridOrigin,
                _initialVal),
            [](CellCenteredScalarGrid2* obj) {
                delete obj;
            });
    }
}  // namespace jet

#endif  // INCLUDE_JET_CELL_CENTERED_SCALAR_GRID2_H_