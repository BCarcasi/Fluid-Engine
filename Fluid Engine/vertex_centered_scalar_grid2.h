#ifndef INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID2_H_
#define INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID2_H_

#include "array2.h"
#include "scalar_grid2.h"
#include <utility>  // just make cpplint happy..

#include "pch.h"
#include "private_helpers.h"

namespace jet {

    //!
    //! \brief 2-D Vertex-centered scalar grid structure.
    //!
    //! This class represents 2-D vertex-centered scalar grid which extends
    //! ScalarGrid3. As its name suggests, the class defines the data point at the
    //! grid vertices (corners). Thus, A x B grid resolution will have (A+1) x (B+1)
    //! data points.
    //!
    class VertexCenteredScalarGrid2 final : public ScalarGrid2 {
    public:
        JET_GRID2_TYPE_NAME(VertexCenteredScalarGrid2)

            class Builder;

        //! Constructs zero-sized grid.
        VertexCenteredScalarGrid2();

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredScalarGrid2(
            size_t resolutionX,
            size_t resolutionY,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double initialValue = 0.0);

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredScalarGrid2(
            const Size2& resolution,
            const Vector2D& gridSpacing = Vector2D(1.0, 1.0),
            const Vector2D& origin = Vector2D(),
            double initialValue = 0.0);

        //! Copy constructor.
        VertexCenteredScalarGrid2(const VertexCenteredScalarGrid2& other);

        //! Returns the actual data point size.
        Size2 dataSize() const override;

        //! Returns data position for the grid point at (0, 0).
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        Vector2D dataOrigin() const override;

        //! Returns the copy of the grid instance.
        std::shared_ptr<ScalarGrid2> clone() const override;

        //!
        //! \brief Swaps the contents with the given \p other grid.
        //!
        //! This function swaps the contents of the grid instance with the given
        //! grid object \p other only if \p other has the same type with this grid.
        //!
        void swap(Grid2* other) override;

        //! Sets the contents with the given \p other grid.
        void set(const VertexCenteredScalarGrid2& other);

        //! Sets the contents with the given \p other grid.
        VertexCenteredScalarGrid2& operator=(
            const VertexCenteredScalarGrid2& other);

        //! Returns builder fox VertexCenteredScalarGrid2.
        static Builder builder();
    };

    //! Shared pointer for the VertexCenteredScalarGrid2 type.
    typedef std::shared_ptr<VertexCenteredScalarGrid2> VertexCenteredScalarGrid2Ptr;


    //!
    //! \brief Front-end to create VertexCenteredScalarGrid2 objects step by step.
    //!
    class VertexCenteredScalarGrid2::Builder final : public ScalarGridBuilder2 {
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

        //! Builds VertexCenteredScalarGrid2 instance.
        VertexCenteredScalarGrid2 build() const;

        //! Builds shared pointer of VertexCenteredScalarGrid2 instance.
        VertexCenteredScalarGrid2Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of VertexCenteredScalarGrid2 instance.
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

    VertexCenteredScalarGrid2::VertexCenteredScalarGrid2() {
    }

    VertexCenteredScalarGrid2::VertexCenteredScalarGrid2(
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

    VertexCenteredScalarGrid2::VertexCenteredScalarGrid2(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& origin,
        double initialValue) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    VertexCenteredScalarGrid2::VertexCenteredScalarGrid2(
        const VertexCenteredScalarGrid2& other) {
        set(other);
    }

    Size2 VertexCenteredScalarGrid2::dataSize() const {
        if (resolution() != Size2(0, 0)) {
            return resolution() + Size2(1, 1);
        }
        else {
            return Size2(0, 0);
        }
    }

    Vector2D VertexCenteredScalarGrid2::dataOrigin() const {
        return origin();
    }

    std::shared_ptr<ScalarGrid2> VertexCenteredScalarGrid2::clone() const {
        return CLONE_W_CUSTOM_DELETER(VertexCenteredScalarGrid2);
    }

    void VertexCenteredScalarGrid2::swap(Grid2* other) {
        VertexCenteredScalarGrid2* sameType
            = dynamic_cast<VertexCenteredScalarGrid2*>(other);
        if (sameType != nullptr) {
            swapScalarGrid(sameType);
        }
    }

    void VertexCenteredScalarGrid2::set(const VertexCenteredScalarGrid2& other) {
        setScalarGrid(other);
    }

    VertexCenteredScalarGrid2&
        VertexCenteredScalarGrid2::operator=(const VertexCenteredScalarGrid2& other) {
        set(other);
        return *this;
    }

    VertexCenteredScalarGrid2::Builder VertexCenteredScalarGrid2::builder() {
        return Builder();
    }


    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withResolution(const Size2& resolution) {
        _resolution = resolution;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withResolution(
            size_t resolutionX, size_t resolutionY) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withGridSpacing(
            const Vector2D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withGridSpacing(
            double gridSpacingX, double gridSpacingY) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withOrigin(const Vector2D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withOrigin(
            double gridOriginX, double gridOriginY) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        return *this;
    }

    VertexCenteredScalarGrid2::Builder&
        VertexCenteredScalarGrid2::Builder::withInitialValue(double initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    VertexCenteredScalarGrid2 VertexCenteredScalarGrid2::Builder::build() const {
        return VertexCenteredScalarGrid2(
            _resolution,
            _gridSpacing,
            _gridOrigin,
            _initialVal);
    }

    VertexCenteredScalarGrid2Ptr
        VertexCenteredScalarGrid2::Builder::makeShared() const {
        return std::shared_ptr<VertexCenteredScalarGrid2>(
            new VertexCenteredScalarGrid2(
                _resolution,
                _gridSpacing,
                _gridOrigin,
                _initialVal),
            [](VertexCenteredScalarGrid2* obj) {
                delete obj;
            });
    }

    ScalarGrid2Ptr VertexCenteredScalarGrid2::Builder::build(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& gridOrigin,
        double initialVal) const {
        return std::shared_ptr<VertexCenteredScalarGrid2>(
            new VertexCenteredScalarGrid2(
                resolution,
                gridSpacing,
                gridOrigin,
                initialVal),
            [](VertexCenteredScalarGrid2* obj) {
                delete obj;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID2_H_