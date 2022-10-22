#ifndef INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID2_H_
#define INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID2_H_

#include "collocated_vector_grid2.h"
#include <utility>  // just make cpplint happy..

#include "array_samplers2.h"
#include "parallel.h"
#include "pch.h"

namespace jet {

    //!
    //! \brief 2-D Vertex-centered vector grid structure.
    //!
    //! This class represents 2-D vertex-centered vector grid which extends
    //! CollocatedVectorGrid2. As its name suggests, the class defines the data
    //! point at the grid vertices (corners). Thus, A x B grid resolution will have
    //! (A+1) x (B+1) data points.
    //!
    class VertexCenteredVectorGrid2 final : public CollocatedVectorGrid2 {
    public:
        JET_GRID2_TYPE_NAME(VertexCenteredVectorGrid2)

            class Builder;

        //! Constructs zero-sized grid.
        VertexCenteredVectorGrid2();

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredVectorGrid2(
            size_t resolutionX,
            size_t resolutionY,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double initialValueU = 0.0,
            double initialValueV = 0.0);

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredVectorGrid2(
            const Size2& resolution,
            const Vector2D& gridSpacing = Vector2D(1.0, 1.0),
            const Vector2D& origin = Vector2D(),
            const Vector2D& initialValue = Vector2D());

        //! Copy constructor.
        VertexCenteredVectorGrid2(const VertexCenteredVectorGrid2& other);

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
        void set(const VertexCenteredVectorGrid2& other);

        //! Sets the contents with the given \p other grid.
        VertexCenteredVectorGrid2& operator=(
            const VertexCenteredVectorGrid2& other);

        //! Fills the grid with given value.
        void fill(const Vector2D& value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Fills the grid with given function.
        void fill(const std::function<Vector2D(const Vector2D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Returns the copy of the grid instance.
        std::shared_ptr<VectorGrid2> clone() const override;

        //! Returns builder fox VertexCenteredVectorGrid2.
        static Builder builder();
    };

    //! Shared pointer for the VertexCenteredVectorGrid2 type.
    typedef std::shared_ptr<VertexCenteredVectorGrid2> VertexCenteredVectorGrid2Ptr;


    //! A grid builder class that returns 2-D vertex-centered vector grid.
    class VertexCenteredVectorGrid2::Builder final : public VectorGridBuilder2 {
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
        Builder& withInitialValue(const Vector2D& initialVal);

        //! Returns builder with initial value.
        Builder& withInitialValue(double initialValX, double initialValY);

        //! Builds VertexCenteredVectorGrid2 instance.
        VertexCenteredVectorGrid2 build() const;

        //! Builds shared pointer of VertexCenteredVectorGrid2 instance.
        VertexCenteredVectorGrid2Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of VertexCenteredVectorGrid2 instance.
        //!
        //! This is an overriding function that implements VectorGridBuilder2.
        //!
        VectorGrid2Ptr build(
            const Size2& resolution,
            const Vector2D& gridSpacing,
            const Vector2D& gridOrigin,
            const Vector2D& initialVal) const override;

    private:
        Size2 _resolution{ 1, 1 };
        Vector2D _gridSpacing{ 1, 1 };
        Vector2D _gridOrigin{ 0, 0 };
        Vector2D _initialVal{ 0, 0 };
    };

    VertexCenteredVectorGrid2::VertexCenteredVectorGrid2() {}

    VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(
        size_t resolutionX, size_t resolutionY, double gridSpacingX,
        double gridSpacingY, double originX, double originY, double initialValueU,
        double initialValueV) {
        resize(resolutionX, resolutionY, gridSpacingX, gridSpacingY, originX,
            originY, initialValueU, initialValueV);
    }

    VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(
        const Size2& resolution, const Vector2D& gridSpacing,
        const Vector2D& origin, const Vector2D& initialValue) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    VertexCenteredVectorGrid2::VertexCenteredVectorGrid2(
        const VertexCenteredVectorGrid2& other) {
        set(other);
    }

    Size2 VertexCenteredVectorGrid2::dataSize() const {
        if (resolution() != Size2(0, 0)) {
            return resolution() + Size2(1, 1);
        }
        else {
            return Size2(0, 0);
        }
    }

    Vector2D VertexCenteredVectorGrid2::dataOrigin() const { return origin(); }

    void VertexCenteredVectorGrid2::swap(Grid2* other) {
        VertexCenteredVectorGrid2* sameType =
            dynamic_cast<VertexCenteredVectorGrid2*>(other);
        if (sameType != nullptr) {
            swapCollocatedVectorGrid(sameType);
        }
    }

    void VertexCenteredVectorGrid2::set(const VertexCenteredVectorGrid2& other) {
        setCollocatedVectorGrid(other);
    }

    VertexCenteredVectorGrid2& VertexCenteredVectorGrid2::operator=(
        const VertexCenteredVectorGrid2& other) {
        set(other);
        return *this;
    }

    void VertexCenteredVectorGrid2::fill(const Vector2D& value,
        ExecutionPolicy policy) {
        Size2 size = dataSize();
        auto acc = dataAccessor();
        parallelFor(kZeroSize, size.x, kZeroSize, size.y,
            [value, &acc](size_t i, size_t j) { acc(i, j) = value; },
            policy);
    }

    void VertexCenteredVectorGrid2::fill(
        const std::function<Vector2D(const Vector2D&)>& func,
        ExecutionPolicy policy) {
        Size2 size = dataSize();
        auto acc = dataAccessor();
        DataPositionFunc pos = dataPosition();
        parallelFor(kZeroSize, size.x, kZeroSize, size.y,
            [&func, &acc, &pos](size_t i, size_t j) {
                acc(i, j) = func(pos(i, j));
            },
            policy);
    }

    std::shared_ptr<VectorGrid2> VertexCenteredVectorGrid2::clone() const {
        return CLONE_W_CUSTOM_DELETER(VertexCenteredVectorGrid2);
    }

    VertexCenteredVectorGrid2::Builder VertexCenteredVectorGrid2::builder() {
        return Builder();
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withResolution(const Size2& resolution) {
        _resolution = resolution;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withResolution(size_t resolutionX,
            size_t resolutionY) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withGridSpacing(
            const Vector2D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withGridSpacing(double gridSpacingX,
            double gridSpacingY) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withOrigin(const Vector2D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withOrigin(double gridOriginX,
            double gridOriginY) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withInitialValue(
            const Vector2D& initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    VertexCenteredVectorGrid2::Builder&
        VertexCenteredVectorGrid2::Builder::withInitialValue(double initialValX,
            double initialValY) {
        _initialVal.x = initialValX;
        _initialVal.y = initialValY;
        return *this;
    }

    VertexCenteredVectorGrid2 VertexCenteredVectorGrid2::Builder::build() const {
        return VertexCenteredVectorGrid2(_resolution, _gridSpacing, _gridOrigin,
            _initialVal);
    }

    VertexCenteredVectorGrid2Ptr VertexCenteredVectorGrid2::Builder::makeShared()
        const {
        return std::shared_ptr<VertexCenteredVectorGrid2>(
            new VertexCenteredVectorGrid2(_resolution, _gridSpacing, _gridOrigin,
                _initialVal),
            [](VertexCenteredVectorGrid2* obj) { delete obj; });
    }

    VectorGrid2Ptr VertexCenteredVectorGrid2::Builder::build(
        const Size2& resolution, const Vector2D& gridSpacing,
        const Vector2D& gridOrigin, const Vector2D& initialVal) const {
        return std::shared_ptr<VertexCenteredVectorGrid2>(
            new VertexCenteredVectorGrid2(resolution, gridSpacing, gridOrigin,
                initialVal),
            [](VertexCenteredVectorGrid2* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID2_H_