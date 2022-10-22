#ifndef INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID3_H_
#define INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID3_H_

#include "collocated_vector_grid3.h"
#include <utility>  // just make cpplint happy..

#include "pch.h"

#include "array_samplers3.h"
#include "parallel.h"



namespace jet {

    //!
    //! \brief 3-D Vertex-centered vector grid structure.
    //!
    //! This class represents 3-D vertex-centered vector grid which extends
    //! CollocatedVectorGrid3. As its name suggests, the class defines the data
    //! point at the grid vertices (corners). Thus, A x B x C grid resolution will
    //! have (A+1) x (B+1) x (C+1) data points.
    //!
    class VertexCenteredVectorGrid3 final : public CollocatedVectorGrid3 {
    public:
        JET_GRID3_TYPE_NAME(VertexCenteredVectorGrid3)

            class Builder;

        //! Constructs zero-sized grid.
        VertexCenteredVectorGrid3();

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredVectorGrid3(
            size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double gridSpacingZ = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double originZ = 0.0,
            double initialValueU = 0.0,
            double initialValueV = 0.0,
            double initialValueW = 0.0);

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredVectorGrid3(
            const Size3& resolution,
            const Vector3D& gridSpacing = Vector3D(1.0, 1.0, 1.0),
            const Vector3D& origin = Vector3D(),
            const Vector3D& initialValue = Vector3D());

        //! Returns the actual data point size.
        Size3 dataSize() const override;

        //! Returns data position for the grid point at (0, 0, 0).
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        Vector3D dataOrigin() const override;

        //!
        //! \brief Swaps the contents with the given \p other grid.
        //!
        //! This function swaps the contents of the grid instance with the given
        //! grid object \p other only if \p other has the same type with this grid.
        //!
        void swap(Grid3* other) override;

        //! Fills the grid with given value.
        void fill(const Vector3D& value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Fills the grid with given function.
        void fill(const std::function<Vector3D(const Vector3D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Returns the copy of the grid instance.
        std::shared_ptr<VectorGrid3> clone() const override;

        //! Sets the contents with the given \p other grid.
        void set(const VertexCenteredVectorGrid3& other);

        //! Sets the contents with the given \p other grid.
        VertexCenteredVectorGrid3& operator=(
            const VertexCenteredVectorGrid3& other);

        //! Returns builder fox VertexCenteredVectorGrid3.
        static Builder builder();
    };

    //! Shared pointer for the VertexCenteredVectorGrid3 type.
    typedef std::shared_ptr<VertexCenteredVectorGrid3> VertexCenteredVectorGrid3Ptr;


    //!
    //! \brief Front-end to create VertexCenteredVectorGrid3 objects step by step.
    //!
    class VertexCenteredVectorGrid3::Builder final : public VectorGridBuilder3 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size3& resolution);

        //! Returns builder with resolution.
        Builder& withResolution(
            size_t resolutionX, size_t resolutionY, size_t resolutionZ);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(const Vector3D& gridSpacing);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(
            double gridSpacingX, double gridSpacingY, double gridSpacingZ);

        //! Returns builder with grid origin.
        Builder& withOrigin(const Vector3D& gridOrigin);

        //! Returns builder with grid origin.
        Builder& withOrigin(
            double gridOriginX, double gridOriginY, double gridOriginZ);

        //! Returns builder with initial value.
        Builder& withInitialValue(const Vector3D& initialVal);

        //! Returns builder with initial value.
        Builder& withInitialValue(
            double initialValX, double initialValY, double initialValZ);

        //! Builds VertexCenteredVectorGrid3 instance.
        VertexCenteredVectorGrid3 build() const;

        //! Builds shared pointer of VertexCenteredVectorGrid3 instance.
        VertexCenteredVectorGrid3Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of VertexCenteredVectorGrid3 instance.
        //!
        //! This is an overriding function that implements VectorGridBuilder3.
        //!
        VectorGrid3Ptr build(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin,
            const Vector3D& initialVal) const override;

    private:
        Size3 _resolution{ 1, 1, 1 };
        Vector3D _gridSpacing{ 1, 1, 1 };
        Vector3D _gridOrigin{ 0, 0, 0 };
        Vector3D _initialVal{ 0, 0, 0 };
    };

    VertexCenteredVectorGrid3::VertexCenteredVectorGrid3() {}

    VertexCenteredVectorGrid3::VertexCenteredVectorGrid3(
        size_t resolutionX, size_t resolutionY, size_t resolutionZ,
        double gridSpacingX, double gridSpacingY, double gridSpacingZ,
        double originX, double originY, double originZ, double initialValueU,
        double initialValueV, double initialValueW) {
        resize(resolutionX, resolutionY, resolutionZ, gridSpacingX, gridSpacingY,
            gridSpacingZ, originX, originY, originZ, initialValueU,
            initialValueV, initialValueW);
    }

    VertexCenteredVectorGrid3::VertexCenteredVectorGrid3(
        const Size3& resolution, const Vector3D& gridSpacing,
        const Vector3D& origin, const Vector3D& initialValue) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    Size3 VertexCenteredVectorGrid3::dataSize() const {
        if (resolution() != Size3(0, 0, 0)) {
            return resolution() + Size3(1, 1, 1);
        }
        else {
            return Size3(0, 0, 0);
        }
    }

    Vector3D VertexCenteredVectorGrid3::dataOrigin() const { return origin(); }

    void VertexCenteredVectorGrid3::swap(Grid3* other) {
        VertexCenteredVectorGrid3* sameType =
            dynamic_cast<VertexCenteredVectorGrid3*>(other);
        if (sameType != nullptr) {
            swapCollocatedVectorGrid(sameType);
        }
    }

    void VertexCenteredVectorGrid3::fill(const Vector3D& value,
        ExecutionPolicy policy) {
        Size3 size = dataSize();
        auto acc = dataAccessor();
        parallelFor(
            kZeroSize, size.x, kZeroSize, size.y, kZeroSize, size.z,
            [value, &acc](size_t i, size_t j, size_t k) { acc(i, j, k) = value; },
            policy);
    }

    void VertexCenteredVectorGrid3::fill(
        const std::function<Vector3D(const Vector3D&)>& func,
        ExecutionPolicy policy) {
        Size3 size = dataSize();
        auto acc = dataAccessor();
        DataPositionFunc pos = dataPosition();
        parallelFor(kZeroSize, size.x, kZeroSize, size.y, kZeroSize, size.z,
            [&func, &acc, &pos](size_t i, size_t j, size_t k) {
                acc(i, j, k) = func(pos(i, j, k));
            },
            policy);
    }

    std::shared_ptr<VectorGrid3> VertexCenteredVectorGrid3::clone() const {
        return CLONE_W_CUSTOM_DELETER(VertexCenteredVectorGrid3);
    }

    void VertexCenteredVectorGrid3::set(const VertexCenteredVectorGrid3& other) {
        setCollocatedVectorGrid(other);
    }

    VertexCenteredVectorGrid3& VertexCenteredVectorGrid3::operator=(
        const VertexCenteredVectorGrid3& other) {
        set(other);
        return *this;
    }

    VertexCenteredVectorGrid3::Builder VertexCenteredVectorGrid3::builder() {
        return Builder();
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withResolution(const Size3& resolution) {
        _resolution = resolution;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withResolution(size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        _resolution.z = resolutionZ;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withGridSpacing(
            const Vector3D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withGridSpacing(double gridSpacingX,
            double gridSpacingY,
            double gridSpacingZ) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        _gridSpacing.z = gridSpacingZ;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withOrigin(const Vector3D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withOrigin(double gridOriginX,
            double gridOriginY,
            double gridOriginZ) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        _gridOrigin.z = gridOriginZ;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withInitialValue(
            const Vector3D& initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    VertexCenteredVectorGrid3::Builder&
        VertexCenteredVectorGrid3::Builder::withInitialValue(double initialValX,
            double initialValY,
            double initialValZ) {
        _initialVal.x = initialValX;
        _initialVal.y = initialValY;
        _initialVal.z = initialValZ;
        return *this;
    }

    VertexCenteredVectorGrid3 VertexCenteredVectorGrid3::Builder::build() const {
        return VertexCenteredVectorGrid3(_resolution, _gridSpacing, _gridOrigin,
            _initialVal);
    }

    VertexCenteredVectorGrid3Ptr VertexCenteredVectorGrid3::Builder::makeShared()
        const {
        return std::shared_ptr<VertexCenteredVectorGrid3>(
            new VertexCenteredVectorGrid3(_resolution, _gridSpacing, _gridOrigin,
                _initialVal),
            [](VertexCenteredVectorGrid3* obj) { delete obj; });
    }

    VectorGrid3Ptr VertexCenteredVectorGrid3::Builder::build(
        const Size3& resolution, const Vector3D& gridSpacing,
        const Vector3D& gridOrigin, const Vector3D& initialVal) const {
        return std::shared_ptr<VertexCenteredVectorGrid3>(
            new VertexCenteredVectorGrid3(resolution, gridSpacing, gridOrigin,
                initialVal),
            [](VertexCenteredVectorGrid3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VERTEX_CENTERED_VECTOR_GRID3_H_