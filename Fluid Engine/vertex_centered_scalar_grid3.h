#ifndef INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID3_H_
#define INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID3_H_

#include "array3.h"
#include "scalar_grid3.h"
#include <utility>  // just make cpplint happy..

#include "pch.h"
#include "private_helpers.h"


namespace jet {

    //!
    //! \brief 3-D Vertex-centered scalar grid structure.
    //!
    //! This class represents 3-D vertex-centered scalar grid which extends
    //! ScalarGrid3. As its name suggests, the class defines the data point at the
    //! grid vertices (corners). Thus, A x B x C grid resolution will have
    //! (A+1) x (B+1) x (C+1) data points.
    //!
    class VertexCenteredScalarGrid3 final : public ScalarGrid3 {
    public:
        JET_GRID3_TYPE_NAME(VertexCenteredScalarGrid3)

            class Builder;

        //! Constructs zero-sized grid.
        VertexCenteredScalarGrid3();

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredScalarGrid3(
            size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double gridSpacingZ = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double originZ = 0.0,
            double initialValue = 0.0);

        //! Constructs a grid with given resolution, grid spacing, origin and
        //! initial value.
        VertexCenteredScalarGrid3(
            const Size3& resolution,
            const Vector3D& gridSpacing = Vector3D(1.0, 1.0, 1.0),
            const Vector3D& origin = Vector3D(),
            double initialValue = 0.0);

        //! Copy constructor.
        VertexCenteredScalarGrid3(const VertexCenteredScalarGrid3& other);

        //! Returns the actual data point size.
        Size3 dataSize() const override;

        //! Returns data position for the grid point at (0, 0, 0).
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        Vector3D dataOrigin() const override;

        //! Returns the copy of the grid instance.
        std::shared_ptr<ScalarGrid3> clone() const override;

        //!
        //! \brief Swaps the contents with the given \p other grid.
        //!
        //! This function swaps the contents of the grid instance with the given
        //! grid object \p other only if \p other has the same type with this grid.
        //!
        void swap(Grid3* other) override;

        //! Sets the contents with the given \p other grid.
        void set(const VertexCenteredScalarGrid3& other);

        //! Sets the contents with the given \p other grid.
        VertexCenteredScalarGrid3& operator=(
            const VertexCenteredScalarGrid3& other);

        //! Returns builder fox VertexCenteredScalarGrid3.
        static Builder builder();
    };

    //! Shared pointer for the VertexCenteredScalarGrid3 type.
    typedef std::shared_ptr<VertexCenteredScalarGrid3> VertexCenteredScalarGrid3Ptr;


    //! A grid builder class that returns 3-D vertex-centered scalar grid.
    class VertexCenteredScalarGrid3::Builder final : public ScalarGridBuilder3 {
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
        Builder& withInitialValue(double initialVal);

        //! Builds VertexCenteredScalarGrid3 instance.
        VertexCenteredScalarGrid3 build() const;

        //! Builds shared pointer of VertexCenteredScalarGrid3 instance.
        VertexCenteredScalarGrid3Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of VertexCenteredScalarGrid3 instance.
        //!
        //! This is an overriding function that implements ScalarGridBuilder3.
        //!
        ScalarGrid3Ptr build(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin,
            double initialVal) const override;

    private:
        Size3 _resolution{ 1, 1, 1 };
        Vector3D _gridSpacing{ 1, 1, 1 };
        Vector3D _gridOrigin{ 0, 0, 0 };
        double _initialVal = 0.0;
    };

    VertexCenteredScalarGrid3::VertexCenteredScalarGrid3() {
    }

    VertexCenteredScalarGrid3::VertexCenteredScalarGrid3(
        size_t resolutionX,
        size_t resolutionY,
        size_t resolutionZ,
        double gridSpacingX,
        double gridSpacingY,
        double gridSpacingZ,
        double originX,
        double originY,
        double originZ,
        double initialValue) {
        resize(
            resolutionX,
            resolutionY,
            resolutionZ,
            gridSpacingX,
            gridSpacingY,
            gridSpacingZ,
            originX,
            originY,
            originZ,
            initialValue);
    }

    VertexCenteredScalarGrid3::VertexCenteredScalarGrid3(
        const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin,
        double initialValue) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    VertexCenteredScalarGrid3::VertexCenteredScalarGrid3(
        const VertexCenteredScalarGrid3& other) {
        set(other);
    }

    Size3 VertexCenteredScalarGrid3::dataSize() const {
        if (resolution() != Size3(0, 0, 0)) {
            return resolution() + Size3(1, 1, 1);
        }
        else {
            return Size3(0, 0, 0);
        }
    }

    Vector3D VertexCenteredScalarGrid3::dataOrigin() const {
        return origin();
    }

    std::shared_ptr<ScalarGrid3> VertexCenteredScalarGrid3::clone() const {
        return CLONE_W_CUSTOM_DELETER(VertexCenteredScalarGrid3);
    }

    void VertexCenteredScalarGrid3::swap(Grid3* other) {
        VertexCenteredScalarGrid3* sameType
            = dynamic_cast<VertexCenteredScalarGrid3*>(other);
        if (sameType != nullptr) {
            swapScalarGrid(sameType);
        }
    }

    void VertexCenteredScalarGrid3::set(const VertexCenteredScalarGrid3& other) {
        setScalarGrid(other);
    }

    VertexCenteredScalarGrid3&
        VertexCenteredScalarGrid3::operator=(const VertexCenteredScalarGrid3& other) {
        set(other);
        return *this;
    }

    VertexCenteredScalarGrid3::Builder VertexCenteredScalarGrid3::builder() {
        return Builder();
    }


    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withResolution(const Size3& resolution) {
        _resolution = resolution;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withResolution(
            size_t resolutionX, size_t resolutionY, size_t resolutionZ) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        _resolution.z = resolutionZ;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withGridSpacing(
            const Vector3D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withGridSpacing(
            double gridSpacingX, double gridSpacingY, double gridSpacingZ) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        _gridSpacing.z = gridSpacingZ;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withOrigin(const Vector3D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withOrigin(
            double gridOriginX, double gridOriginY, double gridOriginZ) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        _gridOrigin.z = gridOriginZ;
        return *this;
    }

    VertexCenteredScalarGrid3::Builder&
        VertexCenteredScalarGrid3::Builder::withInitialValue(double initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    VertexCenteredScalarGrid3 VertexCenteredScalarGrid3::Builder::build() const {
        return VertexCenteredScalarGrid3(
            _resolution,
            _gridSpacing,
            _gridOrigin,
            _initialVal);
    }

    VertexCenteredScalarGrid3Ptr
        VertexCenteredScalarGrid3::Builder::makeShared() const {
        return std::shared_ptr<VertexCenteredScalarGrid3>(
            new VertexCenteredScalarGrid3(
                _resolution,
                _gridSpacing,
                _gridOrigin,
                _initialVal),
            [](VertexCenteredScalarGrid3* obj) {
                delete obj;
            });
    }

    ScalarGrid3Ptr VertexCenteredScalarGrid3::Builder::build(
        const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& gridOrigin,
        double initialVal) const {
        return std::shared_ptr<VertexCenteredScalarGrid3>(
            new VertexCenteredScalarGrid3(
                resolution,
                gridSpacing,
                gridOrigin,
                initialVal),
            [](VertexCenteredScalarGrid3* obj) {
                delete obj;
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VERTEX_CENTERED_SCALAR_GRID3_H_