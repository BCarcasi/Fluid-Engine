#ifndef INCLUDE_JET_VECTOR_GRID2_H_
#define INCLUDE_JET_VECTOR_GRID2_H_

#include "array_accessor2.h"
#include "grid2.h"
#include "parallel.h"
#include "vector_field2.h"

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "fbs_helpers.h"
#include "vector_grid2_generated.h"

#include "array_samplers2.h"

#include "flatbuffers.h"

#include <algorithm>
#include <string>
#include <vector>

#include <memory>
#include <vector>

namespace jet {

    //! Abstract base class for 2-D vector grid structure.
    class VectorGrid2 : public VectorField2, public Grid2 {
    public:
        //! Read-write array accessor type.
        typedef ArrayAccessor2<Vector2D> VectorDataAccessor;

        //! Read-only array accessor type.
        typedef ConstArrayAccessor2<Vector2D> ConstVectorDataAccessor;

        //! Constructs an empty grid.
        VectorGrid2();

        //! Default destructor.
        virtual ~VectorGrid2();

        //! Clears the contents of the grid.
        void clear();

        //! Resizes the grid using given parameters.
        void resize(size_t resolutionX, size_t resolutionY,
            double gridSpacingX = 1.0, double gridSpacingY = 1.0,
            double originX = 0.0, double originY = 0.0,
            double initialValueX = 0.0, double initialValueY = 0.0);

        //! Resizes the grid using given parameters.
        void resize(const Size2& resolution,
            const Vector2D& gridSpacing = Vector2D(1, 1),
            const Vector2D& origin = Vector2D(),
            const Vector2D& initialValue = Vector2D());

        //! Resizes the grid using given parameters.
        void resize(double gridSpacingX, double gridSpacingY, double originX,
            double originY);

        //! Resizes the grid using given parameters.
        void resize(const Vector2D& gridSpacing, const Vector2D& origin);

        //! Fills the grid with given value.
        virtual void fill(const Vector2D& value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) = 0;

        //! Fills the grid with given position-to-value mapping function.
        virtual void fill(const std::function<Vector2D(const Vector2D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) = 0;

        //! Returns the copy of the grid instance.
        virtual std::shared_ptr<VectorGrid2> clone() const = 0;

        //! Serializes the grid instance to the output buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the input buffer to the grid instance.
        void deserialize(const std::vector<uint8_t>& buffer) override;

    protected:
        //!
        //! \brief Invoked when the resizing happens.
        //!
        //! This callback function is called when the grid gets resized. The
        //! overriding class should allocate the internal storage based on its
        //! data layout scheme.
        //!
        virtual void onResize(const Size2& resolution, const Vector2D& gridSpacing,
            const Vector2D& origin,
            const Vector2D& initialValue) = 0;
    };

    //! Shared pointer for the VectorGrid2 type.
    typedef std::shared_ptr<VectorGrid2> VectorGrid2Ptr;

    //! Abstract base class for 2-D vector grid builder.
    class VectorGridBuilder2 {
    public:
        //! Creates a builder.
        VectorGridBuilder2();

        //! Default destructor.
        virtual ~VectorGridBuilder2();

        //! Returns 2-D vector grid with given parameters.
        virtual VectorGrid2Ptr build(const Size2& resolution,
            const Vector2D& gridSpacing,
            const Vector2D& gridOrigin,
            const Vector2D& initialVal) const = 0;
    };

    //! Shared pointer for the VectorGridBuilder2 type.
    typedef std::shared_ptr<VectorGridBuilder2> VectorGridBuilder2Ptr;

    VectorGrid2::VectorGrid2() {
    }

    VectorGrid2::~VectorGrid2() {
    }

    void VectorGrid2::clear() {
        resize(Size2(), gridSpacing(), origin(), Vector2D());
    }

    void VectorGrid2::resize(
        size_t resolutionX,
        size_t resolutionY,
        double gridSpacingX,
        double gridSpacingY,
        double originX,
        double originY,
        double initialValueX,
        double initialValueY) {
        resize(
            Size2(resolutionX, resolutionY),
            Vector2D(gridSpacingX, gridSpacingY),
            Vector2D(originX, originY),
            Vector2D(initialValueX, initialValueY));
    }

    void VectorGrid2::resize(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& origin,
        const Vector2D& initialValue) {
        setSizeParameters(resolution, gridSpacing, origin);

        onResize(resolution, gridSpacing, origin, initialValue);
    }

    void VectorGrid2::resize(
        double gridSpacingX,
        double gridSpacingY,
        double originX,
        double originY) {
        resize(
            Vector2D(gridSpacingX, gridSpacingY),
            Vector2D(originX, originY));
    }

    void VectorGrid2::resize(const Vector2D& gridSpacing, const Vector2D& origin) {
        resize(resolution(), gridSpacing, origin);
    }

    void VectorGrid2::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        auto fbsResolution = jetToFbs(resolution());
        auto fbsGridSpacing = jetToFbs(gridSpacing());
        auto fbsOrigin = jetToFbs(origin());

        std::vector<double> gridData;
        getData(&gridData);
        auto data = builder.CreateVector(gridData.data(), gridData.size());

        auto fbsGrid = fbs::CreateVectorGrid2(
            builder, &fbsResolution, &fbsGridSpacing, &fbsOrigin, data);

        builder.Finish(fbsGrid);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void VectorGrid2::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsGrid = fbs::GetVectorGrid2(buffer.data());

        resize(
            fbsToJet(*fbsGrid->resolution()),
            fbsToJet(*fbsGrid->gridSpacing()),
            fbsToJet(*fbsGrid->origin()));

        auto data = fbsGrid->data();
        std::vector<double> gridData(data->size());
        std::copy(data->begin(), data->end(), gridData.begin());

        setData(gridData);
    }


    VectorGridBuilder2::VectorGridBuilder2() {
    }

    VectorGridBuilder2::~VectorGridBuilder2() {
    }

}  // namespace jet

#endif  // INCLUDE_JET_VECTOR_GRID2_H_