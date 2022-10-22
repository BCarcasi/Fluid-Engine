#ifndef INCLUDE_JET_VECTOR_GRID3_H_
#define INCLUDE_JET_VECTOR_GRID3_H_

#include "array_accessor3.h"
#include "grid3.h"
#include "vector_field3.h"
#include <memory>
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "fbs_helpers.h"
#include "vector_grid3_generated.h"

#include "array_samplers3.h"

#include "flatbuffers.h"

#include <algorithm>
#include <string>

namespace jet {

    //! Abstract base class for 3-D vector grid structure.
    class VectorGrid3 : public VectorField3, public Grid3 {
    public:
        //! Read-write array accessor type.
        typedef ArrayAccessor3<Vector3D> VectorDataAccessor;

        //! Read-only array accessor type.
        typedef ConstArrayAccessor3<Vector3D> ConstVectorDataAccessor;

        //! Constructs an empty grid.
        VectorGrid3();

        //! Default destructor.
        virtual ~VectorGrid3();

        //! Clears the contents of the grid.
        void clear();

        //! Resizes the grid using given parameters.
        void resize(
            size_t resolutionX,
            size_t resolutionY,
            size_t resolutionZ,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double gridSpacingZ = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double originZ = 0.0,
            double initialValueX = 0.0,
            double initialValueY = 0.0,
            double initialValueZ = 0.0);

        //! Resizes the grid using given parameters.
        void resize(
            const Size3& resolution,
            const Vector3D& gridSpacing = Vector3D(1, 1, 1),
            const Vector3D& origin = Vector3D(),
            const Vector3D& initialValue = Vector3D());

        //! Resizes the grid using given parameters.
        void resize(
            double gridSpacingX,
            double gridSpacingY,
            double gridSpacingZ,
            double originX,
            double originY,
            double originZ);

        //! Resizes the grid using given parameters.
        void resize(const Vector3D& gridSpacing, const Vector3D& origin);

        //! Fills the grid with given value.
        virtual void fill(const Vector3D& value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) = 0;

        //! Fills the grid with given position-to-value mapping function.
        virtual void fill(const std::function<Vector3D(const Vector3D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) = 0;

        //! Returns the copy of the grid instance.
        virtual std::shared_ptr<VectorGrid3> clone() const = 0;

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
        virtual void onResize(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& origin,
            const Vector3D& initialValue) = 0;
    };

    //! Shared pointer for the VectorGrid3 type.
    typedef std::shared_ptr<VectorGrid3> VectorGrid3Ptr;

    //! Abstract base class for 3-D vector grid builder.
    class VectorGridBuilder3 {
    public:
        //! Creates a builder.
        VectorGridBuilder3();

        //! Default destructor.
        virtual ~VectorGridBuilder3();

        //! Returns 3-D vector grid with given parameters.
        virtual VectorGrid3Ptr build(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& gridOrigin,
            const Vector3D& initialVal) const = 0;
    };

    //! Shared pointer for the VectorGridBuilder3 type.
    typedef std::shared_ptr<VectorGridBuilder3> VectorGridBuilder3Ptr;

    VectorGrid3::VectorGrid3() {}

    VectorGrid3::~VectorGrid3() {}

    void VectorGrid3::clear() {
        resize(Size3(), gridSpacing(), origin(), Vector3D());
    }

    void VectorGrid3::resize(size_t resolutionX, size_t resolutionY,
        size_t resolutionZ, double gridSpacingX,
        double gridSpacingY, double gridSpacingZ,
        double originX, double originY, double originZ,
        double initialValueX, double initialValueY,
        double initialValueZ) {
        resize(Size3(resolutionX, resolutionY, resolutionZ),
            Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
            Vector3D(originX, originY, originZ),
            Vector3D(initialValueX, initialValueY, initialValueZ));
    }

    void VectorGrid3::resize(const Size3& resolution, const Vector3D& gridSpacing,
        const Vector3D& origin, const Vector3D& initialValue) {
        setSizeParameters(resolution, gridSpacing, origin);

        onResize(resolution, gridSpacing, origin, initialValue);
    }

    void VectorGrid3::resize(double gridSpacingX, double gridSpacingY,
        double gridSpacingZ, double originX, double originY,
        double originZ) {
        resize(Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
            Vector3D(originX, originY, originZ));
    }

    void VectorGrid3::resize(const Vector3D& gridSpacing, const Vector3D& origin) {
        resize(resolution(), gridSpacing, origin);
    }

    void VectorGrid3::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        auto fbsResolution = jetToFbs(resolution());
        auto fbsGridSpacing = jetToFbs(gridSpacing());
        auto fbsOrigin = jetToFbs(origin());

        std::vector<double> gridData;
        getData(&gridData);
        auto data = builder.CreateVector(gridData.data(), gridData.size());

        auto fbsGrid = fbs::CreateVectorGrid3(builder, &fbsResolution,
            &fbsGridSpacing, &fbsOrigin, data);

        builder.Finish(fbsGrid);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void VectorGrid3::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsGrid = fbs::GetVectorGrid3(buffer.data());

        resize(fbsToJet(*fbsGrid->resolution()), fbsToJet(*fbsGrid->gridSpacing()),
            fbsToJet(*fbsGrid->origin()));

        auto data = fbsGrid->data();
        std::vector<double> gridData(data->size());
        std::copy(data->begin(), data->end(), gridData.begin());

        setData(gridData);
    }

    VectorGridBuilder3::VectorGridBuilder3() {}

    VectorGridBuilder3::~VectorGridBuilder3() {}

}  // namespace jet

#endif  // INCLUDE_JET_VECTOR_GRID3_H_