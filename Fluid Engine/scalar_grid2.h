#ifndef INCLUDE_JET_SCALAR_GRID2_H_
#define INCLUDE_JET_SCALAR_GRID2_H_

#include "array2.h"
#include "array_accessor2.h"
#include "array_samplers2.h"
#include "grid2.h"
#include "scalar_field2.h"
#include <memory>
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "fbs_helpers.h"
#include "scalar_grid2_generated.h"

#include "fdm_utils.h"
#include "parallel.h"
#include "scalar_grid2.h"
#include "serial.h"

#include "flatbuffers.h"

#include <algorithm>
#include <string>
#include <utility>  // just make cpplint happy..
#include <vector>

namespace jet {

    //! Abstract base class for 2-D scalar grid structure.
    class ScalarGrid2 : public ScalarField2, public Grid2 {
    public:
        //! Read-write array accessor type.
        typedef ArrayAccessor2<double> ScalarDataAccessor;

        //! Read-only array accessor type.
        typedef ConstArrayAccessor2<double> ConstScalarDataAccessor;

        //! Constructs an empty grid.
        ScalarGrid2();

        //! Default destructor.
        virtual ~ScalarGrid2();

        //!
        //! \brief Returns the size of the grid data.
        //!
        //! This function returns the size of the grid data which is not necessarily
        //! equal to the grid resolution if the data is not stored at cell-center.
        //!
        virtual Size2 dataSize() const = 0;

        //!
        //! \brief Returns the origin of the grid data.
        //!
        //! This function returns data position for the grid point at (0, 0).
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        //!
        virtual Vector2D dataOrigin() const = 0;

        //! Returns the copy of the grid instance.
        virtual std::shared_ptr<ScalarGrid2> clone() const = 0;

        //! Clears the contents of the grid.
        void clear();

        //! Resizes the grid using given parameters.
        void resize(
            size_t resolutionX,
            size_t resolutionY,
            double gridSpacingX = 1.0,
            double gridSpacingY = 1.0,
            double originX = 0.0,
            double originY = 0.0,
            double initialValue = 0.0);

        //! Resizes the grid using given parameters.
        void resize(
            const Size2& resolution,
            const Vector2D& gridSpacing = Vector2D(1, 1),
            const Vector2D& origin = Vector2D(),
            double initialValue = 0.0);

        //! Resizes the grid using given parameters.
        void resize(
            double gridSpacingX,
            double gridSpacingY,
            double originX,
            double originY);

        //! Resizes the grid using given parameters.
        void resize(const Vector2D& gridSpacing, const Vector2D& origin);

        //! Returns the grid data at given data point.
        const double& operator()(size_t i, size_t j) const;

        //! Returns the grid data at given data point.
        double& operator()(size_t i, size_t j);

        //! Returns the gradient vector at given data point.
        Vector2D gradientAtDataPoint(size_t i, size_t j) const;

        //! Returns the Laplacian at given data point.
        double laplacianAtDataPoint(size_t i, size_t j) const;

        //! Returns the read-write data array accessor.
        ScalarDataAccessor dataAccessor();

        //! Returns the read-only data array accessor.
        ConstScalarDataAccessor constDataAccessor() const;

        //! Returns the function that maps data point to its position.
        DataPositionFunc dataPosition() const;

        //! Fills the grid with given value.
        void fill(double value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel);

        //! Fills the grid with given position-to-value mapping function.
        void fill(const std::function<double(const Vector2D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel);

        //!
        //! \brief Invokes the given function \p func for each data point.
        //!
        //! This function invokes the given function object \p func for each data
        //! point in serial manner. The input parameters are i and j indices of a
        //! data point. The order of execution is i-first, j-last.
        //!
        void forEachDataPointIndex(
            const std::function<void(size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each data point
        //! parallelly.
        //!
        //! This function invokes the given function object \p func for each data
        //! point in parallel manner. The input parameters are i and j indices of a
        //! data point. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachDataPointIndex(
            const std::function<void(size_t, size_t)>& func) const;

        // ScalarField2 implementations

        //!
        //! \brief Returns the sampled value at given position \p x.
        //!
        //! This function returns the data sampled at arbitrary position \p x.
        //! The sampling function is linear.
        //!
        double sample(const Vector2D& x) const override;

        //!
        //! \brief Returns the sampler function.
        //!
        //! This function returns the data sampler function object. The sampling
        //! function is linear.
        //!
        std::function<double(const Vector2D&)> sampler() const override;

        //! Returns the gradient vector at given position \p x.
        Vector2D gradient(const Vector2D& x) const override;

        //! Returns the Laplacian at given position \p x.
        double laplacian(const Vector2D& x) const override;

        //! Serializes the grid instance to the output buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes the input buffer to the grid instance.
        void deserialize(const std::vector<uint8_t>& buffer) override;

    protected:
        //! Swaps the data storage and predefined samplers with given grid.
        void swapScalarGrid(ScalarGrid2* other);

        //! Sets the data storage and predefined samplers with given grid.
        void setScalarGrid(const ScalarGrid2& other);

        //! Fetches the data into a continuous linear array.
        void getData(std::vector<double>* data) const override;

        //! Sets the data from a continuous linear array.
        void setData(const std::vector<double>& data) override;

    private:
        Array2<double> _data;
        LinearArraySampler2<double, double> _linearSampler;
        std::function<double(const Vector2D&)> _sampler;

        void resetSampler();
    };

    //! Shared pointer for the ScalarGrid2 type.
    typedef std::shared_ptr<ScalarGrid2> ScalarGrid2Ptr;

    //! Abstract base class for 2-D scalar grid builder.
    class ScalarGridBuilder2 {
    public:
        //! Creates a builder.
        ScalarGridBuilder2();

        //! Default destructor.
        virtual ~ScalarGridBuilder2();

        //! Returns 2-D scalar grid with given parameters.
        virtual ScalarGrid2Ptr build(
            const Size2& resolution,
            const Vector2D& gridSpacing,
            const Vector2D& gridOrigin,
            double initialVal) const = 0;
    };

    //! Shared pointer for the ScalarGridBuilder2 type.
    typedef std::shared_ptr<ScalarGridBuilder2> ScalarGridBuilder2Ptr;

    ScalarGrid2::ScalarGrid2()
        : _linearSampler(LinearArraySampler2<double, double>(
            _data.constAccessor(), Vector2D(1, 1), Vector2D())) {}

    ScalarGrid2::~ScalarGrid2() {}

    void ScalarGrid2::clear() { resize(Size2(), gridSpacing(), origin(), 0.0); }

    void ScalarGrid2::resize(size_t resolutionX, size_t resolutionY,
        double gridSpacingX, double gridSpacingY,
        double originX, double originY, double initialValue) {
        resize(Size2(resolutionX, resolutionY),
            Vector2D(gridSpacingX, gridSpacingY), Vector2D(originX, originY),
            initialValue);
    }

    void ScalarGrid2::resize(const Size2& resolution, const Vector2D& gridSpacing,
        const Vector2D& origin, double initialValue) {
        setSizeParameters(resolution, gridSpacing, origin);

        _data.resize(dataSize(), initialValue);
        resetSampler();
    }

    void ScalarGrid2::resize(double gridSpacingX, double gridSpacingY,
        double originX, double originY) {
        resize(Vector2D(gridSpacingX, gridSpacingY), Vector2D(originX, originY));
    }

    void ScalarGrid2::resize(const Vector2D& gridSpacing, const Vector2D& origin) {
        resize(resolution(), gridSpacing, origin);
    }

    const double& ScalarGrid2::operator()(size_t i, size_t j) const {
        return _data(i, j);
    }

    double& ScalarGrid2::operator()(size_t i, size_t j) { return _data(i, j); }

    Vector2D ScalarGrid2::gradientAtDataPoint(size_t i, size_t j) const {
        return gradient2(_data.constAccessor(), gridSpacing(), i, j);
    }

    double ScalarGrid2::laplacianAtDataPoint(size_t i, size_t j) const {
        return laplacian2(_data.constAccessor(), gridSpacing(), i, j);
    }

    double ScalarGrid2::sample(const Vector2D& x) const { return _sampler(x); }

    std::function<double(const Vector2D&)> ScalarGrid2::sampler() const {
        return _sampler;
    }

    Vector2D ScalarGrid2::gradient(const Vector2D& x) const {
        std::array<Point2UI, 4> indices;
        std::array<double, 4> weights;
        _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);

        Vector2D result;

        for (int i = 0; i < 4; ++i) {
            result += weights[i] * gradientAtDataPoint(indices[i].x, indices[i].y);
        }

        return result;
    }

    double ScalarGrid2::laplacian(const Vector2D& x) const {
        std::array<Point2UI, 4> indices;
        std::array<double, 4> weights;
        _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);

        double result = 0.0;

        for (int i = 0; i < 4; ++i) {
            result += weights[i] * laplacianAtDataPoint(indices[i].x, indices[i].y);
        }

        return result;
    }

    ScalarGrid2::ScalarDataAccessor ScalarGrid2::dataAccessor() {
        return _data.accessor();
    }

    ScalarGrid2::ConstScalarDataAccessor ScalarGrid2::constDataAccessor() const {
        return _data.constAccessor();
    }

    ScalarGrid2::DataPositionFunc ScalarGrid2::dataPosition() const {
        Vector2D o = dataOrigin();
        return [this, o](size_t i, size_t j) -> Vector2D {
            return o + gridSpacing() * Vector2D({ i, j });
        };
    }

    void ScalarGrid2::fill(double value, ExecutionPolicy policy) {
        parallelFor(kZeroSize, _data.width(), kZeroSize, _data.height(),
            [this, value](size_t i, size_t j) { _data(i, j) = value; },
            policy);
    }

    void ScalarGrid2::fill(const std::function<double(const Vector2D&)>& func,
        ExecutionPolicy policy) {
        DataPositionFunc pos = dataPosition();
        parallelFor(kZeroSize, _data.width(), kZeroSize, _data.height(),
            [this, &func, &pos](size_t i, size_t j) {
                _data(i, j) = func(pos(i, j));
            },
            policy);
    }

    void ScalarGrid2::forEachDataPointIndex(
        const std::function<void(size_t, size_t)>& func) const {
        _data.forEachIndex(func);
    }

    void ScalarGrid2::parallelForEachDataPointIndex(
        const std::function<void(size_t, size_t)>& func) const {
        _data.parallelForEachIndex(func);
    }

    void ScalarGrid2::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        auto fbsResolution = jetToFbs(resolution());
        auto fbsGridSpacing = jetToFbs(gridSpacing());
        auto fbsOrigin = jetToFbs(origin());

        std::vector<double> gridData;
        getData(&gridData);
        auto data = builder.CreateVector(gridData.data(), gridData.size());

        auto fbsGrid = fbs::CreateScalarGrid2(builder, &fbsResolution,
            &fbsGridSpacing, &fbsOrigin, data);

        builder.Finish(fbsGrid);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void ScalarGrid2::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsGrid = fbs::GetScalarGrid2(buffer.data());

        resize(fbsToJet(*fbsGrid->resolution()), fbsToJet(*fbsGrid->gridSpacing()),
            fbsToJet(*fbsGrid->origin()));

        auto data = fbsGrid->data();
        std::vector<double> gridData(data->size());
        std::copy(data->begin(), data->end(), gridData.begin());

        setData(gridData);
    }

    void ScalarGrid2::swapScalarGrid(ScalarGrid2* other) {
        swapGrid(other);

        _data.swap(other->_data);
        std::swap(_linearSampler, other->_linearSampler);
        std::swap(_sampler, other->_sampler);
    }

    void ScalarGrid2::setScalarGrid(const ScalarGrid2& other) {
        setGrid(other);

        _data.set(other._data);
        resetSampler();
    }

    void ScalarGrid2::resetSampler() {
        _linearSampler = LinearArraySampler2<double, double>(
            _data.constAccessor(), gridSpacing(), dataOrigin());
        _sampler = _linearSampler.functor();
    }

    void ScalarGrid2::getData(std::vector<double>* data) const {
        size_t size = dataSize().x * dataSize().y;
        data->resize(size);
        std::copy(_data.begin(), _data.end(), data->begin());
    }

    void ScalarGrid2::setData(const std::vector<double>& data) {
        JET_ASSERT(dataSize().x * dataSize().y == data.size());

        std::copy(data.begin(), data.end(), _data.begin());
    }

    ScalarGridBuilder2::ScalarGridBuilder2() {}

    ScalarGridBuilder2::~ScalarGridBuilder2() {}

}  // namespace jet

#endif  // INCLUDE_JET_SCALAR_GRID2_H_