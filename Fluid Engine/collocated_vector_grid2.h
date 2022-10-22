#ifndef INCLUDE_JET_COLLOCATED_VECTOR_GRID2_H_
#define INCLUDE_JET_COLLOCATED_VECTOR_GRID2_H_

#include "array2.h"
#include "array_samplers2.h"
#include "vector_grid2.h"
#include <vector>

#include "pch.h"

#include "parallel.h"
#include "serial.h"

#include <algorithm>
#include <utility> 
#include <vector>
#include "private_helpers.h"

namespace jet {

    //! \brief Abstract base class for 2-D collocated vector grid structure.
    class CollocatedVectorGrid2 : public VectorGrid2 {
    public:
        //! Constructs an empty grid.
        CollocatedVectorGrid2();

        //! Default destructor.
        virtual ~CollocatedVectorGrid2();

        //! Returns the actual data point size.
        virtual Size2 dataSize() const = 0;

        //!
        //! \brief Returns data position for the grid point at (0, 0).
        //!
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        //!
        virtual Vector2D dataOrigin() const = 0;

        //! Returns the grid data at given data point.
        const Vector2D& operator()(size_t i, size_t j) const;

        //! Returns the grid data at given data point.
        Vector2D& operator()(size_t i, size_t j);

        //! Returns divergence at data point location.
        double divergenceAtDataPoint(size_t i, size_t j) const;

        //! Returns curl at data point location.
        double curlAtDataPoint(size_t i, size_t j) const;

        //! Returns the read-write data array accessor.
        VectorDataAccessor dataAccessor();

        //! Returns the read-only data array accessor.
        ConstVectorDataAccessor constDataAccessor() const;

        //! Returns the function that maps data point to its position.
        DataPositionFunc dataPosition() const;

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

        // VectorField2 implementations

        //! Returns sampled value at given position \p x.
        Vector2D sample(const Vector2D& x) const override;

        //! Returns divergence at given position \p x.
        double divergence(const Vector2D& x) const override;

        //! Returns curl at given position \p x.
        double curl(const Vector2D& x) const override;

        //!
        //! \brief Returns the sampler function.
        //!
        //! This function returns the data sampler function object. The sampling
        //! function is linear.
        //!
        std::function<Vector2D(const Vector2D&)> sampler() const override;

    protected:
        //! Swaps the data storage and predefined samplers with given grid.
        void swapCollocatedVectorGrid(CollocatedVectorGrid2* other);

        //! Sets the data storage and predefined samplers with given grid.
        void setCollocatedVectorGrid(const CollocatedVectorGrid2& other);

        //! Fetches the data into a continuous linear array.
        void getData(std::vector<double>* data) const override;

        //! Sets the data from a continuous linear array.
        void setData(const std::vector<double>& data) override;

    private:
        Array2<Vector2D> _data;
        LinearArraySampler2<Vector2D, double> _linearSampler;
        std::function<Vector2D(const Vector2D&)> _sampler;

        void onResize(
            const Size2& resolution,
            const Vector2D& gridSpacing,
            const Vector2D& origin,
            const Vector2D& initialValue) final;

        void resetSampler();
    };

    //! Shared pointer for the CollocatedVectorGrid2 type.
    typedef std::shared_ptr<CollocatedVectorGrid2> CollocatedVectorGrid2Ptr;

    CollocatedVectorGrid2::CollocatedVectorGrid2() :
        _linearSampler(_data.constAccessor(), Vector2D(1, 1), Vector2D()) {
    }

    CollocatedVectorGrid2::~CollocatedVectorGrid2() {
    }

    const Vector2D& CollocatedVectorGrid2::operator()(
        size_t i, size_t j) const {
        return _data(i, j);
    }

    Vector2D& CollocatedVectorGrid2::operator()(size_t i, size_t j) {
        return _data(i, j);
    }

    double CollocatedVectorGrid2::divergenceAtDataPoint(
        size_t i, size_t j) const {
        const Size2 ds = _data.size();
        const Vector2D& gs = gridSpacing();

        JET_ASSERT(i < ds.x&& j < ds.y);

        double left = _data((i > 0) ? i - 1 : i, j).x;
        double right = _data((i + 1 < ds.x) ? i + 1 : i, j).x;
        double down = _data(i, (j > 0) ? j - 1 : j).y;
        double up = _data(i, (j + 1 < ds.y) ? j + 1 : j).y;

        return 0.5 * (right - left) / gs.x
            + 0.5 * (up - down) / gs.y;
    }

    double CollocatedVectorGrid2::curlAtDataPoint(
        size_t i, size_t j) const {
        const Size2 ds = _data.size();
        const Vector2D& gs = gridSpacing();

        JET_ASSERT(i < ds.x&& j < ds.y);

        Vector2D left = _data((i > 0) ? i - 1 : i, j);
        Vector2D right = _data((i + 1 < ds.x) ? i + 1 : i, j);
        Vector2D bottom = _data(i, (j > 0) ? j - 1 : j);
        Vector2D top = _data(i, (j + 1 < ds.y) ? j + 1 : j);

        double Fx_ym = bottom.x;
        double Fx_yp = top.x;

        double Fy_xm = left.y;
        double Fy_xp = right.y;

        return 0.5 * (Fy_xp - Fy_xm) / gs.x - 0.5 * (Fx_yp - Fx_ym) / gs.y;
    }

    Vector2D CollocatedVectorGrid2::sample(const Vector2D& x) const {
        return _sampler(x);
    }

    double CollocatedVectorGrid2::divergence(const Vector2D& x) const {
        std::array<Point2UI, 4> indices;
        std::array<double, 4> weights;
        _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);

        double result = 0.0;

        for (int i = 0; i < 4; ++i) {
            result += weights[i] * divergenceAtDataPoint(
                indices[i].x, indices[i].y);
        }

        return result;
    }

    double CollocatedVectorGrid2::curl(const Vector2D& x) const {
        std::array<Point2UI, 4> indices;
        std::array<double, 4> weights;
        _linearSampler.getCoordinatesAndWeights(x, &indices, &weights);

        double result = 0.0;

        for (int i = 0; i < 4; ++i) {
            result += weights[i] * curlAtDataPoint(
                indices[i].x, indices[i].y);
        }

        return result;
    }

    std::function<Vector2D(const Vector2D&)>
        CollocatedVectorGrid2::sampler() const {
        return _sampler;
    }

    VectorGrid2::VectorDataAccessor CollocatedVectorGrid2::dataAccessor() {
        return _data.accessor();
    }

    VectorGrid2::ConstVectorDataAccessor
        CollocatedVectorGrid2::constDataAccessor() const {
        return _data.constAccessor();
    }

    VectorGrid2::DataPositionFunc CollocatedVectorGrid2::dataPosition() const {
        Vector2D dataOrigin_ = dataOrigin();
        return [this, dataOrigin_](size_t i, size_t j) -> Vector2D {
            return dataOrigin_ + gridSpacing() * Vector2D({ i, j });
        };
    }

    void CollocatedVectorGrid2::forEachDataPointIndex(
        const std::function<void(size_t, size_t)>& func) const {
        _data.forEachIndex(func);
    }

    void CollocatedVectorGrid2::parallelForEachDataPointIndex(
        const std::function<void(size_t, size_t)>& func) const {
        _data.parallelForEachIndex(func);
    }

    void CollocatedVectorGrid2::swapCollocatedVectorGrid(
        CollocatedVectorGrid2* other) {
        swapGrid(other);

        _data.swap(other->_data);
        std::swap(_linearSampler, other->_linearSampler);
        std::swap(_sampler, other->_sampler);
    }

    void CollocatedVectorGrid2::setCollocatedVectorGrid(
        const CollocatedVectorGrid2& other) {
        setGrid(other);

        _data.set(other._data);
        resetSampler();
    }

    void CollocatedVectorGrid2::onResize(
        const Size2& resolution,
        const Vector2D& gridSpacing,
        const Vector2D& origin,
        const Vector2D& initialValue) {
        UNUSED_VARIABLE(resolution);
        UNUSED_VARIABLE(gridSpacing);
        UNUSED_VARIABLE(origin);

        _data.resize(dataSize(), initialValue);
        resetSampler();
    }

    void CollocatedVectorGrid2::resetSampler() {
        _linearSampler = LinearArraySampler2<Vector2D, double>(
            _data.constAccessor(), gridSpacing(), dataOrigin());
        _sampler = _linearSampler.functor();
    }

    void CollocatedVectorGrid2::getData(std::vector<double>* data) const {
        size_t size = 2 * dataSize().x * dataSize().y;
        data->resize(size);
        size_t cnt = 0;
        _data.forEach([&](const Vector2D& value) {
            (*data)[cnt++] = value.x;
            (*data)[cnt++] = value.y;
            });
    }

    void CollocatedVectorGrid2::setData(const std::vector<double>& data) {
        JET_ASSERT(2 * dataSize().x * dataSize().y == data.size());

        size_t cnt = 0;
        _data.forEachIndex([&](size_t i, size_t j) {
            _data(i, j).x = data[cnt++];
            _data(i, j).y = data[cnt++];
            });
    }

}  // namespace jet

#endif  // INCLUDE_JET_COLLOCATED_VECTOR_GRID2_H_