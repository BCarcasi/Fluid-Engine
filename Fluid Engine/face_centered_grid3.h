#ifndef INCLUDE_JET_FACE_CENTERED_GRID3_H_
#define INCLUDE_JET_FACE_CENTERED_GRID3_H_

#include "array3.h"
#include "array_samplers3.h"
#include "vector_grid3.h"
#include <memory>
#include <utility>  // just make cpplint happy..
#include <vector>

#include "pch.h"


#include "parallel.h"
#include "serial.h"

#include <algorithm>
#include "private_helpers.h"


namespace jet {

    //!
    //! \brief 3-D face-centered (a.k.a MAC or staggered) grid.
    //!
    //! This class implements face-centered grid which is also known as
    //! marker-and-cell (MAC) or staggered grid. This vector grid stores each vector
    //! component at face center. Thus, u, v, and w components are not collocated.
    //!
    class FaceCenteredGrid3 final : public VectorGrid3 {
    public:
        JET_GRID3_TYPE_NAME(FaceCenteredGrid3)

            class Builder;

        //! Read-write scalar data accessor type.
        typedef ArrayAccessor3<double> ScalarDataAccessor;

        //! Read-only scalar data accessor type.
        typedef ConstArrayAccessor3<double> ConstScalarDataAccessor;

        //! Constructs empty grid.
        FaceCenteredGrid3();

        //! Resizes the grid using given parameters.
        FaceCenteredGrid3(size_t resolutionX, size_t resolutionY,
            size_t resolutionZ, double gridSpacingX = 1.0,
            double gridSpacingY = 1.0, double gridSpacingZ = 1.0,
            double originX = 0.0, double originY = 0.0,
            double originZ = 0.0, double initialValueU = 0.0,
            double initialValueV = 0.0, double initialValueW = 0.0);

        //! Resizes the grid using given parameters.
        FaceCenteredGrid3(const Size3& resolution,
            const Vector3D& gridSpacing = Vector3D(1.0, 1.0, 1.0),
            const Vector3D& origin = Vector3D(),
            const Vector3D& initialValue = Vector3D());

        //! Copy constructor.
        FaceCenteredGrid3(const FaceCenteredGrid3& other);

        //!
        //! \brief Swaps the contents with the given \p other grid.
        //!
        //! This function swaps the contents of the grid instance with the given
        //! grid object \p other only if \p other has the same type with this grid.
        //!
        void swap(Grid3* other) override;

        //! Sets the contents with the given \p other grid.
        void set(const FaceCenteredGrid3& other);

        //! Sets the contents with the given \p other grid.
        FaceCenteredGrid3& operator=(const FaceCenteredGrid3& other);

        //! Returns u-value at given data point.
        double& u(size_t i, size_t j, size_t k);

        //! Returns u-value at given data point.
        const double& u(size_t i, size_t j, size_t k) const;

        //! Returns v-value at given data point.
        double& v(size_t i, size_t j, size_t k);

        //! Returns v-value at given data point.
        const double& v(size_t i, size_t j, size_t k) const;

        //! Returns w-value at given data point.
        double& w(size_t i, size_t j, size_t k);

        //! Returns w-value at given data point.
        const double& w(size_t i, size_t j, size_t k) const;

        //! Returns interpolated value at cell center.
        Vector3D valueAtCellCenter(size_t i, size_t j, size_t k) const;

        //! Returns divergence at cell-center location.
        double divergenceAtCellCenter(size_t i, size_t j, size_t k) const;

        //! Returns curl at cell-center location.
        Vector3D curlAtCellCenter(size_t i, size_t j, size_t k) const;

        //! Returns u data accessor.
        ScalarDataAccessor uAccessor();

        //! Returns read-only u data accessor.
        ConstScalarDataAccessor uConstAccessor() const;

        //! Returns v data accessor.
        ScalarDataAccessor vAccessor();

        //! Returns read-only v data accessor.
        ConstScalarDataAccessor vConstAccessor() const;

        //! Returns w data accessor.
        ScalarDataAccessor wAccessor();

        //! Returns read-only w data accessor.
        ConstScalarDataAccessor wConstAccessor() const;

        //! Returns function object that maps u data point to its actual position.
        DataPositionFunc uPosition() const;

        //! Returns function object that maps v data point to its actual position.
        DataPositionFunc vPosition() const;

        //! Returns function object that maps w data point to its actual position.
        DataPositionFunc wPosition() const;

        //! Returns data size of the u component.
        Size3 uSize() const;

        //! Returns data size of the v component.
        Size3 vSize() const;

        //! Returns data size of the w component.
        Size3 wSize() const;

        //!
        //! \brief Returns u-data position for the grid point at (0, 0, 0).
        //!
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        //!
        Vector3D uOrigin() const;

        //!
        //! \brief Returns v-data position for the grid point at (0, 0, 0).
        //!
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        //!
        Vector3D vOrigin() const;

        //!
        //! \brief Returns w-data position for the grid point at (0, 0, 0).
        //!
        //! Note that this is different from origin() since origin() returns
        //! the lower corner point of the bounding box.
        //!
        Vector3D wOrigin() const;

        //! Fills the grid with given value.
        void fill(const Vector3D& value,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Fills the grid with given function.
        void fill(const std::function<Vector3D(const Vector3D&)>& func,
            ExecutionPolicy policy = ExecutionPolicy::kParallel) override;

        //! Returns the copy of the grid instance.
        std::shared_ptr<VectorGrid3> clone() const override;

        //!
        //! \brief Invokes the given function \p func for each u-data point.
        //!
        //! This function invokes the given function object \p func for each u-data
        //! point in serial manner. The input parameters are i and j indices of a
        //! u-data point. The order of execution is i-first, j-last.
        //!
        void forEachUIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each u-data point
        //! parallelly.
        //!
        //! This function invokes the given function object \p func for each u-data
        //! point in parallel manner. The input parameters are i and j indices of a
        //! u-data point. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachUIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each v-data point.
        //!
        //! This function invokes the given function object \p func for each v-data
        //! point in serial manner. The input parameters are i and j indices of a
        //! v-data point. The order of execution is i-first, j-last.
        //!
        void forEachVIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each v-data point
        //! parallelly.
        //!
        //! This function invokes the given function object \p func for each v-data
        //! point in parallel manner. The input parameters are i and j indices of a
        //! v-data point. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachVIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each w-data point.
        //!
        //! This function invokes the given function object \p func for each w-data
        //! point in serial manner. The input parameters are i and j indices of a
        //! w-data point. The order of execution is i-first, j-last.
        //!
        void forEachWIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        //!
        //! \brief Invokes the given function \p func for each w-data point
        //! parallelly.
        //!
        //! This function invokes the given function object \p func for each w-data
        //! point in parallel manner. The input parameters are i and j indices of a
        //! w-data point. The order of execution can be arbitrary since it's
        //! multi-threaded.
        //!
        void parallelForEachWIndex(
            const std::function<void(size_t, size_t, size_t)>& func) const;

        // VectorField3 implementations

        //! Returns sampled value at given position \p x.
        Vector3D sample(const Vector3D& x) const override;

        //! Returns divergence at given position \p x.
        double divergence(const Vector3D& x) const override;

        //! Returns curl at given position \p x.
        Vector3D curl(const Vector3D& x) const override;

        //!
        //! \brief Returns the sampler function.
        //!
        //! This function returns the data sampler function object. The sampling
        //! function is linear.
        //!
        std::function<Vector3D(const Vector3D&)> sampler() const override;

        //! Returns builder fox FaceCenteredGrid3.
        static Builder builder();

    protected:
        // VectorGrid3 implementations
        void onResize(const Size3& resolution, const Vector3D& gridSpacing,
            const Vector3D& origin, const Vector3D& initialValue) final;

        //! Fetches the data into a continuous linear array.
        void getData(std::vector<double>* data) const override;

        //! Sets the data from a continuous linear array.
        void setData(const std::vector<double>& data) override;

    private:
        Array3<double> _dataU;
        Array3<double> _dataV;
        Array3<double> _dataW;
        Vector3D _dataOriginU;
        Vector3D _dataOriginV;
        Vector3D _dataOriginW;
        LinearArraySampler3<double, double> _uLinearSampler;
        LinearArraySampler3<double, double> _vLinearSampler;
        LinearArraySampler3<double, double> _wLinearSampler;
        std::function<Vector3D(const Vector3D&)> _sampler;

        void resetSampler();
    };

    //! Shared pointer type for the FaceCenteredGrid3.
    typedef std::shared_ptr<FaceCenteredGrid3> FaceCenteredGrid3Ptr;

    //!
    //! \brief Front-end to create CellCenteredScalarGrid3 objects step by step.
    //!
    class FaceCenteredGrid3::Builder final : public VectorGridBuilder3 {
    public:
        //! Returns builder with resolution.
        Builder& withResolution(const Size3& resolution);

        //! Returns builder with resolution.
        Builder& withResolution(size_t resolutionX, size_t resolutionY,
            size_t resolutionZ);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(const Vector3D& gridSpacing);

        //! Returns builder with grid spacing.
        Builder& withGridSpacing(double gridSpacingX, double gridSpacingY,
            double gridSpacingZ);

        //! Returns builder with grid origin.
        Builder& withOrigin(const Vector3D& gridOrigin);

        //! Returns builder with grid origin.
        Builder& withOrigin(double gridOriginX, double gridOriginY,
            double gridOriginZ);

        //! Returns builder with initial value.
        Builder& withInitialValue(const Vector3D& initialVal);

        //! Returns builder with initial value.
        Builder& withInitialValue(double initialValX, double initialValY,
            double initialValZ);

        //! Builds CellCenteredScalarGrid3 instance.
        FaceCenteredGrid3 build() const;

        //! Builds shared pointer of FaceCenteredGrid3 instance.
        FaceCenteredGrid3Ptr makeShared() const;

        //!
        //! \brief Builds shared pointer of FaceCenteredGrid3 instance.
        //!
        //! This is an overriding function that implements VectorGridBuilder3.
        //!
        VectorGrid3Ptr build(const Size3& resolution, const Vector3D& gridSpacing,
            const Vector3D& gridOrigin,
            const Vector3D& initialVal) const override;

    private:
        Size3 _resolution{ 1, 1, 1 };
        Vector3D _gridSpacing{ 1, 1, 1 };
        Vector3D _gridOrigin{ 0, 0, 0 };
        Vector3D _initialVal{ 0, 0, 0 };
    };

    FaceCenteredGrid3::FaceCenteredGrid3()
        : _dataOriginU(0.0, 0.5, 0.5),
        _dataOriginV(0.5, 0.0, 0.5),
        _dataOriginW(0.5, 0.5, 0.0),
        _uLinearSampler(LinearArraySampler3<double, double>(
            _dataU.constAccessor(), Vector3D(1, 1, 1), _dataOriginU)),
        _vLinearSampler(LinearArraySampler3<double, double>(
            _dataV.constAccessor(), Vector3D(1, 1, 1), _dataOriginV)),
        _wLinearSampler(LinearArraySampler3<double, double>(
            _dataW.constAccessor(), Vector3D(1, 1, 1), _dataOriginW)) {}

    FaceCenteredGrid3::FaceCenteredGrid3(size_t resolutionX, size_t resolutionY,
        size_t resolutionZ, double gridSpacingX,
        double gridSpacingY, double gridSpacingZ,
        double originX, double originY,
        double originZ, double initialValueU,
        double initialValueV, double initialValueW)
        : FaceCenteredGrid3(Size3(resolutionX, resolutionY, resolutionZ),
            Vector3D(gridSpacingX, gridSpacingY, gridSpacingZ),
            Vector3D(originX, originY, originZ),
            Vector3D(initialValueU, initialValueV, initialValueW)) {
    }

    FaceCenteredGrid3::FaceCenteredGrid3(const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin,
        const Vector3D& initialValue)
        : _uLinearSampler(LinearArraySampler3<double, double>(
            _dataU.constAccessor(), Vector3D(1, 1, 1), _dataOriginU)),
        _vLinearSampler(LinearArraySampler3<double, double>(
            _dataV.constAccessor(), Vector3D(1, 1, 1), _dataOriginV)),
        _wLinearSampler(LinearArraySampler3<double, double>(
            _dataW.constAccessor(), Vector3D(1, 1, 1), _dataOriginW)) {
        resize(resolution, gridSpacing, origin, initialValue);
    }

    FaceCenteredGrid3::FaceCenteredGrid3(const FaceCenteredGrid3& other)
        : _uLinearSampler(LinearArraySampler3<double, double>(
            _dataU.constAccessor(), Vector3D(1, 1, 1), _dataOriginU)),
        _vLinearSampler(LinearArraySampler3<double, double>(
            _dataV.constAccessor(), Vector3D(1, 1, 1), _dataOriginV)),
        _wLinearSampler(LinearArraySampler3<double, double>(
            _dataW.constAccessor(), Vector3D(1, 1, 1), _dataOriginW)) {
        set(other);
    }

    void FaceCenteredGrid3::swap(Grid3* other) {
        FaceCenteredGrid3* sameType = dynamic_cast<FaceCenteredGrid3*>(other);

        if (sameType != nullptr) {
            swapGrid(sameType);

            _dataU.swap(sameType->_dataU);
            _dataV.swap(sameType->_dataV);
            _dataW.swap(sameType->_dataW);
            std::swap(_dataOriginU, sameType->_dataOriginU);
            std::swap(_dataOriginV, sameType->_dataOriginV);
            std::swap(_dataOriginW, sameType->_dataOriginW);
            std::swap(_uLinearSampler, sameType->_uLinearSampler);
            std::swap(_vLinearSampler, sameType->_vLinearSampler);
            std::swap(_wLinearSampler, sameType->_wLinearSampler);
            std::swap(_sampler, sameType->_sampler);
        }
    }

    void FaceCenteredGrid3::set(const FaceCenteredGrid3& other) {
        setGrid(other);

        _dataU.set(other._dataU);
        _dataV.set(other._dataV);
        _dataW.set(other._dataW);
        _dataOriginU = other._dataOriginU;
        _dataOriginV = other._dataOriginV;
        _dataOriginW = other._dataOriginW;

        resetSampler();
    }

    FaceCenteredGrid3& FaceCenteredGrid3::operator=(
        const FaceCenteredGrid3& other) {
        set(other);
        return *this;
    }

    double& FaceCenteredGrid3::u(size_t i, size_t j, size_t k) {
        return _dataU(i, j, k);
    }

    const double& FaceCenteredGrid3::u(size_t i, size_t j, size_t k) const {
        return _dataU(i, j, k);
    }

    double& FaceCenteredGrid3::v(size_t i, size_t j, size_t k) {
        return _dataV(i, j, k);
    }

    const double& FaceCenteredGrid3::v(size_t i, size_t j, size_t k) const {
        return _dataV(i, j, k);
    }

    double& FaceCenteredGrid3::w(size_t i, size_t j, size_t k) {
        return _dataW(i, j, k);
    }

    const double& FaceCenteredGrid3::w(size_t i, size_t j, size_t k) const {
        return _dataW(i, j, k);
    }

    Vector3D FaceCenteredGrid3::valueAtCellCenter(size_t i, size_t j,
        size_t k) const {
        JET_ASSERT(i < resolution().x&& j < resolution().y&& k < resolution().z);

        return 0.5 * Vector3D(_dataU(i, j, k) + _dataU(i + 1, j, k),
            _dataV(i, j, k) + _dataV(i, j + 1, k),
            _dataW(i, j, k) + _dataW(i, j, k + 1));
    }

    double FaceCenteredGrid3::divergenceAtCellCenter(size_t i, size_t j,
        size_t k) const {
        JET_ASSERT(i < resolution().x&& j < resolution().y&& k < resolution().z);

        const Vector3D& gs = gridSpacing();

        double leftU = _dataU(i, j, k);
        double rightU = _dataU(i + 1, j, k);
        double bottomV = _dataV(i, j, k);
        double topV = _dataV(i, j + 1, k);
        double backW = _dataW(i, j, k);
        double frontW = _dataW(i, j, k + 1);

        return (rightU - leftU) / gs.x + (topV - bottomV) / gs.y +
            (frontW - backW) / gs.z;
    }

    Vector3D FaceCenteredGrid3::curlAtCellCenter(size_t i, size_t j,
        size_t k) const {
        const Size3& res = resolution();
        const Vector3D& gs = gridSpacing();

        JET_ASSERT(i < res.x&& j < res.y&& k < res.z);

        Vector3D left = valueAtCellCenter((i > 0) ? i - 1 : i, j, k);
        Vector3D right = valueAtCellCenter((i + 1 < res.x) ? i + 1 : i, j, k);
        Vector3D down = valueAtCellCenter(i, (j > 0) ? j - 1 : j, k);
        Vector3D up = valueAtCellCenter(i, (j + 1 < res.y) ? j + 1 : j, k);
        Vector3D back = valueAtCellCenter(i, j, (k > 0) ? k - 1 : k);
        Vector3D front = valueAtCellCenter(i, j, (k + 1 < res.z) ? k + 1 : k);

        double Fx_ym = down.x;
        double Fx_yp = up.x;
        double Fx_zm = back.x;
        double Fx_zp = front.x;

        double Fy_xm = left.y;
        double Fy_xp = right.y;
        double Fy_zm = back.y;
        double Fy_zp = front.y;

        double Fz_xm = left.z;
        double Fz_xp = right.z;
        double Fz_ym = down.z;
        double Fz_yp = up.z;

        return Vector3D(
            0.5 * (Fz_yp - Fz_ym) / gs.y - 0.5 * (Fy_zp - Fy_zm) / gs.z,
            0.5 * (Fx_zp - Fx_zm) / gs.z - 0.5 * (Fz_xp - Fz_xm) / gs.x,
            0.5 * (Fy_xp - Fy_xm) / gs.x - 0.5 * (Fx_yp - Fx_ym) / gs.y);
    }

    FaceCenteredGrid3::ScalarDataAccessor FaceCenteredGrid3::uAccessor() {
        return _dataU.accessor();
    }

    FaceCenteredGrid3::ConstScalarDataAccessor FaceCenteredGrid3::uConstAccessor()
        const {
        return _dataU.constAccessor();
    }

    FaceCenteredGrid3::ScalarDataAccessor FaceCenteredGrid3::vAccessor() {
        return _dataV.accessor();
    }

    FaceCenteredGrid3::ConstScalarDataAccessor FaceCenteredGrid3::vConstAccessor()
        const {
        return _dataV.constAccessor();
    }

    FaceCenteredGrid3::ScalarDataAccessor FaceCenteredGrid3::wAccessor() {
        return _dataW.accessor();
    }

    FaceCenteredGrid3::ConstScalarDataAccessor FaceCenteredGrid3::wConstAccessor()
        const {
        return _dataW.constAccessor();
    }

    VectorGrid3::DataPositionFunc FaceCenteredGrid3::uPosition() const {
        Vector3D h = gridSpacing();

        return [this, h](size_t i, size_t j, size_t k) -> Vector3D {
            return _dataOriginU + h * Vector3D({ i, j, k });
        };
    }

    VectorGrid3::DataPositionFunc FaceCenteredGrid3::vPosition() const {
        Vector3D h = gridSpacing();

        return [this, h](size_t i, size_t j, size_t k) -> Vector3D {
            return _dataOriginV + h * Vector3D({ i, j, k });
        };
    }

    VectorGrid3::DataPositionFunc FaceCenteredGrid3::wPosition() const {
        Vector3D h = gridSpacing();

        return [this, h](size_t i, size_t j, size_t k) -> Vector3D {
            return _dataOriginW + h * Vector3D({ i, j, k });
        };
    }

    Size3 FaceCenteredGrid3::uSize() const { return _dataU.size(); }

    Size3 FaceCenteredGrid3::vSize() const { return _dataV.size(); }

    Size3 FaceCenteredGrid3::wSize() const { return _dataW.size(); }

    Vector3D FaceCenteredGrid3::uOrigin() const { return _dataOriginU; }

    Vector3D FaceCenteredGrid3::vOrigin() const { return _dataOriginV; }

    Vector3D FaceCenteredGrid3::wOrigin() const { return _dataOriginW; }

    void FaceCenteredGrid3::fill(const Vector3D& value, ExecutionPolicy policy) {
        parallelFor(kZeroSize, _dataU.width(), kZeroSize, _dataU.height(),
            kZeroSize, _dataU.depth(),
            [this, value](size_t i, size_t j, size_t k) {
                _dataU(i, j, k) = value.x;
            },
            policy);

        parallelFor(kZeroSize, _dataV.width(), kZeroSize, _dataV.height(),
            kZeroSize, _dataV.depth(),
            [this, value](size_t i, size_t j, size_t k) {
                _dataV(i, j, k) = value.y;
            },
            policy);

        parallelFor(kZeroSize, _dataW.width(), kZeroSize, _dataW.height(),
            kZeroSize, _dataW.depth(),
            [this, value](size_t i, size_t j, size_t k) {
                _dataW(i, j, k) = value.z;
            },
            policy);
    }

    void FaceCenteredGrid3::fill(
        const std::function<Vector3D(const Vector3D&)>& func,
        ExecutionPolicy policy) {
        DataPositionFunc uPos = uPosition();
        parallelFor(kZeroSize, _dataU.width(), kZeroSize, _dataU.height(),
            kZeroSize, _dataU.depth(),
            [this, &func, &uPos](size_t i, size_t j, size_t k) {
                _dataU(i, j, k) = func(uPos(i, j, k)).x;
            },
            policy);
        DataPositionFunc vPos = vPosition();
        parallelFor(kZeroSize, _dataV.width(), kZeroSize, _dataV.height(),
            kZeroSize, _dataV.depth(),
            [this, &func, &vPos](size_t i, size_t j, size_t k) {
                _dataV(i, j, k) = func(vPos(i, j, k)).y;
            },
            policy);
        DataPositionFunc wPos = wPosition();
        parallelFor(kZeroSize, _dataW.width(), kZeroSize, _dataW.height(),
            kZeroSize, _dataW.depth(),
            [this, &func, &wPos](size_t i, size_t j, size_t k) {
                _dataW(i, j, k) = func(wPos(i, j, k)).z;
            },
            policy);
    }

    std::shared_ptr<VectorGrid3> FaceCenteredGrid3::clone() const {
        return CLONE_W_CUSTOM_DELETER(FaceCenteredGrid3);
    }

    void FaceCenteredGrid3::forEachUIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataU.forEachIndex(func);
    }

    void FaceCenteredGrid3::parallelForEachUIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataU.parallelForEachIndex(func);
    }

    void FaceCenteredGrid3::forEachVIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataV.forEachIndex(func);
    }

    void FaceCenteredGrid3::parallelForEachVIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataV.parallelForEachIndex(func);
    }

    void FaceCenteredGrid3::forEachWIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataW.forEachIndex(func);
    }

    void FaceCenteredGrid3::parallelForEachWIndex(
        const std::function<void(size_t, size_t, size_t)>& func) const {
        _dataW.parallelForEachIndex(func);
    }

    Vector3D FaceCenteredGrid3::sample(const Vector3D& x) const {
        return _sampler(x);
    }

    std::function<Vector3D(const Vector3D&)> FaceCenteredGrid3::sampler() const {
        return _sampler;
    }

    double FaceCenteredGrid3::divergence(const Vector3D& x) const {
        Size3 res = resolution();
        ssize_t i, j, k;
        double fx, fy, fz;
        Vector3D cellCenterOrigin = origin() + 0.5 * gridSpacing();

        Vector3D normalizedX = (x - cellCenterOrigin) / gridSpacing();

        getBarycentric(normalizedX.x, 0, static_cast<ssize_t>(res.x) - 1, &i, &fx);
        getBarycentric(normalizedX.y, 0, static_cast<ssize_t>(res.y) - 1, &j, &fy);
        getBarycentric(normalizedX.z, 0, static_cast<ssize_t>(res.z) - 1, &k, &fz);

        std::array<Point3UI, 8> indices;
        std::array<double, 8> weights;

        indices[0] = Point3UI(i, j, k);
        indices[1] = Point3UI(i + 1, j, k);
        indices[2] = Point3UI(i, j + 1, k);
        indices[3] = Point3UI(i + 1, j + 1, k);
        indices[4] = Point3UI(i, j, k + 1);
        indices[5] = Point3UI(i + 1, j, k + 1);
        indices[6] = Point3UI(i, j + 1, k + 1);
        indices[7] = Point3UI(i + 1, j + 1, k + 1);

        weights[0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
        weights[1] = fx * (1.0 - fy) * (1.0 - fz);
        weights[2] = (1.0 - fx) * fy * (1.0 - fz);
        weights[3] = fx * fy * (1.0 - fz);
        weights[4] = (1.0 - fx) * (1.0 - fy) * fz;
        weights[5] = fx * (1.0 - fy) * fz;
        weights[6] = (1.0 - fx) * fy * fz;
        weights[7] = fx * fy * fz;

        double result = 0.0;

        for (int n = 0; n < 8; ++n) {
            result += weights[n] * divergenceAtCellCenter(
                indices[n].x, indices[n].y, indices[n].z);
        }

        return result;
    }

    Vector3D FaceCenteredGrid3::curl(const Vector3D& x) const {
        Size3 res = resolution();
        ssize_t i, j, k;
        double fx, fy, fz;
        Vector3D cellCenterOrigin = origin() + 0.5 * gridSpacing();

        Vector3D normalizedX = (x - cellCenterOrigin) / gridSpacing();

        getBarycentric(normalizedX.x, 0, static_cast<ssize_t>(res.x) - 1, &i, &fx);
        getBarycentric(normalizedX.y, 0, static_cast<ssize_t>(res.y) - 1, &j, &fy);
        getBarycentric(normalizedX.z, 0, static_cast<ssize_t>(res.z) - 1, &k, &fz);

        std::array<Point3UI, 8> indices;
        std::array<double, 8> weights;

        indices[0] = Point3UI(i, j, k);
        indices[1] = Point3UI(i + 1, j, k);
        indices[2] = Point3UI(i, j + 1, k);
        indices[3] = Point3UI(i + 1, j + 1, k);
        indices[4] = Point3UI(i, j, k + 1);
        indices[5] = Point3UI(i + 1, j, k + 1);
        indices[6] = Point3UI(i, j + 1, k + 1);
        indices[7] = Point3UI(i + 1, j + 1, k + 1);

        weights[0] = (1.0 - fx) * (1.0 - fy) * (1.0 - fz);
        weights[1] = fx * (1.0 - fy) * (1.0 - fz);
        weights[2] = (1.0 - fx) * fy * (1.0 - fz);
        weights[3] = fx * fy * (1.0 - fz);
        weights[4] = (1.0 - fx) * (1.0 - fy) * fz;
        weights[5] = fx * (1.0 - fy) * fz;
        weights[6] = (1.0 - fx) * fy * fz;
        weights[7] = fx * fy * fz;

        Vector3D result;

        for (int n = 0; n < 8; ++n) {
            result += weights[n] *
                curlAtCellCenter(indices[n].x, indices[n].y, indices[n].z);
        }

        return result;
    }

    void FaceCenteredGrid3::onResize(const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin,
        const Vector3D& initialValue) {
        if (resolution != Size3(0, 0, 0)) {
            _dataU.resize(resolution + Size3(1, 0, 0), initialValue.x);
            _dataV.resize(resolution + Size3(0, 1, 0), initialValue.y);
            _dataW.resize(resolution + Size3(0, 0, 1), initialValue.z);
        }
        else {
            _dataU.resize(Size3(0, 0, 0));
            _dataV.resize(Size3(0, 0, 0));
            _dataW.resize(Size3(0, 0, 0));
        }
        _dataOriginU = origin + 0.5 * Vector3D(0.0, gridSpacing.y, gridSpacing.z);
        _dataOriginV = origin + 0.5 * Vector3D(gridSpacing.x, 0.0, gridSpacing.z);
        _dataOriginW = origin + 0.5 * Vector3D(gridSpacing.x, gridSpacing.y, 0.0);

        resetSampler();
    }

    void FaceCenteredGrid3::resetSampler() {
        LinearArraySampler3<double, double> uSampler(_dataU.constAccessor(),
            gridSpacing(), _dataOriginU);
        LinearArraySampler3<double, double> vSampler(_dataV.constAccessor(),
            gridSpacing(), _dataOriginV);
        LinearArraySampler3<double, double> wSampler(_dataW.constAccessor(),
            gridSpacing(), _dataOriginW);

        _uLinearSampler = uSampler;
        _vLinearSampler = vSampler;
        _wLinearSampler = wSampler;

        _sampler = [uSampler, vSampler, wSampler](const Vector3D& x) -> Vector3D {
            double u = uSampler(x);
            double v = vSampler(x);
            double w = wSampler(x);
            return Vector3D(u, v, w);
        };
    }

    FaceCenteredGrid3::Builder FaceCenteredGrid3::builder() { return Builder(); }

    void FaceCenteredGrid3::getData(std::vector<double>* data) const {
        size_t size = uSize().x * uSize().y * uSize().z +
            vSize().x * vSize().y * vSize().z +
            wSize().x * wSize().y * wSize().z;
        data->resize(size);
        size_t cnt = 0;
        _dataU.forEach([&](double value) { (*data)[cnt++] = value; });
        _dataV.forEach([&](double value) { (*data)[cnt++] = value; });
        _dataW.forEach([&](double value) { (*data)[cnt++] = value; });
    }

    void FaceCenteredGrid3::setData(const std::vector<double>& data) {
        JET_ASSERT(uSize().x * uSize().y * uSize().z +
            vSize().x * vSize().y * vSize().z +
            wSize().x * wSize().y * wSize().z ==
            data.size());

        size_t cnt = 0;
        _dataU.forEachIndex(
            [&](size_t i, size_t j, size_t k) { _dataU(i, j, k) = data[cnt++]; });
        _dataV.forEachIndex(
            [&](size_t i, size_t j, size_t k) { _dataV(i, j, k) = data[cnt++]; });
        _dataW.forEachIndex(
            [&](size_t i, size_t j, size_t k) { _dataW(i, j, k) = data[cnt++]; });
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withResolution(
        const Size3& resolution) {
        _resolution = resolution;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withResolution(
        size_t resolutionX, size_t resolutionY, size_t resolutionZ) {
        _resolution.x = resolutionX;
        _resolution.y = resolutionY;
        _resolution.z = resolutionZ;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withGridSpacing(
        const Vector3D& gridSpacing) {
        _gridSpacing = gridSpacing;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withGridSpacing(
        double gridSpacingX, double gridSpacingY, double gridSpacingZ) {
        _gridSpacing.x = gridSpacingX;
        _gridSpacing.y = gridSpacingY;
        _gridSpacing.z = gridSpacingZ;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withOrigin(
        const Vector3D& gridOrigin) {
        _gridOrigin = gridOrigin;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withOrigin(
        double gridOriginX, double gridOriginY, double gridOriginZ) {
        _gridOrigin.x = gridOriginX;
        _gridOrigin.y = gridOriginY;
        _gridOrigin.z = gridOriginZ;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withInitialValue(
        const Vector3D& initialVal) {
        _initialVal = initialVal;
        return *this;
    }

    FaceCenteredGrid3::Builder& FaceCenteredGrid3::Builder::withInitialValue(
        double initialValX, double initialValY, double initialValZ) {
        _initialVal.x = initialValX;
        _initialVal.y = initialValY;
        _initialVal.z = initialValZ;
        return *this;
    }

    FaceCenteredGrid3 FaceCenteredGrid3::Builder::build() const {
        return FaceCenteredGrid3(_resolution, _gridSpacing, _gridOrigin,
            _initialVal);
    }

    FaceCenteredGrid3Ptr FaceCenteredGrid3::Builder::makeShared() const {
        return std::shared_ptr<FaceCenteredGrid3>(
            new FaceCenteredGrid3(_resolution, _gridSpacing, _gridOrigin,
                _initialVal),
            [](FaceCenteredGrid3* obj) { delete obj; });
    }

    VectorGrid3Ptr FaceCenteredGrid3::Builder::build(
        const Size3& resolution, const Vector3D& gridSpacing,
        const Vector3D& gridOrigin, const Vector3D& initialVal) const {
        return std::shared_ptr<FaceCenteredGrid3>(
            new FaceCenteredGrid3(resolution, gridSpacing, gridOrigin, initialVal),
            [](FaceCenteredGrid3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_FACE_CENTERED_GRID3_H_