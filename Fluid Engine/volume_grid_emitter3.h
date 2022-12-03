#ifndef INCLUDE_JET_VOLUME_GRID_EMITTER3_H_
#define INCLUDE_JET_VOLUME_GRID_EMITTER3_H_

#include "grid_emitter3.h"
#include "scalar_grid3.h"
#include "vector_grid3.h"

#include <tuple>
#include <vector>

#include "pch.h"

#include "collocated_vector_grid3.h"
#include "face_centered_grid3.h"
#include "level_set_utils.h"
#include "surface_to_implicit3.h"

#include <algorithm>

namespace jet {

    //!
    //! \brief 3-D grid-based volumetric emitter.
    //!
    class VolumeGridEmitter3 final : public GridEmitter3 {
    public:
        class Builder;

        //! Maps to a scalar value for given signed-dist, location, and old value.
        typedef std::function<double(double, const Vector3D&, double)>
            ScalarMapper;

        //! Maps to a vector value for given signed-dist, location, and old value.
        typedef std::function<Vector3D(double, const Vector3D&, const Vector3D&)>
            VectorMapper;

        //!
        //! \brief      Constructs an emitter with a source and is-one-shot flag.
        //!
        //! \param[in]  sourceRegion    Emitting region given by the SDF.
        //! \param[in]  isOneShot       True if emitter gets disabled after one shot.
        //!
        explicit VolumeGridEmitter3(
            const ImplicitSurface3Ptr& sourceRegion,
            bool isOneShot = true);

        //! Destructor.
        virtual ~VolumeGridEmitter3();

        //! Adds signed-distance target to the scalar grid.
        void addSignedDistanceTarget(const ScalarGrid3Ptr& scalarGridTarget);

        //!
        //! \brief      Adds step function target to the scalar grid.
        //!
        //! \param[in]  scalarGridTarget The scalar grid target.
        //! \param[in]  minValue         The minimum value of the step function.
        //! \param[in]  maxValue         The maximum value of the step function.
        //!
        void addStepFunctionTarget(
            const ScalarGrid3Ptr& scalarGridTarget,
            double minValue,
            double maxValue);

        //!
        //! \brief      Adds a scalar grid target.
        //!
        //! This function adds a custom target to the emitter. The first parameter
        //! defines which grid should it write to. The second parameter,
        //! \p customMapper, defines how to map signed-distance field from the
        //! volume geometry and location of the point to the final scalar value that
        //! is going to be written to the target grid. The third parameter defines
        //! how to blend the old value from the target grid and the new value from
        //! the mapper function.
        //!
        //! \param[in]  scalarGridTarget The scalar grid target
        //! \param[in]  customMapper     The custom mapper.
        //!
        void addTarget(
            const ScalarGrid3Ptr& scalarGridTarget,
            const ScalarMapper& customMapper);

        //!
        //! \brief      Adds a vector grid target.
        //!
        //! This function adds a custom target to the emitter. The first parameter
        //! defines which grid should it write to. The second parameter,
        //! \p customMapper, defines how to map sigend-distance field from the
        //! volume geometry and location of the point to the final vector value that
        //! is going to be written to the target grid. The third parameter defines
        //! how to blend the old value from the target grid and the new value from
        //! the mapper function.
        //!
        //! \param[in]  scalarGridTarget The vector grid target
        //! \param[in]  customMapper     The custom mapper.
        //!
        void addTarget(
            const VectorGrid3Ptr& vectorGridTarget,
            const VectorMapper& customMapper);

        //! Returns implicit surface which defines the source region.
        const ImplicitSurface3Ptr& sourceRegion() const;

        //! Returns true if this emits only once.
        bool isOneShot() const;

        //! Returns builder fox VolumeGridEmitter3.
        static Builder builder();

    private:
        typedef std::tuple<ScalarGrid3Ptr, ScalarMapper> ScalarTarget;
        typedef std::tuple<VectorGrid3Ptr, VectorMapper> VectorTarget;

        ImplicitSurface3Ptr _sourceRegion;
        bool _isOneShot = true;
        bool _hasEmitted = false;
        std::vector<ScalarTarget> _customScalarTargets;
        std::vector<VectorTarget> _customVectorTargets;

        void onUpdate(
            double currentTimeInSeconds,
            double timeIntervalInSeconds) override;

        void emit();
    };

    //! Shared pointer type for the VolumeGridEmitter3.
    typedef std::shared_ptr<VolumeGridEmitter3> VolumeGridEmitter3Ptr;


    //!
    //! \brief Front-end to create VolumeGridEmitter3 objects step by step.
    //!
    class VolumeGridEmitter3::Builder final {
    public:
        //! Returns builder with surface defining source region.
        Builder& withSourceRegion(const Surface3Ptr& sourceRegion);

        //! Returns builder with one-shot flag.
        Builder& withIsOneShot(bool isOneShot);

        //! Builds VolumeGridEmitter3.
        VolumeGridEmitter3 build() const;

        //! Builds shared pointer of VolumeGridEmitter3 instance.
        VolumeGridEmitter3Ptr makeShared() const;

    private:
        ImplicitSurface3Ptr _sourceRegion;
        bool _isOneShot = true;
    };

    VolumeGridEmitter3::VolumeGridEmitter3(const ImplicitSurface3Ptr& sourceRegion,
        bool isOneShot)
        : _sourceRegion(sourceRegion), _isOneShot(isOneShot) {}

    VolumeGridEmitter3::~VolumeGridEmitter3() {}

    void VolumeGridEmitter3::addSignedDistanceTarget(
        const ScalarGrid3Ptr& scalarGridTarget) {
        auto mapper = [](double sdf, const Vector3D&, double oldVal) {
            return std::min(oldVal, sdf);
        };
        addTarget(scalarGridTarget, mapper);
    }

    void VolumeGridEmitter3::addStepFunctionTarget(
        const ScalarGrid3Ptr& scalarGridTarget, double minValue, double maxValue) {
        double smoothingWidth = scalarGridTarget->gridSpacing().min();
        auto mapper = [minValue, maxValue, smoothingWidth, scalarGridTarget](
            double sdf, const Vector3D&, double oldVal) {
                double step = 1.0 - smearedHeavisideSdf(sdf / smoothingWidth);
                return std::max(oldVal, (maxValue - minValue) * step + minValue);
        };
        addTarget(scalarGridTarget, mapper);
    }

    void VolumeGridEmitter3::addTarget(const ScalarGrid3Ptr& scalarGridTarget,
        const ScalarMapper& customMapper) {
        _customScalarTargets.emplace_back(scalarGridTarget, customMapper);
    }

    void VolumeGridEmitter3::addTarget(const VectorGrid3Ptr& vectorGridTarget,
        const VectorMapper& customMapper) {
        _customVectorTargets.emplace_back(vectorGridTarget, customMapper);
    }

    const ImplicitSurface3Ptr& VolumeGridEmitter3::sourceRegion() const {
        return _sourceRegion;
    }

    bool VolumeGridEmitter3::isOneShot() const { return _isOneShot; }

    void VolumeGridEmitter3::onUpdate(double currentTimeInSeconds,
        double timeIntervalInSeconds) {
        UNUSED_VARIABLE(currentTimeInSeconds);
        UNUSED_VARIABLE(timeIntervalInSeconds);

        if (!isEnabled()) {
            return;
        }

        emit();

        if (_isOneShot) {
            setIsEnabled(false);
        }

        _hasEmitted = true;
    }

    void VolumeGridEmitter3::emit() {
        if (!_sourceRegion) {
            return;
        }

        _sourceRegion->updateQueryEngine();

        for (const auto& target : _customScalarTargets) {
            const auto& grid = std::get<0>(target);
            const auto& mapper = std::get<1>(target);

            auto pos = grid->dataPosition();
            grid->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
                Vector3D gx = pos(i, j, k);
                double sdf = sourceRegion()->signedDistance(gx);
                (*grid)(i, j, k) = mapper(sdf, gx, (*grid)(i, j, k));
                });
        }

        for (const auto& target : _customVectorTargets) {
            const auto& grid = std::get<0>(target);
            const auto& mapper = std::get<1>(target);

            CollocatedVectorGrid3Ptr collocated =
                std::dynamic_pointer_cast<CollocatedVectorGrid3>(grid);
            if (collocated != nullptr) {
                auto pos = collocated->dataPosition();
                collocated->parallelForEachDataPointIndex(
                    [&](size_t i, size_t j, size_t k) {
                        Vector3D gx = pos(i, j, k);
                        double sdf = sourceRegion()->signedDistance(gx);
                        if (isInsideSdf(sdf)) {
                            (*collocated)(i, j, k) =
                                mapper(sdf, gx, (*collocated)(i, j, k));
                        }
                    });
                continue;
            }

            FaceCenteredGrid3Ptr faceCentered =
                std::dynamic_pointer_cast<FaceCenteredGrid3>(grid);
            if (faceCentered != nullptr) {
                auto uPos = faceCentered->uPosition();
                auto vPos = faceCentered->vPosition();
                auto wPos = faceCentered->wPosition();

                faceCentered->parallelForEachUIndex(
                    [&](size_t i, size_t j, size_t k) {
                        Vector3D gx = uPos(i, j, k);
                        double sdf = sourceRegion()->signedDistance(gx);
                        Vector3D oldVal = faceCentered->sample(gx);
                        Vector3D newVal = mapper(sdf, gx, oldVal);
                        faceCentered->u(i, j, k) = newVal.x;
                    });
                faceCentered->parallelForEachVIndex(
                    [&](size_t i, size_t j, size_t k) {
                        Vector3D gx = vPos(i, j, k);
                        double sdf = sourceRegion()->signedDistance(gx);
                        Vector3D oldVal = faceCentered->sample(gx);
                        Vector3D newVal = mapper(sdf, gx, oldVal);
                        faceCentered->v(i, j, k) = newVal.y;
                    });
                faceCentered->parallelForEachWIndex(
                    [&](size_t i, size_t j, size_t k) {
                        Vector3D gx = wPos(i, j, k);
                        double sdf = sourceRegion()->signedDistance(gx);
                        Vector3D oldVal = faceCentered->sample(gx);
                        Vector3D newVal = mapper(sdf, gx, oldVal);
                        faceCentered->w(i, j, k) = newVal.z;
                    });
                continue;
            }
        }
    }

    VolumeGridEmitter3::Builder VolumeGridEmitter3::builder() { return Builder(); }

    VolumeGridEmitter3::Builder& VolumeGridEmitter3::Builder::withSourceRegion(
        const Surface3Ptr& sourceRegion) {
        auto implicit = std::dynamic_pointer_cast<ImplicitSurface3>(sourceRegion);
        if (implicit != nullptr) {
            _sourceRegion = implicit;
        }
        else {
            _sourceRegion = std::make_shared<SurfaceToImplicit3>(sourceRegion);
        }
        return *this;
    }

    VolumeGridEmitter3::Builder& VolumeGridEmitter3::Builder::withIsOneShot(
        bool isOneShot) {
        _isOneShot = isOneShot;
        return *this;
    }

    VolumeGridEmitter3 VolumeGridEmitter3::Builder::build() const {
        return VolumeGridEmitter3(_sourceRegion, _isOneShot);
    }

    VolumeGridEmitter3Ptr VolumeGridEmitter3::Builder::makeShared() const {
        return std::shared_ptr<VolumeGridEmitter3>(
            new VolumeGridEmitter3(_sourceRegion, _isOneShot),
            [](VolumeGridEmitter3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_VOLUME_GRID_EMITTER3_H_