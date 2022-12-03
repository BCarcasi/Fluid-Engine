#ifndef INCLUDE_JET_CUBIC_SEMI_LAGRANGIAN3_H_
#define INCLUDE_JET_CUBIC_SEMI_LAGRANGIAN3_H_

#include "semi_lagrangian3.h"
#include "pch.h"
#include "array_samplers3.h"

namespace jet {

    //!
    //! \brief Implementation of 3-D cubic semi-Lagrangian advection solver.
    //!
    //! This class implements 3rd-order cubic 3-D semi-Lagrangian advection solver.
    //!
    class CubicSemiLagrangian3 final : public SemiLagrangian3 {
    public:
        CubicSemiLagrangian3();

    protected:
        //!
        //! \brief Returns spatial interpolation function object for given scalar
        //! grid.
        //!
        //! This function overrides the original function with cubic interpolation.
        //!
        std::function<double(const Vector3D&)> getScalarSamplerFunc(
            const ScalarGrid3& source) const override;

        //!
        //! \brief Returns spatial interpolation function object for given
        //! collocated vector grid.
        //!
        //! This function overrides the original function with cubic interpolation.
        //!
        std::function<Vector3D(const Vector3D&)> getVectorSamplerFunc(
            const CollocatedVectorGrid3& source) const override;

        //!
        //! \brief Returns spatial interpolation function object for given
        //! face-centered vector grid.
        //!
        //! This function overrides the original function with cubic interpolation.
        //!
        std::function<Vector3D(const Vector3D&)> getVectorSamplerFunc(
            const FaceCenteredGrid3& source) const override;
    };

    typedef std::shared_ptr<CubicSemiLagrangian3> CubicSemiLagrangian3Ptr;


    CubicSemiLagrangian3::CubicSemiLagrangian3() {
    }

    std::function<double(const Vector3D&)>
        CubicSemiLagrangian3::getScalarSamplerFunc(const ScalarGrid3& source) const {
        auto sourceSampler = CubicArraySampler3<double, double>(
            source.constDataAccessor(),
            source.gridSpacing(),
            source.dataOrigin());
        return sourceSampler.functor();
    }

    std::function<Vector3D(const Vector3D&)>
        CubicSemiLagrangian3::getVectorSamplerFunc(
            const CollocatedVectorGrid3& source) const {
        auto sourceSampler = CubicArraySampler3<Vector3D, double>(
            source.constDataAccessor(),
            source.gridSpacing(),
            source.dataOrigin());
        return sourceSampler.functor();
    }

    std::function<Vector3D(const Vector3D&)>
        CubicSemiLagrangian3::getVectorSamplerFunc(
            const FaceCenteredGrid3& source) const {
        auto uSourceSampler = CubicArraySampler3<double, double>(
            source.uConstAccessor(),
            source.gridSpacing(),
            source.uOrigin());
        auto vSourceSampler = CubicArraySampler3<double, double>(
            source.vConstAccessor(),
            source.gridSpacing(),
            source.vOrigin());
        auto wSourceSampler = CubicArraySampler3<double, double>(
            source.wConstAccessor(),
            source.gridSpacing(),
            source.wOrigin());
        return
            [uSourceSampler, vSourceSampler, wSourceSampler](const Vector3D& x) {
            return Vector3D(
                uSourceSampler(x), vSourceSampler(x), wSourceSampler(x));
        };
    }

}  // namespace jet

#endif  // INCLUDE_JET_CUBIC_SEMI_LAGRANGIAN3_H_