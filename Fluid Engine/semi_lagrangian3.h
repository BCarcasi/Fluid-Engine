#ifndef INCLUDE_JET_SEMI_LAGRANGIAN3_H_
#define INCLUDE_JET_SEMI_LAGRANGIAN3_H_

#include "advection_solver3.h"
#include <limits>
#include "pch.h"
#include "array_samplers3.h"
#include "parallel.h"
#include <algorithm>

namespace jet {

    //!
    //! \brief Implementation of 3-D semi-Lagrangian advection solver.
    //!
    //! This class implements 3-D semi-Lagrangian advection solver. By default, the
    //! class implements 1st-order (linear) algorithm for the spatial interpolation.
    //! For the back-tracing, this class uses 2nd-order mid-point rule with adaptive
    //! time-stepping (CFL <= 1).
    //! To extend the class using higher-order spatial interpolation, the inheriting
    //! classes can override SemiLagrangian2::getScalarSamplerFunc and
    //! SemiLagrangian2::getVectorSamplerFunc. See CubicSemiLagrangian2 for example.
    //!
    class SemiLagrangian3 : public AdvectionSolver3 {
    public:
        SemiLagrangian3();

        virtual ~SemiLagrangian3();

        //!
        //! \brief Computes semi-Langian for given scalar grid.
        //!
        //! This function computes semi-Lagrangian method to solve advection
        //! equation for given scalar field \p input and underlying vector field
        //! \p flow that carries the input field. The solution after solving the
        //! equation for given time-step \p dt should be stored in scalar field
        //! \p output. The boundary interface is given by a signed-distance field.
        //! The field is negative inside the boundary. By default, a constant field
        //! with max double value (kMaxD) is used, meaning no boundary.
        //!
        //! \param input Input scalar grid.
        //! \param flow Vector field that advects the input field.
        //! \param dt Time-step for the advection.
        //! \param output Output scalar grid.
        //! \param boundarySdf Boundary interface defined by signed-distance
        //!     field.
        //!
        void advect(const ScalarGrid3& input, const VectorField3& flow, double dt,
            ScalarGrid3* output,
            const ScalarField3& boundarySdf = ConstantScalarField3(
                std::numeric_limits<double>::max())) final;

        //!
        //! \brief Computes semi-Langian for given collocated vector grid.
        //!
        //! This function computes semi-Lagrangian method to solve advection
        //! equation for given collocated vector grid \p input and underlying vector
        //! field \p flow that carries the input field. The solution after solving
        //! the equation for given time-step \p dt should be stored in scalar field
        //! \p output. The boundary interface is given by a signed-distance field.
        //! The field is negative inside the boundary. By default, a constant field
        //! with max double value (kMaxD) is used, meaning no boundary.
        //!
        //! \param input Input vector grid.
        //! \param flow Vector field that advects the input field.
        //! \param dt Time-step for the advection.
        //! \param output Output vector grid.
        //! \param boundarySdf Boundary interface defined by signed-distance
        //!     field.
        //!
        void advect(const CollocatedVectorGrid3& input, const VectorField3& flow,
            double dt, CollocatedVectorGrid3* output,
            const ScalarField3& boundarySdf = ConstantScalarField3(
                std::numeric_limits<double>::max())) final;

        //!
        //! \brief Computes semi-Langian for given face-centered vector grid.
        //!
        //! This function computes semi-Lagrangian method to solve advection
        //! equation for given face-centered vector grid \p input and underlying
        //! vector field \p flow that carries the input field. The solution after
        //! solving the equation for given time-step \p dt should be stored in
        //! vector field \p output. The boundary interface is given by a
        //! signed-distance field. The field is negative inside the boundary. By
        //! default, a constant field with max double value (kMaxD) is used, meaning
        //! no boundary.
        //!
        //! \param input Input vector grid.
        //! \param flow Vector field that advects the input field.
        //! \param dt Time-step for the advection.
        //! \param output Output vector grid.
        //! \param boundarySdf Boundary interface defined by signed-distance
        //!     field.
        //!
        void advect(const FaceCenteredGrid3& input, const VectorField3& flow,
            double dt, FaceCenteredGrid3* output,
            const ScalarField3& boundarySdf = ConstantScalarField3(
                std::numeric_limits<double>::max())) final;

    protected:
        //!
        //! \brief Returns spatial interpolation function object for given scalar
        //! grid.
        //!
        //! This function returns spatial interpolation function (sampler) for given
        //! scalar grid \p input. By default, this function returns linear
        //! interpolation function. Override this function to have custom
        //! interpolation for semi-Lagrangian process.
        //!
        virtual std::function<double(const Vector3D&)> getScalarSamplerFunc(
            const ScalarGrid3& input) const;

        //!
        //! \brief Returns spatial interpolation function object for given
        //! collocated vector grid.
        //!
        //! This function returns spatial interpolation function (sampler) for given
        //! collocated vector grid \p input. By default, this function returns
        //! linear interpolation function. Override this function to have custom
        //! interpolation for semi-Lagrangian process.
        //!
        virtual std::function<Vector3D(const Vector3D&)> getVectorSamplerFunc(
            const CollocatedVectorGrid3& input) const;

        //!
        //! \brief Returns spatial interpolation function object for given
        //! face-centered vector grid.
        //!
        //! This function returns spatial interpolation function (sampler) for given
        //! face-centered vector grid \p input. By default, this function returns
        //! linear interpolation function. Override this function to have custom
        //! interpolation for semi-Lagrangian process.
        //!
        virtual std::function<Vector3D(const Vector3D&)> getVectorSamplerFunc(
            const FaceCenteredGrid3& input) const;

    private:
        Vector3D backTrace(const VectorField3& flow, double dt, double h,
            const Vector3D& pt0, const ScalarField3& boundarySdf);
    };

    typedef std::shared_ptr<SemiLagrangian3> SemiLagrangian3Ptr;

    SemiLagrangian3::SemiLagrangian3() {
    }

    SemiLagrangian3::~SemiLagrangian3() {
    }

    void SemiLagrangian3::advect(
        const ScalarGrid3& input,
        const VectorField3& flow,
        double dt,
        ScalarGrid3* output,
        const ScalarField3& boundarySdf) {
        auto outputDataPos = output->dataPosition();
        auto outputDataAcc = output->dataAccessor();
        auto inputSamplerFunc = getScalarSamplerFunc(input);
        auto inputDataPos = input.dataPosition();

        double h = min3(
            output->gridSpacing().x,
            output->gridSpacing().y,
            output->gridSpacing().z);

        output->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            if (boundarySdf.sample(inputDataPos(i, j, k)) > 0.0) {
                Vector3D pt = backTrace(
                    flow, dt, h, outputDataPos(i, j, k), boundarySdf);
                outputDataAcc(i, j, k) = inputSamplerFunc(pt);
            }
            });
    }

    void SemiLagrangian3::advect(
        const CollocatedVectorGrid3& input,
        const VectorField3& flow,
        double dt,
        CollocatedVectorGrid3* output,
        const ScalarField3& boundarySdf) {
        auto inputSamplerFunc = getVectorSamplerFunc(input);

        double h = min3(
            output->gridSpacing().x,
            output->gridSpacing().y,
            output->gridSpacing().z);

        auto outputDataPos = output->dataPosition();
        auto outputDataAcc = output->dataAccessor();
        auto inputDataPos = input.dataPosition();

        output->parallelForEachDataPointIndex([&](size_t i, size_t j, size_t k) {
            if (boundarySdf.sample(inputDataPos(i, j, k)) > 0.0) {
                Vector3D pt = backTrace(
                    flow, dt, h, outputDataPos(i, j, k), boundarySdf);
                outputDataAcc(i, j, k) = inputSamplerFunc(pt);
            }
            });
    }

    void SemiLagrangian3::advect(
        const FaceCenteredGrid3& input,
        const VectorField3& flow,
        double dt,
        FaceCenteredGrid3* output,
        const ScalarField3& boundarySdf) {
        auto inputSamplerFunc = getVectorSamplerFunc(input);

        double h = min3(
            output->gridSpacing().x,
            output->gridSpacing().y,
            output->gridSpacing().z);

        auto uTargetDataPos = output->uPosition();
        auto uTargetDataAcc = output->uAccessor();
        auto uSourceDataPos = input.uPosition();

        output->parallelForEachUIndex([&](size_t i, size_t j, size_t k) {
            if (boundarySdf.sample(uSourceDataPos(i, j, k)) > 0.0) {
                Vector3D pt = backTrace(
                    flow, dt, h, uTargetDataPos(i, j, k), boundarySdf);
                uTargetDataAcc(i, j, k) = inputSamplerFunc(pt).x;
            }
            });

        auto vTargetDataPos = output->vPosition();
        auto vTargetDataAcc = output->vAccessor();
        auto vSourceDataPos = input.vPosition();

        output->parallelForEachVIndex([&](size_t i, size_t j, size_t k) {
            if (boundarySdf.sample(vSourceDataPos(i, j, k)) > 0.0) {
                Vector3D pt = backTrace(
                    flow, dt, h, vTargetDataPos(i, j, k), boundarySdf);
                vTargetDataAcc(i, j, k) = inputSamplerFunc(pt).y;
            }
            });

        auto wTargetDataPos = output->wPosition();
        auto wTargetDataAcc = output->wAccessor();
        auto wSourceDataPos = input.wPosition();

        output->parallelForEachWIndex([&](size_t i, size_t j, size_t k) {
            if (boundarySdf.sample(wSourceDataPos(i, j, k)) > 0.0) {
                Vector3D pt = backTrace(
                    flow, dt, h, wTargetDataPos(i, j, k), boundarySdf);
                wTargetDataAcc(i, j, k) = inputSamplerFunc(pt).z;
            }
            });
    }

    Vector3D SemiLagrangian3::backTrace(
        const VectorField3& flow,
        double dt,
        double h,
        const Vector3D& startPt,
        const ScalarField3& boundarySdf) {

        double remainingT = dt;
        Vector3D pt0 = startPt;
        Vector3D pt1 = startPt;

        while (remainingT > kEpsilonD) {
            // Adaptive time-stepping
            Vector3D vel0 = flow.sample(pt0);
            double numSubSteps
                = std::max(std::ceil(vel0.length() * remainingT / h), 1.0);
            dt = remainingT / numSubSteps;

            // Mid-point rule
            Vector3D midPt = pt0 - 0.5 * dt * vel0;
            Vector3D midVel = flow.sample(midPt);
            pt1 = pt0 - dt * midVel;

            // Boundary handling
            double phi0 = boundarySdf.sample(pt0);
            double phi1 = boundarySdf.sample(pt1);

            if (phi0 * phi1 < 0.0) {
                double w = std::fabs(phi1) / (std::fabs(phi0) + std::fabs(phi1));
                pt1 = w * pt0 + (1.0 - w) * pt1;
                break;
            }

            remainingT -= dt;
            pt0 = pt1;
        }

        return pt1;
    }

    std::function<double(const Vector3D&)>
        SemiLagrangian3::getScalarSamplerFunc(const ScalarGrid3& input) const {
        return input.sampler();
    }

    std::function<Vector3D(const Vector3D&)>
        SemiLagrangian3::getVectorSamplerFunc(
            const CollocatedVectorGrid3& input) const {
        return input.sampler();
    }

    std::function<Vector3D(const Vector3D&)>
        SemiLagrangian3::getVectorSamplerFunc(const FaceCenteredGrid3& input) const {
        return input.sampler();
    }

}  // namespace jet

#endif  // INCLUDE_JET_SEMI_LAGRANGIAN3_H_