#ifndef INCLUDE_JET_FDM_MG_SOLVER3_H_
#define INCLUDE_JET_FDM_MG_SOLVER3_H_

#include "fdm_linear_system_solver3.h"
#include "fdm_mg_linear_system3.h"
#include "mg.h"

#include "pch.h"

#include "fdm_gauss_seidel_solver3.h"

namespace jet {

    //! \brief 3-D finite difference-type linear system solver using Multigrid.
    class FdmMgSolver3 : public FdmLinearSystemSolver3 {
    public:
        FdmMgSolver3() = default;

        virtual ~FdmMgSolver3() = default;

        //! Constructs the solver with given parameters.
        FdmMgSolver3(size_t maxNumberOfLevels,
            unsigned int numberOfRestrictionIter = 5,
            unsigned int numberOfCorrectionIter = 5,
            unsigned int numberOfCoarsestIter = 20,
            unsigned int numberOfFinalIter = 20,
            double maxTolerance = 1e-9, double sorFactor = 1.5,
            bool useRedBlackOrdering = false);

        //! Returns the Multigrid parameters.
        const MgParameters<FdmBlas3>& params() const;

        //! Returns the SOR (Successive Over Relaxation) factor.
        double sorFactor() const;

        //! Returns true if red-black ordering is enabled.
        bool useRedBlackOrdering() const;

        //! No-op. Multigrid-type solvers do not solve FdmLinearSystem3.
        bool solve(FdmLinearSystem3* system) final;

        //! Solves Multigrid linear system.
        virtual bool solve(FdmMgLinearSystem3* system);

    private:
        MgParameters<FdmBlas3> _mgParams;
        double _sorFactor;
        bool _useRedBlackOrdering;
    };

    //! Shared pointer type for the FdmMgSolver3.
    typedef std::shared_ptr<FdmMgSolver3> FdmMgSolver3Ptr;

    FdmMgSolver3::FdmMgSolver3(size_t maxNumberOfLevels,
        unsigned int numberOfRestrictionIter,
        unsigned int numberOfCorrectionIter,
        unsigned int numberOfCoarsestIter,
        unsigned int numberOfFinalIter, double maxTolerance,
        double sorFactor, bool useRedBlackOrdering) {
        _mgParams.maxNumberOfLevels = maxNumberOfLevels;
        _mgParams.numberOfRestrictionIter = numberOfRestrictionIter;
        _mgParams.numberOfCorrectionIter = numberOfCorrectionIter;
        _mgParams.numberOfCoarsestIter = numberOfCoarsestIter;
        _mgParams.numberOfFinalIter = numberOfFinalIter;
        _mgParams.maxTolerance = maxTolerance;
        if (useRedBlackOrdering) {
            _mgParams.relaxFunc = [sorFactor](
                const FdmMatrix3& A, const FdmVector3& b,
                unsigned int numberOfIterations, double maxTolerance, FdmVector3* x,
                FdmVector3* buffer) {
                    UNUSED_VARIABLE(buffer);
                    UNUSED_VARIABLE(maxTolerance);

                    for (unsigned int iter = 0; iter < numberOfIterations; ++iter) {
                        FdmGaussSeidelSolver3::relaxRedBlack(A, b, sorFactor, x);
                    }
            };
        }
        else {
            _mgParams.relaxFunc = [sorFactor](
                const FdmMatrix3& A, const FdmVector3& b,
                unsigned int numberOfIterations, double maxTolerance, FdmVector3* x,
                FdmVector3* buffer) {
                    UNUSED_VARIABLE(buffer);
                    UNUSED_VARIABLE(maxTolerance);

                    for (unsigned int iter = 0; iter < numberOfIterations; ++iter) {
                        FdmGaussSeidelSolver3::relax(A, b, sorFactor, x);
                    }
            };
        }
        _mgParams.restrictFunc = FdmMgUtils3::restrict;
        _mgParams.correctFunc = FdmMgUtils3::correct;

        _sorFactor = sorFactor;
        _useRedBlackOrdering = useRedBlackOrdering;
    }

    const MgParameters<FdmBlas3>& FdmMgSolver3::params() const { return _mgParams; }

    double FdmMgSolver3::sorFactor() const { return _sorFactor; }

    bool FdmMgSolver3::useRedBlackOrdering() const { return _useRedBlackOrdering; }

    bool FdmMgSolver3::solve(FdmLinearSystem3* system) {
        UNUSED_VARIABLE(system);
        return false;
    }

    bool FdmMgSolver3::solve(FdmMgLinearSystem3* system) {
        FdmMgVector3 buffer = system->x;
        auto result =
            mgVCycle(system->A, _mgParams, &system->x, &system->b, &buffer);
        return result.lastResidualNorm < _mgParams.maxTolerance;
    }

}  // namespace jet

#endif  // INCLUDE_JET_FDM_MG_SOLVER3_H_