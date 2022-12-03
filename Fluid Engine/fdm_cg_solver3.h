#ifndef INCLUDE_JET_FDM_CG_SOLVER3_H_
#define INCLUDE_JET_FDM_CG_SOLVER3_H_

#include "fdm_linear_system_solver3.h"

#include "cg.h"
#include "constants.h"
#include "pch.h"

namespace jet {

    //! \brief 3-D finite difference-type linear system solver using conjugate
    //!        gradient.
    class FdmCgSolver3 final : public FdmLinearSystemSolver3 {
    public:
        //! Constructs the solver with given parameters.
        FdmCgSolver3(unsigned int maxNumberOfIterations, double tolerance);

        //! Solves the given linear system.
        bool solve(FdmLinearSystem3* system) override;

        //! Solves the given compressed linear system.
        bool solveCompressed(FdmCompressedLinearSystem3* system) override;

        //! Returns the max number of CG iterations.
        unsigned int maxNumberOfIterations() const;

        //! Returns the last number of CG iterations the solver made.
        unsigned int lastNumberOfIterations() const;

        //! Returns the max residual tolerance for the CG method.
        double tolerance() const;

        //! Returns the last residual after the CG iterations.
        double lastResidual() const;

    private:
        unsigned int _maxNumberOfIterations;
        unsigned int _lastNumberOfIterations;
        double _tolerance;
        double _lastResidual;

        // Uncompressed vectors
        FdmVector3 _r;
        FdmVector3 _d;
        FdmVector3 _q;
        FdmVector3 _s;

        // Compressed vectors
        VectorND _rComp;
        VectorND _dComp;
        VectorND _qComp;
        VectorND _sComp;

        void clearUncompressedVectors();
        void clearCompressedVectors();
    };

    //! Shared pointer type for the FdmCgSolver3.
    typedef std::shared_ptr<FdmCgSolver3> FdmCgSolver3Ptr;

    FdmCgSolver3::FdmCgSolver3(unsigned int maxNumberOfIterations, double tolerance)
        : _maxNumberOfIterations(maxNumberOfIterations),
        _lastNumberOfIterations(0),
        _tolerance(tolerance),
        _lastResidual(kMaxD) {}

    bool FdmCgSolver3::solve(FdmLinearSystem3* system) {
        FdmMatrix3& matrix = system->A;
        FdmVector3& solution = system->x;
        FdmVector3& rhs = system->b;

        JET_ASSERT(matrix.size() == rhs.size());
        JET_ASSERT(matrix.size() == solution.size());

        clearCompressedVectors();

        Size3 size = matrix.size();
        _r.resize(size);
        _d.resize(size);
        _q.resize(size);
        _s.resize(size);

        system->x.set(0.0);
        _r.set(0.0);
        _d.set(0.0);
        _q.set(0.0);
        _s.set(0.0);

        cg<FdmBlas3>(matrix, rhs, _maxNumberOfIterations, _tolerance, &solution,
            &_r, &_d, &_q, &_s, &_lastNumberOfIterations, &_lastResidual);

        return _lastResidual <= _tolerance ||
            _lastNumberOfIterations < _maxNumberOfIterations;
    }

    bool FdmCgSolver3::solveCompressed(FdmCompressedLinearSystem3* system) {
        MatrixCsrD& matrix = system->A;
        VectorND& solution = system->x;
        VectorND& rhs = system->b;

        clearUncompressedVectors();

        size_t size = solution.size();
        _rComp.resize(size);
        _dComp.resize(size);
        _qComp.resize(size);
        _sComp.resize(size);

        system->x.set(0.0);
        _rComp.set(0.0);
        _dComp.set(0.0);
        _qComp.set(0.0);
        _sComp.set(0.0);

        cg<FdmCompressedBlas3>(matrix, rhs, _maxNumberOfIterations, _tolerance,
            &solution, &_rComp, &_dComp, &_qComp, &_sComp,
            &_lastNumberOfIterations, &_lastResidual);

        return _lastResidual <= _tolerance ||
            _lastNumberOfIterations < _maxNumberOfIterations;
    }

    unsigned int FdmCgSolver3::maxNumberOfIterations() const {
        return _maxNumberOfIterations;
    }

    unsigned int FdmCgSolver3::lastNumberOfIterations() const {
        return _lastNumberOfIterations;
    }

    double FdmCgSolver3::tolerance() const { return _tolerance; }

    double FdmCgSolver3::lastResidual() const { return _lastResidual; }

    void FdmCgSolver3::clearUncompressedVectors() {
        _r.clear();
        _d.clear();
        _q.clear();
        _s.clear();
    }

    void FdmCgSolver3::clearCompressedVectors() {
        _rComp.clear();
        _dComp.clear();
        _qComp.clear();
        _sComp.clear();
    }

}  // namespace jet

#endif  // INCLUDE_JET_FDM_CG_SOLVER3_H_