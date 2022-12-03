#ifndef INCLUDE_JET_FDM_LINEAR_SYSTEM3_H_
#define INCLUDE_JET_FDM_LINEAR_SYSTEM3_H_

#include "array1.h"
#include "array3.h"
#include "matrix_csr.h"
#include "vector_n.h"

#include "pch.h"
#include "math_utils.h"
#include "parallel.h"

namespace jet {

    //! The row of FdmMatrix3 where row corresponds to (i, j, k) grid point.
    struct FdmMatrixRow3 {
        //! Diagonal component of the matrix (row, row).
        double center = 0.0;

        //! Off-diagonal element where colum refers to (i+1, j, k) grid point.
        double right = 0.0;

        //! Off-diagonal element where column refers to (i, j+1, k) grid point.
        double up = 0.0;

        //! OFf-diagonal element where column refers to (i, j, k+1) grid point.
        double front = 0.0;
    };

    //! Vector type for 3-D finite differencing.
    typedef Array3<double> FdmVector3;

    //! Matrix type for 3-D finite differencing.
    typedef Array3<FdmMatrixRow3> FdmMatrix3;

    //! Linear system (Ax=b) for 3-D finite differencing.
    struct FdmLinearSystem3 {
        //! System matrix.
        FdmMatrix3 A;

        //! Solution vector.
        FdmVector3 x;

        //! RHS vector.
        FdmVector3 b;

        //! Clears all the data.
        void clear();

        //! Resizes the arrays with given grid size.
        void resize(const Size3& size);
    };

    //! Compressed linear system (Ax=b) for 3-D finite differencing.
    struct FdmCompressedLinearSystem3 {
        //! System matrix.
        MatrixCsrD A;

        //! Solution vector.
        VectorND x;

        //! RHS vector.
        VectorND b;

        //! Clears all the data.
        void clear();
    };

    //! BLAS operator wrapper for 3-D finite differencing.
    struct FdmBlas3 {
        typedef double ScalarType;
        typedef FdmVector3 VectorType;
        typedef FdmMatrix3 MatrixType;

        //! Sets entire element of given vector \p result with scalar \p s.
        static void set(ScalarType s, VectorType* result);

        //! Copies entire element of given vector \p result with other vector \p v.
        static void set(const VectorType& v, VectorType* result);

        //! Sets entire element of given matrix \p result with scalar \p s.
        static void set(ScalarType s, MatrixType* result);

        //! Copies entire element of given matrix \p result with other matrix \p v.
        static void set(const MatrixType& m, MatrixType* result);

        //! Performs dot product with vector \p a and \p b.
        static double dot(const VectorType& a, const VectorType& b);

        //! Performs ax + y operation where \p a is a matrix and \p x and \p y are
        //! vectors.
        static void axpy(double a, const VectorType& x, const VectorType& y,
            VectorType* result);

        //! Performs matrix-vector multiplication.
        static void mvm(const MatrixType& m, const VectorType& v,
            VectorType* result);

        //! Computes residual vector (b - ax).
        static void residual(const MatrixType& a, const VectorType& x,
            const VectorType& b, VectorType* result);

        //! Returns L2-norm of the given vector \p v.
        static ScalarType l2Norm(const VectorType& v);

        //! Returns Linf-norm of the given vector \p v.
        static ScalarType lInfNorm(const VectorType& v);
    };

    //! BLAS operator wrapper for compressed 3-D finite differencing.
    struct FdmCompressedBlas3 {
        typedef double ScalarType;
        typedef VectorND VectorType;
        typedef MatrixCsrD MatrixType;

        //! Sets entire element of given vector \p result with scalar \p s.
        static void set(ScalarType s, VectorType* result);

        //! Copies entire element of given vector \p result with other vector \p v.
        static void set(const VectorType& v, VectorType* result);

        //! Sets entire element of given matrix \p result with scalar \p s.
        static void set(ScalarType s, MatrixType* result);

        //! Copies entire element of given matrix \p result with other matrix \p v.
        static void set(const MatrixType& m, MatrixType* result);

        //! Performs dot product with vector \p a and \p b.
        static double dot(const VectorType& a, const VectorType& b);

        //! Performs ax + y operation where \p a is a matrix and \p x and \p y are
        //! vectors.
        static void axpy(double a, const VectorType& x, const VectorType& y,
            VectorType* result);

        //! Performs matrix-vector multiplication.
        static void mvm(const MatrixType& m, const VectorType& v,
            VectorType* result);

        //! Computes residual vector (b - ax).
        static void residual(const MatrixType& a, const VectorType& x,
            const VectorType& b, VectorType* result);

        //! Returns L2-norm of the given vector \p v.
        static ScalarType l2Norm(const VectorType& v);

        //! Returns Linf-norm of the given vector \p v.
        static ScalarType lInfNorm(const VectorType& v);
    };

    void FdmLinearSystem3::clear() {
        A.clear();
        x.clear();
        b.clear();
    }

    void FdmLinearSystem3::resize(const Size3& size) {
        A.resize(size);
        x.resize(size);
        b.resize(size);
    }

    //

    void FdmCompressedLinearSystem3::clear() {
        A.clear();
        x.clear();
        b.clear();
    }

    //

    void FdmBlas3::set(double s, FdmVector3* result) { result->set(s); }

    void FdmBlas3::set(const FdmVector3& v, FdmVector3* result) { result->set(v); }

    void FdmBlas3::set(double s, FdmMatrix3* result) {
        FdmMatrixRow3 row;
        row.center = row.right = row.up = row.front = s;
        result->set(row);
    }

    void FdmBlas3::set(const FdmMatrix3& m, FdmMatrix3* result) { result->set(m); }

    double FdmBlas3::dot(const FdmVector3& a, const FdmVector3& b) {
        Size3 size = a.size();

        JET_THROW_INVALID_ARG_IF(size != b.size());

        double result = 0.0;

        for (size_t k = 0; k < size.z; ++k) {
            for (size_t j = 0; j < size.y; ++j) {
                for (size_t i = 0; i < size.x; ++i) {
                    result += a(i, j, k) * b(i, j, k);
                }
            }
        }

        return result;
    }

    void FdmBlas3::axpy(double a, const FdmVector3& x, const FdmVector3& y,
        FdmVector3* result) {
        Size3 size = x.size();

        JET_THROW_INVALID_ARG_IF(size != y.size());
        JET_THROW_INVALID_ARG_IF(size != result->size());

        x.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            (*result)(i, j, k) = a * x(i, j, k) + y(i, j, k);
            });
    }

    void FdmBlas3::mvm(const FdmMatrix3& m, const FdmVector3& v,
        FdmVector3* result) {
        Size3 size = m.size();

        JET_THROW_INVALID_ARG_IF(size != v.size());
        JET_THROW_INVALID_ARG_IF(size != result->size());

        m.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            (*result)(i, j, k) =
                m(i, j, k).center * v(i, j, k) +
                ((i > 0) ? m(i - 1, j, k).right * v(i - 1, j, k) : 0.0) +
                ((i + 1 < size.x) ? m(i, j, k).right * v(i + 1, j, k) : 0.0) +
                ((j > 0) ? m(i, j - 1, k).up * v(i, j - 1, k) : 0.0) +
                ((j + 1 < size.y) ? m(i, j, k).up * v(i, j + 1, k) : 0.0) +
                ((k > 0) ? m(i, j, k - 1).front * v(i, j, k - 1) : 0.0) +
                ((k + 1 < size.z) ? m(i, j, k).front * v(i, j, k + 1) : 0.0);
            });
    }

    void FdmBlas3::residual(const FdmMatrix3& a, const FdmVector3& x,
        const FdmVector3& b, FdmVector3* result) {
        Size3 size = a.size();

        JET_THROW_INVALID_ARG_IF(size != x.size());
        JET_THROW_INVALID_ARG_IF(size != b.size());
        JET_THROW_INVALID_ARG_IF(size != result->size());

        a.parallelForEachIndex([&](size_t i, size_t j, size_t k) {
            (*result)(i, j, k) =
                b(i, j, k) - a(i, j, k).center * x(i, j, k) -
                ((i > 0) ? a(i - 1, j, k).right * x(i - 1, j, k) : 0.0) -
                ((i + 1 < size.x) ? a(i, j, k).right * x(i + 1, j, k) : 0.0) -
                ((j > 0) ? a(i, j - 1, k).up * x(i, j - 1, k) : 0.0) -
                ((j + 1 < size.y) ? a(i, j, k).up * x(i, j + 1, k) : 0.0) -
                ((k > 0) ? a(i, j, k - 1).front * x(i, j, k - 1) : 0.0) -
                ((k + 1 < size.z) ? a(i, j, k).front * x(i, j, k + 1) : 0.0);
            });
    }

    double FdmBlas3::l2Norm(const FdmVector3& v) { return std::sqrt(dot(v, v)); }

    double FdmBlas3::lInfNorm(const FdmVector3& v) {
        Size3 size = v.size();

        double result = 0.0;

        for (size_t k = 0; k < size.z; ++k) {
            for (size_t j = 0; j < size.y; ++j) {
                for (size_t i = 0; i < size.x; ++i) {
                    result = absmax(result, v(i, j, k));
                }
            }
        }

        return std::fabs(result);
    }

    //

    void FdmCompressedBlas3::set(double s, VectorND* result) { result->set(s); }

    void FdmCompressedBlas3::set(const VectorND& v, VectorND* result) {
        result->set(v);
    }

    void FdmCompressedBlas3::set(double s, MatrixCsrD* result) { result->set(s); }

    void FdmCompressedBlas3::set(const MatrixCsrD& m, MatrixCsrD* result) {
        result->set(m);
    }

    double FdmCompressedBlas3::dot(const VectorND& a, const VectorND& b) {
        return a.dot(b);
    }

    void FdmCompressedBlas3::axpy(double a, const VectorND& x, const VectorND& y,
        VectorND* result) {
        *result = a * x + y;
    }

    void FdmCompressedBlas3::mvm(const MatrixCsrD& m, const VectorND& v,
        VectorND* result) {
        const auto rp = m.rowPointersBegin();
        const auto ci = m.columnIndicesBegin();
        const auto nnz = m.nonZeroBegin();

        v.parallelForEachIndex([&](size_t i) {
            const size_t rowBegin = rp[i];
            const size_t rowEnd = rp[i + 1];

            double sum = 0.0;

            for (size_t jj = rowBegin; jj < rowEnd; ++jj) {
                size_t j = ci[jj];
                sum += nnz[jj] * v[j];
            }

            (*result)[i] = sum;
            });
    }

    void FdmCompressedBlas3::residual(const MatrixCsrD& a, const VectorND& x,
        const VectorND& b, VectorND* result) {
        const auto rp = a.rowPointersBegin();
        const auto ci = a.columnIndicesBegin();
        const auto nnz = a.nonZeroBegin();

        x.parallelForEachIndex([&](size_t i) {
            const size_t rowBegin = rp[i];
            const size_t rowEnd = rp[i + 1];

            double sum = 0.0;

            for (size_t jj = rowBegin; jj < rowEnd; ++jj) {
                size_t j = ci[jj];
                sum += nnz[jj] * x[j];
            }

            (*result)[i] = b[i] - sum;
            });
    }

    double FdmCompressedBlas3::l2Norm(const VectorND& v) {
        return std::sqrt(v.dot(v));
    }

    double FdmCompressedBlas3::lInfNorm(const VectorND& v) {
        return std::fabs(v.absmax());
    }

}  // namespace jet

#endif  // INCLUDE_JET_FDM_LINEAR_SYSTEM3_H_