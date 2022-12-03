#ifndef INCLUDE_JET_FDM_MG_LINEAR_SYSTEM3_H_
#define INCLUDE_JET_FDM_MG_LINEAR_SYSTEM3_H_

#include "face_centered_grid3.h"
#include "fdm_linear_system3.h"
#include "mg.h"

namespace jet {

    //! Multigrid-style 3-D FDM matrix.
    typedef MgMatrix<FdmBlas3> FdmMgMatrix3;

    //! Multigrid-style 3-D FDM vector.
    typedef MgVector<FdmBlas3> FdmMgVector3;

    //! Multigrid-syle 3-D linear system.
    struct FdmMgLinearSystem3 {
        //! The system matrix.
        FdmMgMatrix3 A;

        //! The solution vector.
        FdmMgVector3 x;

        //! The RHS vector.
        FdmMgVector3 b;

        //! Clears the linear system.
        void clear();

        //! Returns the number of multigrid levels.
        size_t numberOfLevels() const;

        //! Resizes the system with the coarsest resolution and number of levels.
        void resizeWithCoarsest(const Size3& coarsestResolution,
            size_t numberOfLevels);

        //!
        //! \brief Resizes the system with the finest resolution and max number of
        //! levels.
        //!
        //! This function resizes the system with multiple levels until the
        //! resolution is divisible with 2^(level-1).
        //!
        //! \param finestResolution - The finest grid resolution.
        //! \param maxNumberOfLevels - Maximum number of multigrid levels.
        //!
        void resizeWithFinest(const Size3& finestResolution,
            size_t maxNumberOfLevels);
    };

    //! Multigrid utilities for 2-D FDM system.
    class FdmMgUtils3 {
    public:
        //! Restricts given finer grid to the coarser grid.
        static void restrict(const FdmVector3& finer, FdmVector3* coarser);

        //! Corrects given coarser grid to the finer grid.
        static void correct(const FdmVector3& coarser, FdmVector3* finer);

        //! Resizes the array with the coarsest resolution and number of levels.
        template <typename T>
        static void resizeArrayWithCoarsest(const Size3& coarsestResolution,
            size_t numberOfLevels,
            std::vector<Array3<T>>* levels);

        //!
        //! \brief Resizes the array with the finest resolution and max number of
        //! levels.
        //!
        //! This function resizes the system with multiple levels until the
        //! resolution is divisible with 2^(level-1).
        //!
        //! \param finestResolution - The finest grid resolution.
        //! \param maxNumberOfLevels - Maximum number of multigrid levels.
        //!
        template <typename T>
        static void resizeArrayWithFinest(const Size3& finestResolution,
            size_t maxNumberOfLevels,
            std::vector<Array3<T>>* levels);
    };

    template <typename T>
    void FdmMgUtils3::resizeArrayWithCoarsest(const Size3& coarsestResolution,
        size_t numberOfLevels,
        std::vector<Array3<T>>* levels) {
        numberOfLevels = std::max(numberOfLevels, kOneSize);

        levels->resize(numberOfLevels);

        // Level 0 is the finest level, thus takes coarsestResolution ^
        // numberOfLevels.
        // Level numberOfLevels - 1 is the coarsest, taking coarsestResolution.
        Size3 res = coarsestResolution;
        for (size_t level = 0; level < numberOfLevels; ++level) {
            (*levels)[numberOfLevels - level - 1].resize(res);
            res.x = res.x << 1;
            res.y = res.y << 1;
            res.z = res.z << 1;
        }
    }

    template <typename T>
    void FdmMgUtils3::resizeArrayWithFinest(const Size3& finestResolution,
        size_t maxNumberOfLevels,
        std::vector<Array3<T>>* levels) {
        Size3 res = finestResolution;
        size_t i = 1;
        for (; i < maxNumberOfLevels; ++i) {
            if (res.x % 2 == 0 && res.y % 2 == 0 && res.z % 2 == 0) {
                res.x = res.x >> 1;
                res.y = res.y >> 1;
                res.z = res.z >> 1;
            }
            else {
                break;
            }
        }
        resizeArrayWithCoarsest(res, i, levels);
    }

    void FdmMgLinearSystem3::clear() {
        A.levels.clear();
        x.levels.clear();
        b.levels.clear();
    }

    size_t FdmMgLinearSystem3::numberOfLevels() const { return A.levels.size(); }

    void FdmMgLinearSystem3::resizeWithCoarsest(const Size3& coarsestResolution,
        size_t numberOfLevels) {
        FdmMgUtils3::resizeArrayWithCoarsest(coarsestResolution, numberOfLevels,
            &A.levels);
        FdmMgUtils3::resizeArrayWithCoarsest(coarsestResolution, numberOfLevels,
            &x.levels);
        FdmMgUtils3::resizeArrayWithCoarsest(coarsestResolution, numberOfLevels,
            &b.levels);
    }

    void FdmMgLinearSystem3::resizeWithFinest(const Size3& finestResolution,
        size_t maxNumberOfLevels) {
        FdmMgUtils3::resizeArrayWithFinest(finestResolution, maxNumberOfLevels,
            &A.levels);
        FdmMgUtils3::resizeArrayWithFinest(finestResolution, maxNumberOfLevels,
            &x.levels);
        FdmMgUtils3::resizeArrayWithFinest(finestResolution, maxNumberOfLevels,
            &b.levels);
    }

    void FdmMgUtils3::restrict(const FdmVector3& finer, FdmVector3* coarser) {
        JET_ASSERT(finer.size().x == 2 * coarser->size().x);
        JET_ASSERT(finer.size().y == 2 * coarser->size().y);
        JET_ASSERT(finer.size().z == 2 * coarser->size().z);

        // --*--|--*--|--*--|--*--
        //  1/8   3/8   3/8   1/8
        //           to
        // -----|-----*-----|-----
        static const std::array<double, 4> kernel = { {0.125, 0.375, 0.375, 0.125} };

        const Size3 n = coarser->size();
        parallelRangeFor(
            kZeroSize, n.x, kZeroSize, n.y, kZeroSize, n.z,
            [&](size_t iBegin, size_t iEnd, size_t jBegin, size_t jEnd,
                size_t kBegin, size_t kEnd) {
                    std::array<size_t, 4> kIndices;

                    for (size_t k = kBegin; k < kEnd; ++k) {
                        kIndices[0] = (k > 0) ? 2 * k - 1 : 2 * k;
                        kIndices[1] = 2 * k;
                        kIndices[2] = 2 * k + 1;
                        kIndices[3] = (k + 1 < n.z) ? 2 * k + 2 : 2 * k + 1;

                        std::array<size_t, 4> jIndices;

                        for (size_t j = jBegin; j < jEnd; ++j) {
                            jIndices[0] = (j > 0) ? 2 * j - 1 : 2 * j;
                            jIndices[1] = 2 * j;
                            jIndices[2] = 2 * j + 1;
                            jIndices[3] = (j + 1 < n.y) ? 2 * j + 2 : 2 * j + 1;

                            std::array<size_t, 4> iIndices;
                            for (size_t i = iBegin; i < iEnd; ++i) {
                                iIndices[0] = (i > 0) ? 2 * i - 1 : 2 * i;
                                iIndices[1] = 2 * i;
                                iIndices[2] = 2 * i + 1;
                                iIndices[3] = (i + 1 < n.x) ? 2 * i + 2 : 2 * i + 1;

                                double sum = 0.0;
                                for (size_t z = 0; z < 4; ++z) {
                                    for (size_t y = 0; y < 4; ++y) {
                                        for (size_t x = 0; x < 4; ++x) {
                                            double w =
                                                kernel[x] * kernel[y] * kernel[z];
                                            sum += w * finer(iIndices[x], jIndices[y],
                                                kIndices[z]);
                                        }
                                    }
                                }
                                (*coarser)(i, j, k) = sum;
                            }
                        }
                    }
            });
    }

    void FdmMgUtils3::correct(const FdmVector3& coarser, FdmVector3* finer) {
        JET_ASSERT(finer->size().x == 2 * coarser.size().x);
        JET_ASSERT(finer->size().y == 2 * coarser.size().y);
        JET_ASSERT(finer->size().z == 2 * coarser.size().z);

        // -----|-----*-----|-----
        //           to
        //  1/4   3/4   3/4   1/4
        // --*--|--*--|--*--|--*--
        const Size3 n = finer->size();
        parallelRangeFor(
            kZeroSize, n.x, kZeroSize, n.y, kZeroSize, n.z,
            [&](size_t iBegin, size_t iEnd, size_t jBegin, size_t jEnd,
                size_t kBegin, size_t kEnd) {
                    for (size_t k = kBegin; k < kEnd; ++k) {
                        for (size_t j = jBegin; j < jEnd; ++j) {
                            for (size_t i = iBegin; i < iEnd; ++i) {
                                std::array<size_t, 2> iIndices;
                                std::array<size_t, 2> jIndices;
                                std::array<size_t, 2> kIndices;
                                std::array<double, 2> iWeights;
                                std::array<double, 2> jWeights;
                                std::array<double, 2> kWeights;

                                const size_t ci = i / 2;
                                const size_t cj = j / 2;
                                const size_t ck = k / 2;

                                if (i % 2 == 0) {
                                    iIndices[0] = (i > 1) ? ci - 1 : ci;
                                    iIndices[1] = ci;
                                    iWeights[0] = 0.25;
                                    iWeights[1] = 0.75;
                                }
                                else {
                                    iIndices[0] = ci;
                                    iIndices[1] = (i + 1 < n.x) ? ci + 1 : ci;
                                    iWeights[0] = 0.75;
                                    iWeights[1] = 0.25;
                                }

                                if (j % 2 == 0) {
                                    jIndices[0] = (j > 1) ? cj - 1 : cj;
                                    jIndices[1] = cj;
                                    jWeights[0] = 0.25;
                                    jWeights[1] = 0.75;
                                }
                                else {
                                    jIndices[0] = cj;
                                    jIndices[1] = (j + 1 < n.y) ? cj + 1 : cj;
                                    jWeights[0] = 0.75;
                                    jWeights[1] = 0.25;
                                }

                                if (k % 2 == 0) {
                                    kIndices[0] = (k > 1) ? ck - 1 : ck;
                                    kIndices[1] = ck;
                                    kWeights[0] = 0.25;
                                    kWeights[1] = 0.75;
                                }
                                else {
                                    kIndices[0] = ck;
                                    kIndices[1] = (k + 1 < n.y) ? ck + 1 : ck;
                                    kWeights[0] = 0.75;
                                    kWeights[1] = 0.25;
                                }

                                for (size_t z = 0; z < 2; ++z) {
                                    for (size_t y = 0; y < 2; ++y) {
                                        for (size_t x = 0; x < 2; ++x) {
                                            double w = iWeights[x] * jWeights[y] *
                                                kWeights[z] *
                                                coarser(iIndices[x], jIndices[y],
                                                    kIndices[z]);
                                            (*finer)(i, j, k) += w;
                                        }
                                    }
                                }
                            }
                        }
                    }
            });
    }

}  // namespace jet


#endif  // INCLUDE_JET_FDM_MG_LINEAR_SYSTEM3_H_