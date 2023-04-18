#ifndef INCLUDE_JET_ZHU_BRIDSON_POINTS_TO_IMPLICIT3_H_
#define INCLUDE_JET_ZHU_BRIDSON_POINTS_TO_IMPLICIT3_H_

#include "points_to_implicit3.h"


#include "pch.h"

#include "fmm_level_set_solver3.h"
#include "particle_system_data3.h"


namespace jet {

    //!
    //! \brief 3-D points-to-implicit converter based on Zhu and Bridson's method.
    //!
    //! \see Zhu, Yongning, and Robert Bridson. "Animating sand as a fluid."
    //!      ACM Transactions on Graphics (TOG). Vol. 24. No. 3. ACM, 2005.
    //!
    class ZhuBridsonPointsToImplicit3 final : public PointsToImplicit3 {
    public:
        //! Constructs the converter with given kernel radius and cut-off threshold.
        ZhuBridsonPointsToImplicit3(double kernelRadius = 1.0,
            double cutOffThreshold = 0.25,
            bool isOutputSdf = true);

        //! Converts the given points to implicit surface scalar field.
        void convert(const ConstArrayAccessor1<Vector3D>& points,
            ScalarGrid3* output) const override;

    private:
        double _kernelRadius = 1.0;
        double _cutOffThreshold = 0.25;
        bool _isOutputSdf = true;
    };

    //! Shared pointer type for ZhuBridsonPointsToImplicit3 class
    typedef std::shared_ptr<ZhuBridsonPointsToImplicit3> ZhuBridsonPointsToImplicit3Ptr;

    inline double k(double s) { return std::max(0.0, cubic(1.0 - s * s)); }

    ZhuBridsonPointsToImplicit3::ZhuBridsonPointsToImplicit3(double kernelRadius,
        double cutOffThreshold,
        bool isOutputSdf)
        : _kernelRadius(kernelRadius),
        _cutOffThreshold(cutOffThreshold),
        _isOutputSdf(isOutputSdf) {}

    void ZhuBridsonPointsToImplicit3::convert(
        const ConstArrayAccessor1<Vector3D>& points, ScalarGrid3* output) const {
        if (output == nullptr) {
            JET_WARN << "Null scalar grid output pointer provided.";
            return;
        }

        const auto res = output->resolution();
        if (res.x * res.y * res.z == 0) {
            JET_WARN << "Empty grid is provided.";
            return;
        }

        const auto bbox = output->boundingBox();
        if (bbox.isEmpty()) {
            JET_WARN << "Empty domain is provided.";
            return;
        }

        ParticleSystemData3 particles;
        particles.addParticles(points);
        particles.buildNeighborSearcher(_kernelRadius);

        const auto neighborSearcher = particles.neighborSearcher();
        const double isoContValue = _cutOffThreshold * _kernelRadius;

        auto temp = output->clone();
        temp->fill([&](const Vector3D& x) -> double {
            Vector3D xAvg;
            double wSum = 0.0;
            const auto func = [&](size_t, const Vector3D& xi) {
                const double wi = k((x - xi).length() / _kernelRadius);
                wSum += wi;
                xAvg += wi * xi;
            };
            neighborSearcher->forEachNearbyPoint(x, _kernelRadius, func);

            if (wSum > 0.0) {
                xAvg /= wSum;
                return (x - xAvg).length() - isoContValue;
            }
            else {
                return output->boundingBox().diagonalLength();
            }
            });

        if (_isOutputSdf) {
            FmmLevelSetSolver3 solver;
            solver.reinitialize(*temp, kMaxD, output);
        }
        else {
            temp->swap(output);
        }
    }

}  // namespace jet

#endif  // INCLUDE_JET_ZHU_BRIDSON_POINTS_TO_IMPLICIT3_H_