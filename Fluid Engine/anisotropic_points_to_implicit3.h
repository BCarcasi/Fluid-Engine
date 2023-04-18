#ifndef INCLUDE_JET_ANISOTROPIC_POINTS_TO_IMPLICIT3_H_
#define INCLUDE_JET_ANISOTROPIC_POINTS_TO_IMPLICIT3_H_

#include "pch.h"

#include "fmm_level_set_solver3.h"
#include "point_kdtree_searcher3.h"
#include "sph_kernels3.h"
#include "sph_system_data3.h"
#include "svd.h"

#include "points_to_implicit3.h"

namespace jet {

    //!
    //! \brief 3-D points-to-implicit converter using Anisotropic kernels.
    //!
    //! This class converts 3-D points to implicit surface using anisotropic kernels
    //! so that the kernels are oriented and stretched to reflect the point
    //! distribution more naturally (thus less bumps). The implementation is based
    //! on Yu and Turk's 2013 paper with some modifications.
    //!
    //! \see Yu, Jihun, and Greg Turk. "Reconstructing surfaces of particle-based
    //!      fluids using anisotropic kernels." ACM Transactions on Graphics (TOG)
    //!      32.1 (2013): 5.
    //!
    class AnisotropicPointsToImplicit3 final : public PointsToImplicit3 {
    public:
        //!
        //! \brief Constructs the converter with given parameters.
        //!
        //! \param kernelRadius Kernel radius for interpolations.
        //! \param cutOffDensity Iso-contour density value.
        //! \param positionSmoothingFactor Position smoothing factor.
        //! \param minNumNeighbors Minimum number of neighbors to enable anisotropic
        //!                        kernel.
        //!
        AnisotropicPointsToImplicit3(double kernelRadius = 1.0,
            double cutOffDensity = 0.5,
            double positionSmoothingFactor = 0.5,
            size_t minNumNeighbors = 25,
            bool isOutputSdf = true);

        //! Converts the given points to implicit surface scalar field.
        void convert(const ConstArrayAccessor1<Vector3D>& points,
            ScalarGrid3* output) const override;

    private:
        double _kernelRadius = 1.0;
        double _cutOffDensity = 0.5;
        double _positionSmoothingFactor = 0.0;
        size_t _minNumNeighbors = 25;
        bool _isOutputSdf = true;
    };

    //! Shared pointer for the AnisotropicPointsToImplicit3 type.
    typedef std::shared_ptr<AnisotropicPointsToImplicit3>
        AnisotropicPointsToImplicit3Ptr;

    inline double p(double distance) {
        const double distanceSquared = distance * distance;

        if (distanceSquared >= 1.0) {
            return 0.0;
        }
        else {
            const double x = 1.0 - distanceSquared;
            return x * x * x;
        }
    }

    inline double wij(double distance, double r) {
        if (distance < r) {
            return 1.0 - cubic(distance / r);
        }
        else {
            return 0.0;
        }
    }

    inline Matrix3x3D vvt(const Vector3D& v) {
        return Matrix3x3D(v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y,
            v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z);
    }

    inline double w(const Vector3D& r, const Matrix3x3D& g, double gDet) {
        static const double sigma = 315.0 / (64 * kPiD);
        return sigma * gDet * p((g * r).length());
    }

    //

    AnisotropicPointsToImplicit3::AnisotropicPointsToImplicit3(
        double kernelRadius, double cutOffDensity, double positionSmoothingFactor,
        size_t minNumNeighbors, bool isOutputSdf)
        : _kernelRadius(kernelRadius),
        _cutOffDensity(cutOffDensity),
        _positionSmoothingFactor(positionSmoothingFactor),
        _minNumNeighbors(minNumNeighbors),
        _isOutputSdf(isOutputSdf) {}

    void AnisotropicPointsToImplicit3::convert(
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

        JET_INFO << "Start converting points to implicit surface.";

        const double h = _kernelRadius;
        const double invH = 1 / h;
        const double r = 2.0 * h;

        // Mean estimator for cov. mat.
        const auto meanNeighborSearcher =
            PointKdTreeSearcher3::builder().makeShared();
        meanNeighborSearcher->build(points);

        JET_INFO << "Built neighbor searcher.";

        SphSystemData3 meanParticles;
        meanParticles.addParticles(points);
        meanParticles.setNeighborSearcher(meanNeighborSearcher);
        meanParticles.setKernelRadius(r);

        // Compute G and xMean
        std::vector<Matrix3x3D> gs(points.size());
        Array1<Vector3D> xMeans(points.size());

        parallelFor(kZeroSize, points.size(), [&](size_t i) {
            const auto& x = points[i];

            // Compute xMean
            Vector3D xMean;
            double wSum = 0.0;
            size_t numNeighbors = 0;
            const auto getXMean = [&](size_t, const Vector3D& xj) {
                const double wj = wij((x - xj).length(), r);
                wSum += wj;
                xMean += wj * xj;
                ++numNeighbors;
            };
            meanNeighborSearcher->forEachNearbyPoint(x, r, getXMean);


            JET_ASSERT(wSum > 0.0);
            xMean /= wSum;

            xMeans[i] = lerp(x, xMean, _positionSmoothingFactor);

            if (numNeighbors < _minNumNeighbors) {
                const auto g = Matrix3x3D::makeScaleMatrix(invH, invH, invH);
                gs[i] = g;
       
            }
            else {
                // Compute covariance matrix
                // We start with small scale matrix (h*h) in order to
                // prevent zero covariance matrix when points are all
                // perfectly lined up.
                auto cov = Matrix3x3D::makeScaleMatrix(h * h, h * h, h * h);
                wSum = 0.0;
                const auto getCov = [&](size_t, const Vector3D& xj) {
                    const double wj = wij((xMean - xj).length(), r);
                    wSum += wj;
                    cov += wj * vvt(xj - xMean);
                };
                meanNeighborSearcher->forEachNearbyPoint(x, r, getCov);

                cov /= wSum;

                // SVD
                Matrix3x3D u;
                Vector3D v;
                Matrix3x3D w;
                svd(cov, u, v, w);
           

                // Take off the sign
                v.x = std::fabs(v.x);
                v.y = std::fabs(v.y);
                v.z = std::fabs(v.z);

                // Constrain Sigma
                const double maxSingularVal = v.max();
                const double kr = 4.0;
                v.x = std::max(v.x, maxSingularVal / kr);
                v.y = std::max(v.y, maxSingularVal / kr);
                v.z = std::max(v.z, maxSingularVal / kr);

                const auto invSigma = Matrix3x3D::makeScaleMatrix(1.0 / v);

                // Compute G
                const double scale =
                    std::pow(v.x * v.y * v.z, 1.0 / 3.0);  // volume preservation
                const Matrix3x3D g = invH * scale * (w * invSigma * u.transposed());
                gs[i] = g;

            }
            });

        JET_INFO << "Computed G and means.";

        // SPH estimator
        meanParticles.setKernelRadius(h);
        meanParticles.updateDensities();
        const auto d = meanParticles.densities();
        const double m = meanParticles.mass();

        PointKdTreeSearcher3 meanNeighborSearcher2;
        meanNeighborSearcher2.build(xMeans);

        // Compute SDF
        auto temp = output->clone();
        temp->fill([&](const Vector3D& x) {
            double sum = 0.0;
            meanNeighborSearcher2.forEachNearbyPoint(
                x, r, [&](size_t i, const Vector3D& neighborPosition) {
                    sum += m / d[i] *
                        w(neighborPosition - x, gs[i], gs[i].determinant());
                });

            return _cutOffDensity - sum;
            });

        JET_INFO << "Computed SDF.";

        if (_isOutputSdf) {
            FmmLevelSetSolver3 solver;
            solver.reinitialize(*temp, kMaxD, output);

            JET_INFO << "Completed einitialization.";
        }
        else {
            temp->swap(output);
        }

        JET_INFO << "Done converting points to implicit surface.";
    }

}  // namespace jet

#endif  // INCLUDE_JET_ANISOTROPIC_POINTS_TO_IMPLICIT3_H_