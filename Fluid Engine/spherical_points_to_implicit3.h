#ifndef INCLUDE_JET_SPHERICAL_POINTS_TO_IMPLICIT3_H_
#define INCLUDE_JET_SPHERICAL_POINTS_TO_IMPLICIT3_H_

#include "points_to_implicit3.h"

#include "pch.h"

#include "fmm_level_set_solver3.h"
#include "particle_system_data3.h"


namespace jet {

    //!
    //! \brief 3-D points-to-implicit converter based on simple sphere model.
    //!
    class SphericalPointsToImplicit3 final : public PointsToImplicit3 {
    public:
        //! Constructs the converter with given sphere radius.
        SphericalPointsToImplicit3(double radius = 1.0, bool isOutputSdf = true);

        //! Converts the given points to implicit surface scalar field.
        void convert(const ConstArrayAccessor1<Vector3D>& points,
            ScalarGrid3* output) const override;

    private:
        double _radius = 1.0;
        bool _isOutputSdf = true;
    };

    //! Shared pointer type for SphericalPointsToImplicit3.
    typedef std::shared_ptr<SphericalPointsToImplicit3>
        SphericalPointsToImplicit3Ptr;

    SphericalPointsToImplicit3::SphericalPointsToImplicit3(double radius,
        bool isOutputSdf)
        : _radius(radius), _isOutputSdf(isOutputSdf) {}

    void SphericalPointsToImplicit3::convert(
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
        particles.buildNeighborSearcher(2.0 * _radius);

        const auto neighborSearcher = particles.neighborSearcher();

        auto temp = output->clone();
        temp->fill([&](const Vector3D& x) {
            double minDist = 2.0 * _radius;
            neighborSearcher->forEachNearbyPoint(
                x, 2.0 * _radius, [&](size_t, const Vector3D& xj) {
                    minDist = std::min(minDist, (x - xj).length());
                });

            return minDist - _radius;
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

#endif  // INCLUDE_JET_SPHERICAL_POINTS_TO_IMPLICIT3_H_