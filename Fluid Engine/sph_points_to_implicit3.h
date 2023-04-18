#ifndef INCLUDE_JET_SPH_POINTS_TO_IMPLICIT3_H_
#define INCLUDE_JET_SPH_POINTS_TO_IMPLICIT3_H_

#include "points_to_implicit3.h"

#include "pch.h"

#include "fmm_level_set_solver3.h"
#include "sph_points_to_implicit3.h"
#include "sph_system_data3.h"

namespace jet {

    //!
    //! \brief 3-D points-to-implicit converter based on standard SPH kernel.
    //!
    //! \see Müller, Matthias, David Charypar, and Markus Gross.
    //!      "Particle-based fluid simulation for interactive applications."
    //!      Proceedings of the 2003 ACM SIGGRAPH/Eurographics symposium on Computer
    //!      animation. Eurographics Association, 2003.
    //!
    class SphPointsToImplicit3 final : public PointsToImplicit3 {
    public:
        //! Constructs the converter with given kernel radius and cut-off density.
        SphPointsToImplicit3(double kernelRadius = 1.0, double cutOffDensity = 0.5,
            bool isOutputSdf = true);

        //! Converts the given points to implicit surface scalar field.
        void convert(const ConstArrayAccessor1<Vector3D>& points,
            ScalarGrid3* output) const override;

    private:
        double _kernelRadius = 1.0;
        double _cutOffDensity = 0.5;
        bool _isOutputSdf = true;
    };

    //! Shared pointer type for SphPointsToImplicit3 class.
    typedef std::shared_ptr<SphPointsToImplicit3> SphPointsToImplicit3Ptr;

    SphPointsToImplicit3::SphPointsToImplicit3(double kernelRadius,
        double cutOffDensity,
        bool isOutputSdf)
        : _kernelRadius(kernelRadius),
        _cutOffDensity(cutOffDensity),
        _isOutputSdf(isOutputSdf) {}

    void SphPointsToImplicit3::convert(const ConstArrayAccessor1<Vector3D>& points,
        ScalarGrid3* output) const {
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

        SphSystemData3 sphParticles;
        sphParticles.addParticles(points);
        sphParticles.setKernelRadius(_kernelRadius);
        sphParticles.buildNeighborSearcher();
        sphParticles.updateDensities();

        Array1<double> constData(sphParticles.numberOfParticles(), 1.0);
        auto temp = output->clone();
        temp->fill([&](const Vector3D& x) {
            double d = sphParticles.interpolate(x, constData);
            return _cutOffDensity - d;
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

#endif  // INCLUDE_JET_SPH_POINTS_TO_IMPLICIT3_H_