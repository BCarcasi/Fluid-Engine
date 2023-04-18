#ifndef INCLUDE_JET_GRID_POINT_GENERATOR3_H_
#define INCLUDE_JET_GRID_POINT_GENERATOR3_H_

#include "pch.h"
#include "point_generator3.h"

namespace jet {

    //!
    //! \brief 3-D regular-grid point generator.
    //!
    class GridPointGenerator3 final : public PointGenerator3 {
    public:
        //!
        //! \brief Invokes \p callback function for each regular grid points inside
        //! \p boundingBox.
        //!
        //! This function iterates every regular grid points inside \p boundingBox
        //! where \p spacing is the size of the unit cell of regular grid structure.
        //!
        void forEachPoint(
            const BoundingBox3D& boundingBox,
            double spacing,
            const std::function<bool(const Vector3D&)>& callback) const;
    };

    //! Shared pointer type for the GridPointGenerator3.
    typedef std::shared_ptr<GridPointGenerator3> GridPointGenerator3Ptr;

    void GridPointGenerator3::forEachPoint(
        const BoundingBox3D& boundingBox,
        double spacing,
        const std::function<bool(const Vector3D&)>& callback) const {
        Vector3D position;
        double boxWidth = boundingBox.width();
        double boxHeight = boundingBox.height();
        double boxDepth = boundingBox.depth();

        bool shouldQuit = false;
        for (int k = 0; k * spacing <= boxDepth && !shouldQuit; ++k) {
            position.z = k * spacing + boundingBox.lowerCorner.z;

            for (int j = 0; j * spacing <= boxHeight && !shouldQuit; ++j) {
                position.y = j * spacing + boundingBox.lowerCorner.y;

                for (int i = 0; i * spacing <= boxWidth && !shouldQuit; ++i) {
                    position.x = i * spacing + boundingBox.lowerCorner.x;
                    if (!callback(position)) {
                        shouldQuit = true;
                        break;
                    }
                }
            }
        }
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_POINT_GENERATOR3_H_