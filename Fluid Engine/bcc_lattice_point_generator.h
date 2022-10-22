#ifndef INCLUDE_JET_BCC_LATTICE_POINT_GENERATOR_H_
#define INCLUDE_JET_BCC_LATTICE_POINT_GENERATOR_H_

#include "point_generator3.h"

#include "pch.h"

namespace jet {

    //!
    //! \brief Body-centered lattice points generator.
    //!
    //! \see http://en.wikipedia.org/wiki/Cubic_crystal_system
    //!      http://mathworld.wolfram.com/CubicClosePacking.html
    //!
    class BccLatticePointGenerator final : public PointGenerator3 {
    public:
        //!
        //! \brief Invokes \p callback function for each BCC-lattice points inside
        //! \p boundingBox.
        //!
        //! This function iterates every BCC-lattice points inside \p boundingBox
        //! where \p spacing is the size of the unit cell of BCC structure.
        //!
        void forEachPoint(
            const BoundingBox3D& boundingBox,
            double spacing,
            const std::function<bool(const Vector3D&)>& callback) const override;
    };

    //! Shared pointer type for the BccLatticePointGenerator.
    typedef std::shared_ptr<BccLatticePointGenerator> BccLatticePointGeneratorPtr;

    void BccLatticePointGenerator::forEachPoint(
        const BoundingBox3D& boundingBox,
        double spacing,
        const std::function<bool(const Vector3D&)>& callback) const {
        double halfSpacing = spacing / 2.0;
        double boxWidth = boundingBox.width();
        double boxHeight = boundingBox.height();
        double boxDepth = boundingBox.depth();

        Vector3D position;
        bool hasOffset = false;
        bool shouldQuit = false;
        for (int k = 0; k * halfSpacing <= boxDepth && !shouldQuit; ++k) {
            position.z = k * halfSpacing + boundingBox.lowerCorner.z;

            double offset = (hasOffset) ? halfSpacing : 0.0;

            for (int j = 0; j * spacing + offset <= boxHeight && !shouldQuit; ++j) {
                position.y = j * spacing + offset + boundingBox.lowerCorner.y;

                for (int i = 0; i * spacing + offset <= boxWidth; ++i) {
                    position.x = i * spacing + offset + boundingBox.lowerCorner.x;
                    if (!callback(position)) {
                        shouldQuit = true;
                        break;
                    }
                }
            }

            hasOffset = !hasOffset;
        }
    }


}  // namespace jet

#endif  // INCLUDE_JET_BCC_LATTICE_POINT_GENERATOR_H_