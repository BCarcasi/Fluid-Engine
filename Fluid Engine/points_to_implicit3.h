#ifndef INCLUDE_JET_POINTS_TO_IMPLICIT3_H_
#define INCLUDE_JET_POINTS_TO_IMPLICIT3_H_

#include "array_accessor1.h"
#include "scalar_grid3.h"
#include "vector3.h"

#include <memory>

namespace jet {

    //! Abstract base class for 3-D points-to-implicit converters.
    class PointsToImplicit3 {
    public:
        //! Default constructor.
        PointsToImplicit3();

        //! Default destructor.
        virtual ~PointsToImplicit3();

        //! Converts the given points to implicit surface scalar field.
        virtual void convert(const ConstArrayAccessor1<Vector3D>& points,
            ScalarGrid3* output) const = 0;
    };

    //! Shared pointer for the PointsToImplicit3 type.
    typedef std::shared_ptr<PointsToImplicit3> PointsToImplicit3Ptr;

    PointsToImplicit3::PointsToImplicit3() {}

    PointsToImplicit3::~PointsToImplicit3() {}

}  // namespace jet

#endif  // INCLUDE_JET_POINTS_TO_IMPLICIT3_H_