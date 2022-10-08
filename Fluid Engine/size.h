#ifndef INCLUDE_JET_SIZE_H_
#define INCLUDE_JET_SIZE_H_

#include "point.h"

namespace jet {

	// N-D size type.
	template <size_t N> using Size = Point<size_t, N>;

}  // namespace jet

// #include "detail/size-inl.h"

#endif  // INCLUDE_JET_SIZE_H_