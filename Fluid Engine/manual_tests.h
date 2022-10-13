#ifndef SRC_TESTS_MANUAL_TESTS_MANUAL_TESTS_H_
#define SRC_TESTS_MANUAL_TESTS_MANUAL_TESTS_H_

#include "array_accessor1.h"
#include "array_accessor2.h"
#include "array_accessor3.h"
#include "triangle_mesh3.h"

#include <xtensor/xarray.hpp>
#include <pystring/pystring.h>
#include <xtensor/xio.hpp>
#include <xtensor-io/xhighfive.hpp>

#include <fstream>
#include <memory>
#include <string>
#include <vector>
#ifdef JET_WINDOWS
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace jet {
		template <typename T>
		void saveData(
			const ConstArrayAccessor1<T>& data,
			const std::string& name) {
			xt::dump_hdf5(name, "", data.data());
		}


}

#endif  // SRC_TESTS_MANUAL_TESTS_MANUAL_TESTS_H_