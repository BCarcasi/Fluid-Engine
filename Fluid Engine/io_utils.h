#pragma once
#include "jet.h"
#include <pystring/pystring.h>

inline bool readPositions(const std::string& filename,
    jet::Array1<jet::Vector3D>& positions) {
    if (pystring::endswith(filename, ".xyz")) {
        std::ifstream xyzFile(filename.c_str());
        while (xyzFile) {
            jet::Vector3D xyz;

            xyzFile >> xyz.x >> xyz.y >> xyz.z;

            if (xyzFile.fail() || xyzFile.eof()) {
                break;
            }

            positions.append(xyz);
        }

        xyzFile.close();
        return true;
    }
    else {
        std::ifstream posFile(filename.c_str(), std::ifstream::binary);
        if (posFile) {
            std::vector<uint8_t> buffer(
                (std::istreambuf_iterator<char>(posFile)),
                (std::istreambuf_iterator<char>()));
            deserialize(buffer, &positions);
            posFile.close();
            return true;
        }
        else {
            return false;
        }
    }
}