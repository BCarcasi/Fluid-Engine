#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <gl/GL.h>

namespace jet {
    struct point3 {
        GLfloat x, y, z;
        friend std::istream& operator>>(std::istream& in, point3& xyz) {
            in >> xyz.x >> xyz.y >> xyz.z;
            return in;
        }
    };


    class Coordinates
    {
    public:
        std::vector<point3> values;
        std::vector<point3> getValues() {
            return values;
        }


        friend std::istream& operator>>(std::istream& in, Coordinates& c) {
            // skip comment lines
            point3 xyz;
            while (in >> xyz) {
                c.values.push_back(xyz);
            }
            return in;
        }


    };
}