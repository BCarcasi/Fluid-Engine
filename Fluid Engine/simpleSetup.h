#pragma once

#include "simple_spring.h"

#include <iostream>



namespace jet {

    inline int generateFile() {
        int springs, frames;
        float windX, windY;
        std::cout << "Enter a number from 1-20 for number of Springs: ";
        std::cin >> springs;
        while (springs <= 0 || springs > 20) {
            std::cout << "Enter a number from 1-20 for number of Springs: ";
            std::cin >> springs;
        }
        std::cout << "Enter a nonzero number of Frames (The simulation runs at 60 frames a second): ";
        std::cin >> frames;
        std::cout << "set the wind in the x direction: ";
        std::cin >> windX;
        std::cout << "set the wind in the y direction: ";
        std::cin >> windY;




        Array1<double> x;
        Array1<double> y;


        SimpleMassSpringAnimation anim;
        anim.makeChain(springs);
        anim.wind = std::make_shared<ConstantVectorField3>(Vector3D(windX, windY, 0.0));
        anim.constraints.push_back(SimpleMassSpringAnimation::Constraint{ 0, Vector3D(), Vector3D() });
        anim.exportStates(x, y);
        anim.openFile("simpleSpring.txt", std::to_string(1.0 / 60.0) + " " + std::to_string(frames));



        for (Frame frame(0, 1.0 / 60.0); frame.index < frames; frame.advance())
        {
            anim.update(frame);
            anim.exportStates(x, y);

            anim.exportToFile();


        }

        return 1;

    }

}