#pragma once

#include <GLFW/glfw3.h>
#include <iostream>
#include<fstream>
namespace {
  

    const int windowWidth = 640, windowHeight = 480;


    inline int graphicsTest(void)
    {

        std::ifstream inputfile("simpleSpring.txt");
        double framerate;
        int frames = 0;
        int springs = 0;
        inputfile >> springs >> framerate >> frames;
        

        std::istream::streampos p = inputfile.tellg(); // or, in C++11: auto p = f.tellg();
        // f.clear() here if there's a possibility that the stream is in a bad state
       
        

        GLFWwindow* window;

        /* Initialize the library */
        if (!glfwInit())
            return -1;

        /* Create a windowed mode window and its OpenGL context */
        window = glfwCreateWindow(windowWidth, windowHeight, "Simple Spring", NULL, NULL);
        if (!window)
        {
            glfwTerminate();
            return -1;
        }

        /* Make the window's context current */
        glfwMakeContextCurrent(window);
        glfwSwapInterval(3);
        /* Loop until the user closes the window */
        while (!glfwWindowShouldClose(window))
        {
            /* Render here */


            glClearColor(1.0f, 1.0f, 1.0f, 0);
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
            glLoadIdentity();

            if (frames < 0) {
                inputfile.clear();
                inputfile.seekg(0, std::ios::beg);
                inputfile >> springs >> framerate >> frames;


            }
              
            frames--;
            float x1  = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0;
            inputfile >> x1 >> y1;
            for (int i = 1; i < (springs -1) * 11; i++) {
                inputfile >> x2 >> y2;
                glColor3f(0.0f, 0.0f, 0.0f);
                glBegin(GL_LINES);
                glVertex2f(x1 / 10.0f, y1 / 3.5f + 1.0f);
                glVertex2f(x2 / 10.0f, y2 / 3.5f + 1.0f);
                glEnd();
                x1 = x2; y1 = y2;
               

            }
          
            

            /* Swap front and back buffers */
            glfwSwapBuffers(window);

            /* Poll for and process events */
            glfwPollEvents();
        }

        glfwTerminate();
        return 0;
    }
}