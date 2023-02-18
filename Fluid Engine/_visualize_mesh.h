#pragma once
#include <pystring/pystring.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#include <GLFW/glfw3.h>
GLFWwindow* window;

// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
using namespace glm;

#include "shader.h"
#include "controls.h"
#include "mesh_extractor.h"


#define APP_NAME "level_set_liquid_sim"




int AnimationLength = 100;







int main(int argc, char** argv) {

	// init GLUT and create window

	if (!glfwInit())
	{
		fprintf(stderr, "Failed to initialize GLFW\n");
		getchar();
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(1024, 768, "Mesh Renderer", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		getchar();
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		getchar();
		glfwTerminate();
		return -1;
	}

	// Ensure we can capture the escape key being pressed below
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	// Hide the mouse and enable unlimited mouvement
	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

	// Set the mouse at the center of the screen
	glfwPollEvents();
	//glfwSetCursorPos(window, 1024 / 2, 768 / 2);

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	// Cull triangles which normal is not towards the camera
	//glEnable(GL_CULL_FACE);

	std::cout << "error 5" << std::endl;

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);


	std::cout << "error 6" << std::endl;
	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("TransformVertexShader.vertexshader", "TextureFragmentShader.fragmentshader");

	std::cout << "error 11" << std::endl;

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");

	std::cout << "error 7" << std::endl;

	// Load the texture
	//GLuint Texture = loadDDS("uvmap.DDS");

	std::cout << "error 8" << std::endl;

	// Get a handle for our "myTextureSampler" uniform
	GLuint TextureID = glGetUniformLocation(programID, "myTextureSampler");
	std::vector<std::vector<glm::vec3>> vertices;
	std::vector<std::vector<glm::vec2>> uvs;
	std::vector<std::vector<glm::vec3>> normals;



	std::string outputDir = APP_NAME "_output3";
	char basename[256];

	std::cout << "error 1" << std::endl;
	for (int i = 0; i < AnimationLength; i++) {
		snprintf(basename, sizeof(basename), "frame_%06d.obj", i);
		std::string filename = pystring::os::path::join(outputDir, basename);
		std::vector<glm::vec3> verticesBase; std::vector<glm::vec2> uvsBase; std::vector<glm::vec3> normalsBase;
		bool res = loadOBJ(filename.c_str(), verticesBase, uvsBase, normalsBase);
		vertices.push_back(verticesBase); uvs.push_back(uvsBase); normals.push_back(normalsBase);
	}

	std::cout << "error 12" << std::endl;
	std::vector<GLuint> vertexbuffer;
	for (int i = 0; i < AnimationLength; i++) {
		GLuint vertexbufferbase;
		glGenBuffers(1, &vertexbufferbase);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferbase);
		glBufferData(GL_ARRAY_BUFFER, vertices[i].size() * sizeof(glm::vec3), &vertices[i][0], GL_STATIC_DRAW);
		vertexbuffer.push_back(vertexbufferbase);
	}

	std::cout << "error 13" << std::endl;
	std::vector<GLuint> uvbuffer;
	for (int i = 0; i < AnimationLength; i++) {
		GLuint uvbufferbase;
		glGenBuffers(1, &uvbufferbase);
		glBindBuffer(GL_ARRAY_BUFFER, uvbufferbase);
		glBufferData(GL_ARRAY_BUFFER, uvs.size() * sizeof(glm::vec2), &uvs[i][0], GL_STATIC_DRAW);
		uvbuffer.push_back(uvbufferbase);
	}
	std::cout << "error 14" << std::endl;
	float frame = 0;
	do {

		// Clear the screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Use our shader
		glUseProgram(programID);

		// Compute the MVP matrix from keyboard and mouse input
		//computeMatricesFromInputs(window);
		//glm::mat4 ProjectionMatrix = getProjectionMatrix();
		//glm::mat4 ViewMatrix = getViewMatrix();
		//glm::mat4 ModelMatrix = glm::mat4(1.0);
		//glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, 
		// in the "MVP" uniform
		//glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

		/*// Bind our texture in Texture Unit 0
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, Texture);
		// Set our "myTextureSampler" sampler to use Texture Unit 0
		glUniform1i(TextureID, 0);
		*/

		// 1rst attribute buffer : vertices
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer[floor(frame)]);
		glVertexAttribPointer(
			0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);

		// 2nd attribute buffer : UVs
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, uvbuffer[floor(frame)]);
		glVertexAttribPointer(
			1,                                // attribute
			2,                                // size
			GL_FLOAT,                         // type
			GL_FALSE,                         // normalized?
			0,                                // stride
			(void*)0                          // array buffer offset
		);

		frame += 0.003f;
		if (floor(frame) > AnimationLength) {
			frame = 0;
		}

		// Draw the triangle !
		glDrawArrays(GL_TRIANGLES, 0, vertices[floor(frame)].size());

		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);

		glEnd();
		glFlush();

		// Swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents();

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	for (int i = 0; i < AnimationLength; i++) {
		glDeleteBuffers(1, &vertexbuffer[i]);
		glDeleteBuffers(1, &uvbuffer[i]);
	}
	glDeleteProgram(programID);
	//glDeleteTextures(1, &Texture);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;
}