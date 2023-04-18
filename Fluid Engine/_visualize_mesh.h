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
#include "constants.h"
#include "texture.h"
#include "vboindexer.h"
#include <numeric>
#include <tbb/parallel_for.h>


#define APP_NAME "hybrid_liquid_sim"





const int AnimationLength = 100;
float currentTime = 0;
float lastTime = 0;
float framerate = 60.0;
float frame = 0;



int main(void)
{
	// Initialise GLFW
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
	window = glfwCreateWindow(1024, 768, "meshRenderer", NULL, NULL);
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
	glfwSetCursorPos(window, 1024 / 2, 768 / 2);

	// Dark blue background
	glClearColor(0.8f, 0.8f, 0.8f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	// Cull triangles which normal is not towards the camera
	//glEnable(GL_CULL_FACE); // Not this time !

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	GLuint programID = LoadShaders("TransformVertexShader.vertexshader", "TextureFragmentShader.fragmentshader");

	// Get a handle for our "MVP" uniform
	GLuint MatrixID = glGetUniformLocation(programID, "MVP");
	GLuint ViewMatrixID = glGetUniformLocation(programID, "V");
	GLuint ModelMatrixID = glGetUniformLocation(programID, "M");

	// Load the texture
	GLuint Texture = loadDDS("uvmap.DDS");

	// Get a handle for our "myTextureSampler" uniform
	GLuint TextureID = glGetUniformLocation(programID, "myTextureSampler");

	// Read our .obj file
	std::vector<std::vector<glm::vec3>> vertices;
	std::vector < std::vector<glm::vec2>> uvs;
	std::vector < std::vector<glm::vec3>> normals;
	std::string outputDir = APP_NAME "_output2";
	char basename[256];


	for (int i = 0; i < AnimationLength; i++) {
		snprintf(basename, sizeof(basename), "frame_%06d.obj", i);
		std::string filename = pystring::os::path::join(outputDir, basename);
		std::vector<glm::vec3> verticesBase; std::vector<glm::vec2> uvsBase; std::vector<glm::vec3> normalsBase;
		bool res = loadOBJ(filename.c_str(), verticesBase, uvsBase, normalsBase);



		vertices.push_back(verticesBase); uvs.push_back(uvsBase); normals.push_back(normalsBase);
	}





	//bool res = loadOBJ("suzanne.obj", vertices, uvs, normals);

	std::vector < std::vector<unsigned int>> indices;
	std::vector < std::vector<glm::vec3>> indexed_vertices;
	std::vector < std::vector<glm::vec2>> indexed_uvs;
	std::vector < std::vector<glm::vec3>> indexed_normals;



	for (int i = 0; i < AnimationLength; i++) {
		std::vector<unsigned int> indicesBuffer;
		std::vector<glm::vec3> indexed_verticesBuffer;
		std::vector<glm::vec2> indexed_uvsBuffer;
		std::vector<glm::vec3> indexed_normalsBuffer;
		indexVBO(vertices[i], uvs[i], normals[i], indicesBuffer, indexed_verticesBuffer, indexed_uvsBuffer, indexed_normalsBuffer);


		indices.push_back(indicesBuffer); indexed_vertices.push_back(indexed_verticesBuffer); indexed_uvs.push_back(indexed_uvsBuffer); indexed_normals.push_back(indexed_normalsBuffer);
	}




	// Load it into a VBO





	GLuint vertexbuffers[AnimationLength] = { 0 };
	GLuint uvbuffers[AnimationLength] = { 0 };
	GLuint normalbuffers[AnimationLength] = { 0 };
	GLuint elementbuffers[AnimationLength] = { 0 };


	glGenBuffers(AnimationLength, vertexbuffers);
	glGenBuffers(AnimationLength, uvbuffers);
	glGenBuffers(AnimationLength, normalbuffers);
	glGenBuffers(AnimationLength, elementbuffers);






	// Get a handle for our "LightPosition" uniform
	glUseProgram(programID);
	GLuint LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

	// For speed computation
	double lastTime = glfwGetTime();
	int nbFrames = 0;

	// Enable blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	bool first = true;

	do {


		currentTime = glfwGetTime();

		if (currentTime - lastTime > 1.0 / framerate) {
			lastTime = currentTime;

			// Clear the screen
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			// Use our shader
			glUseProgram(programID);



			// Compute the MVP matrix from keyboard and mouse input
			computeMatricesFromInputs(window);
			glm::mat4 ProjectionMatrix = getProjectionMatrix();
			glm::mat4 ViewMatrix = getViewMatrix();
			glm::mat4 ModelMatrix = glm::mat4(1.0);
			glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

			/*
			if (first) {
				for (int i = 0; i < indices[0].size()/3; i++) {
					std::cout << indices[0][i] << " " << indices[0][i + 1] << " " << indices[0][i + 2] << std::endl;
				}
			}
			*/

			//std::cout << "end" << endl;

			std::vector<vec3> sortArray(indices[floor(frame)].size() / 3);
			for (int i = 0; i < sortArray.size(); i++) {
				sortArray[i].x = indices[floor(frame)][i * 3];
				sortArray[i].y = indices[floor(frame)][i * 3 + 1];
				sortArray[i].z = indices[floor(frame)][i * 3 + 2];
			}

			std::sort(sortArray.begin(), sortArray.end(), [&](vec3 x, vec3 y) {
				vec4 temp1 = MVP * vec4(indexed_vertices[floor(frame)][x.x], 1);
				vec4 temp2 = MVP * vec4(indexed_vertices[floor(frame)][y.x], 1);
				return(temp1.z / temp1.w) > (temp2.z / temp2.w);
				});


			for (int i = 0; i < sortArray.size(); i++) {
				indices[floor(frame)][i * 3] = sortArray[i].x;
				indices[floor(frame)][i * 3 + 1] = sortArray[i].y;
				indices[floor(frame)][i * 3 + 2] = sortArray[i].z;
			}
			/*
						if (first) {
							for (int i = 0; i < indices[0].size()/3; i++) {
								std::cout << indices[0][i*3]<< " " << indices[0][i*3 + 1] << " " << indices[0][i*3 +2] << std::endl;
							}
						}
						first = false;
						*/



			for (int i = 0; i < AnimationLength; i++) {
				glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[i]);
				glBufferData(GL_ARRAY_BUFFER, indexed_vertices[i].size() * sizeof(glm::vec3), &indexed_vertices[i][0], GL_STATIC_DRAW);



				glBindBuffer(GL_ARRAY_BUFFER, uvbuffers[i]);
				glBufferData(GL_ARRAY_BUFFER, indexed_uvs[i].size() * sizeof(glm::vec2), &indexed_uvs[i][0], GL_STATIC_DRAW);


				glBindBuffer(GL_ARRAY_BUFFER, normalbuffers[i]);
				glBufferData(GL_ARRAY_BUFFER, indexed_normals[i].size() * sizeof(glm::vec3), &indexed_normals[i][0], GL_STATIC_DRAW);

				// Generate a buffer for the indices as well

				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffers[i]);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices[i].size() * sizeof(unsigned int), &indices[i][0], GL_STATIC_DRAW);

			}













			// Send our transformation to the currently bound shader, 
			// in the "MVP" uniform
			glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
			glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
			glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);

			glm::vec3 lightPos = glm::vec3(4, 4, 4);
			glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

			// Bind our texture in Texture Unit 0
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D, Texture);
			// Set our "myTextureSampler" sampler to use Texture Unit 0
			glUniform1i(TextureID, 0);

			// 1rst attribute buffer : vertices
			glEnableVertexAttribArray(0);
			glBindBuffer(GL_ARRAY_BUFFER, vertexbuffers[(int)floor(frame)]);
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
			glBindBuffer(GL_ARRAY_BUFFER, uvbuffers[(int)floor(frame)]);
			glVertexAttribPointer(
				1,                                // attribute
				2,                                // size
				GL_FLOAT,                         // type
				GL_FALSE,                         // normalized?
				0,                                // stride
				(void*)0                          // array buffer offset
			);

			// 3rd attribute buffer : normals
			glEnableVertexAttribArray(2);
			glBindBuffer(GL_ARRAY_BUFFER, normalbuffers[(int)floor(frame)]);
			glVertexAttribPointer(
				2,                                // attribute
				3,                                // size
				GL_FLOAT,                         // type
				GL_FALSE,                         // normalized?
				0,                                // stride
				(void*)0                          // array buffer offset
			);

			// Index buffer
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffers[(int)floor(frame)]);

			// Draw the triangles !
			glDrawElements(
				GL_TRIANGLES,      // mode
				indices[floor(frame)].size(),    // count
				GL_UNSIGNED_INT, // type
				(void*)0           // element array buffer offset
			);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			glDisableVertexAttribArray(2);

			frame += 1.0f;

			if (floor(frame) >= AnimationLength) {
				frame = 0;
			}

			// Swap buffers
			glfwSwapBuffers(window);
			glfwPollEvents();
		}

	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	// Cleanup VBO and shader
	for (int i = 0; i < AnimationLength; i++) {
		glDeleteBuffers(1, &vertexbuffers[i]);
		glDeleteBuffers(1, &uvbuffers[i]);
		glDeleteBuffers(1, &normalbuffers[i]);
		glDeleteBuffers(1, &elementbuffers[i]);
	}

	glDeleteProgram(programID);
	glDeleteTextures(1, &Texture);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();

	return 0;

}