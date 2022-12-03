
#include <pystring/pystring.h>
#include <cmath>
#include <gl/glew.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <gl/GLUT.h>


#include <iostream>
#include <vector>

#include "constants.h"
#include <fstream>
#include "xyz_extractor.h"
#include <algorithm>

#define APP_NAME "sph_sim"

// angle of rotation for the camera direction
float angle = 0.0;
float angleY = -0.5f;
// actual vector representing the camera's direction

// XZ position of the camera
float camX = 1.5f, camY = 0.0f, z = 1.0f;
float lx = sin(angle) * cos(angleY), ly = sin(angleY), lz = -cos(angle) * cos(angleY), dis = 5.0f;
bool freeze = false;
int numOfFrames = 100;


int numberOfPoints = 0, frame = 0, count = 0, framesBetweenSwitch = 7;



jet::Coordinates loadXYZFromFile(const std::string& rootDir, int frameCnt) {
	char basename[256];
	snprintf(basename, sizeof(basename), "frame_%06d.xyz", frameCnt);
	std::string filename = pystring::os::path::join(rootDir, basename);
	std::fstream file(filename.c_str());
	jet::Coordinates coordinates;
	if (file) {
		file >> coordinates;
		file.close();
	}
	return coordinates;
}

std::vector<jet::point3> Vertices;

struct Point
{
	float x, y, z;
	unsigned char r, g, b, a;
};
std::vector< Point > points;
std::vector<std::vector<Point>> arrayOfPoints;

bool basicCompareOpenGLPoint(Point p1, Point p2) {

	return p1.y < p2.y;
}

bool trueCompareOpenGLPoint(Point p1, Point p2) {

	return std::pow(p1.x - (camX - lx * dis), 2) + std::pow(p1.y - (camY - ly * dis), 2) + std::pow(p1.z - (z - lz * dis), 2)
		< std::pow(p2.x - (camX - lx * dis), 2) + std::pow(p2.y - (camY - ly * dis), 2) + std::pow(p2.z - (z - lz * dis), 2);
}

bool testTrueCompareOpenGLPoint(Point p1, Point p2) {
	return ((p1.x - (camX - lx * dis)) * sin(angle) * cos(angleY) + (p1.y - (camY - ly * dis)) * -sin(angleY) + (p1.z - (z - lz * dis)) * cos(angle) * cos(angleY))
		> ((p2.x - (camX - lx * dis)) * sin(angle) * cos(angleY) + (p2.y - (camY - ly * dis)) * -sin(angleY) + (p2.z - (z - lz * dis)) * cos(angle) * cos(angleY));
}



void changeSize(int w, int h) {

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (h == 0)
		h = 1;
	float ratio = w * 1.0 / h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, w, h);

	// Set the correct perspective.
	gluPerspective(45.0f, ratio, 0.1f, 100.0f);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
}

void processNormalKeys(unsigned char key, int x, int y) {

	if (key == 27)
		exit(0);
	if (key == 'e' && dis > 0.5)
		dis -= 0.1f;
	if (key == 'q')
		dis += 0.1f;
	if (key == 'w') {
		camX += sin(angle) * 0.2f;
		z -= cos(angle) * 0.2f;
	}
	if (key == 'a') {
		camX -= cos(angle) * 0.2f;
		z -= sin(angle) * 0.2f;
	}
	if (key == 's') {
		camX -= sin(angle) * 0.2f;
		z += cos(angle) * 0.2f;
	}
	if (key == 'd') {
		camX += cos(angle) * 0.2f;
		z += sin(angle) * 0.2f;
	}
	if (key == 'z') camY += 0.2f;
	if (key == 'x') camY -= 0.2f;
	if (key == 'k' && framesBetweenSwitch > 1) framesBetweenSwitch -= 1;
	if (key == 'j') framesBetweenSwitch += 1;
	if (key == 'n') {
		for (int i = 0; i < arrayOfPoints.size(); i++) {
			points = arrayOfPoints.at(i);
			std::sort(points.begin(), points.end(), basicCompareOpenGLPoint);
			std::cout << i << std::endl;
		}
	}
	if (key == 'b') {
		for (int i = 0; i < arrayOfPoints.size(); i++) {
			points = arrayOfPoints.at(i);
			std::sort(points.begin(), points.end(), trueCompareOpenGLPoint);
			std::cout << i << std::endl;
		}
	}
	if (key == 'm') std::cout << camX - lx * dis << " " << camY - ly * dis << camX - lx * dis << " " << z - lz * dis << " "
		<< cos(angle) << " " << sin(angle) << " " << sin(angleY) << std::endl;
	if (key == 'i') freeze = (1 - freeze);
	if (key == 'r') {
		angle = 0.0;
		angleY = -0.5f;
		// actual vector representing the camera's direction

		// XZ position of the camera
		camX = 1.5f, camY = 0.0f, z = 1.0f;
		lx = sin(angle) * cos(angleY), ly = sin(angleY), lz = -cos(angle) * cos(angleY), dis = 5.0f;
		freeze = false;
		framesBetweenSwitch = 7;
	}

}



void renderScene(void) {

	// Clear Color and Depth Buffers

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformations
	glLoadIdentity();
	// Set the camera
	gluLookAt(camX - lx * dis, camY - ly * dis, z - lz * dis,
		camX, camY, z,
		0.0f, 1.0f, 0.0f);

	// Draw ground
	glColor3f(0.9f, 0.9f, 0.9f);
	glBegin(GL_QUADS);
	glVertex3f(-100.0f, 0.0f, -100.0f);
	glVertex3f(-100.0f, 0.0f, 100.0f);
	glVertex3f(100.0f, 0.0f, 100.0f);
	glVertex3f(100.0f, 0.0f, -100.0f);
	glEnd();


	points = arrayOfPoints.at(frame);


	glColor3ub(255, 255, 255);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glVertexPointer(3, GL_FLOAT, sizeof(Point), &points[0].x);
	glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(Point), &points[0].r);
	glPointSize(2.0 + pow(3.0f / (dis + 0.8), 2));
	glDrawArrays(GL_POINTS, 0, points.size());
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glFlush();

	glutSwapBuffers();
	if (count == 0 && !freeze) {
		frame = (frame + 1) % numOfFrames;
	}
	count = (count + 1) % framesBetweenSwitch;
}

void processSpecialKeys(int key, int xx, int yy) {

	float fraction = 0.1f;

	switch (key) {
	case GLUT_KEY_LEFT:
		angle += 0.02f;
		lx = sin(angle) * cos(angleY);
		lz = -cos(angle) * cos(angleY);
		break;
	case GLUT_KEY_RIGHT:
		angle -= 0.02f;
		lx = sin(angle) * cos(angleY);
		lz = -cos(angle) * cos(angleY);
		break;
	case GLUT_KEY_UP:
		if (angleY > -0.99f * jet::kPiF / 2.0f) {
			angleY -= 0.02f;
		}
		lx = sin(angle) * cos(angleY);
		ly = sin(angleY);
		lz = -cos(angle) * cos(angleY);
		break;
	case GLUT_KEY_DOWN:
		if (angleY < -0.1f) {
			angleY += 0.02f;
		}
		lx = sin(angle) * cos(angleY);
		ly = sin(angleY);
		lz = -cos(angle) * cos(angleY);
		break;
	case GLUT_ACTIVE_SHIFT:
		dis -= 0.1f;
		break;
	case GLUT_ACTIVE_CTRL:
		dis += 0.1f;
		break;
	}
}



int main(int argc, char** argv) {

	// init GLUT and create window

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(640, 640);
	glutCreateWindow("RealTimeRenderingTest");
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);



	// register callbacks
	glutDisplayFunc(renderScene);
	glutReshapeFunc(changeSize);
	glutIdleFunc(renderScene);
	glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);

	// OpenGL init
	glEnable(GL_DEPTH_TEST);

	std::string outputDir = APP_NAME "_output3";
	jet::Coordinates coordinates;
	std::vector<bool> redBlue;
	coordinates = loadXYZFromFile(outputDir, 0);
	for (size_t i = 0; i < coordinates.values.size(); i++) {
		redBlue.push_back(coordinates.values.at(i).y > 0.5);
		//std::cout << i << std::endl;
	}

	for (int i = 0; i < numOfFrames; i++) {
		coordinates = loadXYZFromFile(outputDir, i);
		numberOfPoints = coordinates.values.size();
		std::cout << i << std::endl;
		points.clear();


		for (size_t i = 0; i < coordinates.values.size(); i++)
		{
			Point pt;
			pt.x = coordinates.values.at(i).x;
			pt.y = coordinates.values.at(i).y;
			pt.z = coordinates.values.at(i).z;
			//std::cout << i << std::endl;
			pt.r = 1;
			pt.g = 1;
			pt.b = 255;
			pt.a = 8;
			points.push_back(pt);
		}
		sort(points.begin(), points.end(), basicCompareOpenGLPoint);
		arrayOfPoints.push_back(points);
	}



	// enter GLUT event processing cycle
	glutMainLoop();

	return 1;
}