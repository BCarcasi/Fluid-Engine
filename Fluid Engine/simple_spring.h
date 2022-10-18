#include "physics_animation.h"
#include "vector3.h"
#include "constant_vector_field3.h"
#include "array1.h"
#include <iostream>
#include <fstream>


using namespace jet;

class SimpleMassSpringAnimation : public PhysicsAnimation
{
public:
    struct Edge
    {
        size_t first;
        size_t second;
    };

    struct Constraint
    {
        size_t pointIndex;
        Vector3D fixedPosition;
        Vector3D fixedVelocity;
    };

    std::vector<Vector3D> positions;
    std::vector<Vector3D> velocities;
    std::vector<Vector3D> forces;
    std::vector<Edge> edges;

    double mass = 1.0;
    Vector3D gravity = Vector3D(0.0, -9.8, 0.0);
    double stiffness = 500.0;
    double restLength = 1.0;
    double dampingCoefficient = 1.0;
    double dragCoefficient = 0.1;

    double floorPositionY = -7.0;
    double restitutionCoefficient = 0.3;

    VectorField3Ptr wind;

    std::vector<Constraint> constraints;

    const char * fileName;

    SimpleMassSpringAnimation() {}

    void makeChain(size_t numberOfPoints)
    {
        if (numberOfPoints == 0)
        {
            return;
        }

        size_t numberOfEdges = numberOfPoints - 1;

        positions.resize(numberOfPoints);
        velocities.resize(numberOfPoints);
        forces.resize(numberOfPoints);
        edges.resize(numberOfEdges);

        for (size_t i = 0; i < numberOfPoints; ++i)
        {
            positions[i].x = -static_cast<double>(i);
        }

        for (size_t i = 0; i < numberOfEdges; ++i)
        {
            edges[i] = Edge{ i, i + 1 };
        }
    }

    void exportStates(Array1<double>& x, Array1<double>& y) const
    {
        x.resize(positions.size());
        y.resize(positions.size());

        for (size_t i = 0; i < positions.size(); ++i)
        {
            x[i] = positions[i].x;
            y[i] = positions[i].y;
        }
    }

    
    void openFile(const char* fileName, std::string initString = "") {
        this->fileName = fileName;
        std::ofstream outfile;
        outfile.open(fileName, std::fstream::out | std::ofstream::trunc);
        outfile << positions.size() << " ";
        if (!initString.empty()) {
            outfile << initString << std::endl;
        }
        else {
            outfile << std::to_string(1.0 / 60.0) << " 360" << std::endl;
        }
        outfile.close();
    }

    void exportToFile() const {


        std::ofstream outfile;
        outfile.open(fileName, std::fstream::out | std::ofstream::app);
        for (int i = 1; i < positions.size(); i++) {
            Vector2F springBase(positions[i - 1].x, positions[i - 1].y);
            outfile << springBase.x << " " << springBase.y << " ";



            Vector2F springVisual(positions[i].y - positions[i - 1].y, -(positions[i].x - positions[i - 1].x));
            Vector2F springSpacing(positions[i].x - positions[i - 1].x, positions[i].y - positions[i - 1].y);
            springVisual.normalize();
            springVisual = springVisual / 4.0f;
            for (int j = 1; j < 11; j++){
            Vector2F finalSpring = springBase + 1.0f / 11.0f * j * springSpacing + springVisual * static_cast<float>(((j % 2) * 2) - 1);
            outfile << finalSpring.x << " " << finalSpring.y << " ";
            }
        }
        outfile << "\n";
        outfile.close();

    }

    

protected:

    void onAdvanceTimeStep(double timeIntervalInSeconds) override
    {
        size_t numberOfPoints = positions.size();
        size_t numberOfEdges = edges.size();

        // Compute forces
        for (size_t i = 0; i < numberOfPoints; ++i)
        {
            // Gravity force
            forces[i] = mass * gravity;

            // Air drag force
            Vector3D relativeVel = velocities[i];
            if (wind != nullptr)
            {
                relativeVel -= wind->sample(positions[i]);
            }
            forces[i] += -dragCoefficient * relativeVel;
        }

        for (size_t i = 0; i < numberOfEdges; ++i)
        {
            size_t pointIndex0 = edges[i].first;
            size_t pointIndex1 = edges[i].second;

            // Compute spring force
            Vector3D pos0 = positions[pointIndex0];
            Vector3D pos1 = positions[pointIndex1];
            Vector3D r = pos0 - pos1;
            double distance = r.length();
            if (distance > 0.0)
            {
                Vector3D force = -stiffness * (distance - restLength) * r.normalized();
                forces[pointIndex0] += force;
                forces[pointIndex1] -= force;
            }

            // Add damping force
            Vector3D vel0 = velocities[pointIndex0];
            Vector3D vel1 = velocities[pointIndex1];
            Vector3D relativeVel0 = vel0 - vel1;
            Vector3D damping = -dampingCoefficient * relativeVel0;
            forces[pointIndex0] += damping;
            forces[pointIndex1] -= damping;
        }

        // Update states
        for (size_t i = 0; i < numberOfPoints; ++i)
        {
            // Compute new states
            Vector3D newAcceleration = forces[i] / mass;
            Vector3D newVelocity = velocities[i] + timeIntervalInSeconds * newAcceleration;
            Vector3D newPosition = positions[i] + timeIntervalInSeconds * newVelocity;

            // Collision
            if (newPosition.y < floorPositionY)
            {
                newPosition.y = floorPositionY;

                if (newVelocity.y < 0.0)
                {
                    newVelocity.y *= -restitutionCoefficient;
                    newPosition.y += timeIntervalInSeconds * newVelocity.y;
                }
            }

            // Update states
            velocities[i] = newVelocity;
            positions[i] = newPosition;
        }

        // Apply constraints
        for (size_t i = 0; i < constraints.size(); ++i)
        {
            size_t pointIndex = constraints[i].pointIndex;
            positions[pointIndex] = constraints[i].fixedPosition;
            velocities[pointIndex] = constraints[i].fixedVelocity;
        }
    }
};