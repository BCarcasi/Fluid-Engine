

#include "jet.h"
#include <pystring/pystring.h>

#include "clara_utils.h"
#include "io_utils.h"
#include "clara.h"



#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include "obj_simplifier.h"



using namespace jet;

const std::string kSpherical = "spherical";
const std::string kSph = "sph";
const std::string kZhuBridson = "zhu_bridson";
const std::string kAnisotropic = "anisotropic";

double sSphCutOffDensity = 0.5;
double sZhuBridsonCutOffThreshold = 0.25;
double sAnisoCutOffDensity = 0.5;
double sAnisoPositionSmoothingFactor = 0.5;
size_t sAnisoMinNumNeighbors = 25;
int animationLength = 100;

#define APP_NAME "hybrid_liquid_sim"
std::string outputDir = APP_NAME "_output2";

void printInfo(const Size3& resolution, const BoundingBox3D& domain,
    const Vector3D& gridSpacing, size_t numberOfParticles,
    const std::string& method) {
    printf("Resolution: %zu x %zu x %zu\n", resolution.x, resolution.y,
        resolution.z);
    printf("Domain: [%f, %f, %f] x [%f, %f, %f]\n", domain.lowerCorner.x,
        domain.lowerCorner.y, domain.lowerCorner.z, domain.upperCorner.x,
        domain.upperCorner.y, domain.upperCorner.z);
    printf("Grid spacing: [%f, %f, %f]\n", gridSpacing.x, gridSpacing.y,
        gridSpacing.z);
    printf("Number of particles: %zu\n", numberOfParticles);
    printf("Reconstruction method: %s\n", method.c_str());
}

void triangulateAndSave(const ScalarGrid3& sdf,
    const std::string& objFilename) {
    TriangleMesh3 mesh;
    FaceCenteredGrid3Ptr blank = FaceCenteredGrid3Ptr();
    marchingCubes(sdf.constDataAccessor(), sdf.gridSpacing(), sdf.dataOrigin(),
        &mesh, blank, 0.0, kDirectionAll);

    std::ofstream file(objFilename.c_str());
    if (file) {
        printf("Writing %s...\n", objFilename.c_str());
        mesh.writeObj(&file);
        file.close();
    }
    else {
        printf("Cannot write file %s.\n", objFilename.c_str());
        exit(EXIT_FAILURE);
    }
}

void particlesToObj(const Array1<Vector3D>& positions, const Size3& resolution,
    const Vector3D& gridSpacing, const Vector3D& origin,
    double kernelRadius, const std::string& method,
    const std::string& objFilename) {
    PointsToImplicit3Ptr converter;
    if (method == kSpherical) {
        converter =
            std::make_shared<SphericalPointsToImplicit3>(kernelRadius, false);
    }
    else if (method == kSph) {
        converter = std::make_shared<SphPointsToImplicit3>(
            kernelRadius, sSphCutOffDensity, false);
    }
    else if (method == kZhuBridson) {
        converter = std::make_shared<ZhuBridsonPointsToImplicit3>(
            kernelRadius, sZhuBridsonCutOffThreshold, false);
    }
    else {
        converter = std::make_shared<AnisotropicPointsToImplicit3>(
            kernelRadius, sAnisoCutOffDensity, sAnisoPositionSmoothingFactor,
            sAnisoMinNumNeighbors, false);
    }

    VertexCenteredScalarGrid3 sdf(resolution, gridSpacing, origin);
    printInfo(resolution, sdf.boundingBox(), gridSpacing, positions.size(),
        method);

    converter->convert(positions, &sdf);

    triangulateAndSave(sdf, objFilename);


}

bool LoadOBJ(const std::string& fileName,
    std::vector< MeshDecimation::Vec3<MeshDecimation::Float> >& points,
    std::vector< MeshDecimation::Vec3<int> >& triangles);
bool SaveOBJ(const std::string& fileName,
    const std::vector< MeshDecimation::Vec3<MeshDecimation::Float> >& points,
    const std::vector< MeshDecimation::Vec3<int> >& triangles);

bool LoadOBJ(const std::string& fileName,
    std::vector< MeshDecimation::Vec3<MeshDecimation::Float> >& points,
    std::vector< MeshDecimation::Vec3<int> >& triangles)
{
    const char ObjDelimiters[] = " /";
    const unsigned int BufferSize = 1024;
    FILE* fid = fopen(fileName.c_str(), "r");

    if (fid)
    {
        char buffer[BufferSize];
        MeshDecimation::Vec3<int> ip;
        MeshDecimation::Vec3<int> in;
        MeshDecimation::Vec3<int> it;
        char* pch;
        char* str;
        size_t nn = 0;
        size_t nt = 0;
        MeshDecimation::Vec3<MeshDecimation::Float> x;
        while (!feof(fid))
        {
            if (!fgets(buffer, BufferSize, fid))
            {
                break;
            }
            else if (buffer[0] == 'v')
            {
                if (buffer[1] == ' ')
                {
                    str = buffer + 2;
                    for (int k = 0; k < 3; ++k)
                    {
                        pch = strtok(str, " ");
                        if (pch) x[k] = static_cast<MeshDecimation::Float>(atof(pch));
                        else
                        {
                            return false;
                        }
                        str = NULL;
                    }
                    points.push_back(x);
                }
                else if (buffer[1] == 'n')
                {
                    ++nn;
                }
                else if (buffer[1] == 't')
                {
                    ++nt;
                }
            }
            else if (buffer[0] == 'f')
            {

                str = buffer + 2;
                for (int k = 0; k < 3; ++k)
                {
                    pch = strtok(str, ObjDelimiters);
                    if (pch) ip[k] = atoi(pch) - 1;
                    else
                    {
                        return false;
                    }
                    str = NULL;
                    if (nt > 0)
                    {
                        pch = strtok(NULL, ObjDelimiters);
                        if (pch)  it[k] = atoi(pch) - 1;
                        else
                        {
                            return false;
                        }
                    }
                    if (nn > 0)
                    {
                        pch = strtok(NULL, ObjDelimiters);
                        if (pch)  in[k] = atoi(pch) - 1;
                        else
                        {
                            return false;
                        }
                    }
                }
                triangles.push_back(ip);
            }
        }
        fclose(fid);
    }
    else
    {
        std::cout << "File not found" << std::endl;
        return false;
    }
    return true;
}
bool SaveOBJ(const std::string& fileName,
    const std:: vector< MeshDecimation::Vec3<MeshDecimation::Float> >& points,
    const std::vector< MeshDecimation::Vec3<int> >& triangles)
{
    std::cout << "Saving " << fileName << std::endl;
    std::ofstream fout(fileName.c_str());
    if (fout.is_open())
    {
        const size_t nV = points.size();
        const size_t nT = triangles.size();
        for (size_t v = 0; v < nV; v++)
        {
            fout << "v " << points[v][0] << " "
                << points[v][1] << " "
                << points[v][2] << std::endl;
        }
        for (size_t f = 0; f < nT; f++)
        {
            fout << "f " << triangles[f][0] + 1 << " "
                << triangles[f][1] + 1 << " "
                << triangles[f][2] + 1 << std::endl;
        }
        fout.close();
        return true;
    }
    return false;
}

void CallBack(const char* msg)
{
    std::cout << msg;
}


int main(int argc, char* argv[]) {
    bool showHelp = false;
    Size3 resolution(50*3, 50*2, 50*1.5);
    Vector3D gridSpacing(0.02, 0.02, 0.02);
    Vector3D origin;
    std::string method = "sph";
    double kernelRadius = 0.04;

    std::string strResolution;
    std::string strGridSpacing;
    std::string strOrigin;
    std::string strMethod;

    // Parsing
    auto parser =
        clara::Help(showHelp) |
        clara::Opt(strResolution, "resolution")["-r"]["--resolution"](
            "grid resolution in CSV format (default is 100,100,100)") |
        clara::Opt(strGridSpacing, "gridSpacing")["-g"]["--grid_spacing"](
            "grid spacing in CSV format (default is 0.01,0.01,0.01)") |
        clara::Opt(strOrigin, "origin")["-n"]["--origin"](
            "domain origin in CSV format (default is 0,0,0)") |
        clara::Opt(method, "method")["-m"]["--method"](
            "spherical, sph, zhu_bridson, and anisotropic "
            "followed by optional method-dependent parameters (default is "
            "anisotropic)") |
        clara::Opt(kernelRadius, "kernelRadius")["-k"]["--kernel"](
            "interpolation kernel radius (default is 0.2)");



    auto result = parser.parse(clara::Args(argc, argv));
    if (!result) {
        std::cerr << "Error in command line: " << result.errorMessage() << '\n';
        exit(EXIT_FAILURE);
    }

    if (showHelp) {
        std::cout << toString(parser) << '\n';
        exit(EXIT_SUCCESS);
    }

    // Resolution
    if (!strResolution.empty()) {
        std::vector<std::string> tokens;
        pystring::split(strResolution, tokens, ",");

        if (tokens.size() == 1) {
            resolution.x = resolution.y = resolution.z =
                static_cast<size_t>(atoi(strResolution.c_str()));
        }
        else if (tokens.size() == 3) {
            resolution.x = static_cast<size_t>(atoi(tokens[0].c_str()));
            resolution.y = static_cast<size_t>(atoi(tokens[1].c_str()));
            resolution.z = static_cast<size_t>(atoi(tokens[2].c_str()));
        }
    }

    // Grid spacing
    if (!strGridSpacing.empty()) {
        std::vector<std::string> tokens;
        pystring::split(strGridSpacing, tokens, ",");

        if (tokens.size() == 1) {
            gridSpacing.x = gridSpacing.y = gridSpacing.z =
                atof(strGridSpacing.c_str());
        }
        else if (tokens.size() == 3) {
            gridSpacing.x = atof(tokens[0].c_str());
            gridSpacing.y = atof(tokens[1].c_str());
            gridSpacing.z = atof(tokens[2].c_str());
        }
    }

    // Origin
    if (!strOrigin.empty()) {
        std::vector<std::string> tokens;
        pystring::split(strOrigin, tokens, ",");

        if (tokens.size() == 1) {
            origin.x = origin.y = origin.z = atof(strOrigin.c_str());
        }
        else if (tokens.size() == 3) {
            origin.x = atof(tokens[0].c_str());
            origin.y = atof(tokens[1].c_str());
            origin.z = atof(tokens[2].c_str());
        }
    }

    // Method
    if (!strMethod.empty()) {
        std::vector<std::string> tokens;
        pystring::split(strMethod, tokens, ",");

        method = tokens[0];

        if (method == kSpherical) {
            // No other options accepted
        }
        else if (method == kSph) {
            if (tokens.size() > 1) {
                sSphCutOffDensity = atof(tokens[1].c_str());
            }
        }
        else if (method == kZhuBridson) {
            if (tokens.size() > 1) {
                sZhuBridsonCutOffThreshold = atof(tokens[1].c_str());
            }
        }
        else if (method == kAnisotropic) {
            if (tokens.size() > 1) {
                sAnisoCutOffDensity = atof(tokens[1].c_str());
            }
            if (tokens.size() > 2) {
                sAnisoPositionSmoothingFactor = atof(tokens[2].c_str());
            }
            if (tokens.size() > 3) {
                sAnisoMinNumNeighbors =
                    static_cast<size_t>(atoi(tokens[3].c_str()));
            }
        }
        else {
            fprintf(stderr, "Unknown method %s.\n", method.c_str());
            std::cout << toString(parser) << '\n';
            exit(EXIT_FAILURE);
        }
    }

    for (int x = 0; x < animationLength; x++) {
        char basename[256];
        snprintf(basename, sizeof(basename), "frame_%06d.xyz", x);
        std::string filename = pystring::os::path::join(outputDir, basename);
        std::string inputFilename(filename);
        snprintf(basename, sizeof(basename), "frame_%06d.obj", x);
        filename = pystring::os::path::join(outputDir, basename);
        std::string outputFilename(filename);

        if (inputFilename.empty() || outputFilename.empty()) {
            std::cout << toString(parser) << '\n';
            exit(EXIT_FAILURE);
        }


        // Read particle positions
        Array1<Vector3D> positions;
        if (!readPositions(inputFilename, positions)) {
            printf("Cannot read file %s.\n", inputFilename.c_str());
            exit(EXIT_FAILURE);
        }

        // Run marching cube and save it to the disk
        particlesToObj(positions, resolution, gridSpacing, origin, kernelRadius,
            method, outputFilename);

        /*
        std::vector< MeshDecimation::Vec3<MeshDecimation::Float> > points;
        std::vector< MeshDecimation::Vec3<int> > triangles;


        LoadOBJ(outputFilename, points, triangles);
                MeshDecimation::MeshDecimator myMDecimator;
        myMDecimator.SetCallBack(&CallBack);
        myMDecimator.Initialize(points.size(), triangles.size(), &points[0], &triangles[0]);
        myMDecimator.Decimate(triangles.size()/2,
            points.size()/2,
            1.0);

        // allocate memory for decimated mesh
        std::vector< MeshDecimation::Vec3<MeshDecimation::Float> > decimatedPoints;
        std::vector< MeshDecimation::Vec3<int> > decimatedtriangles;
        decimatedPoints.resize(myMDecimator.GetNVertices());
        decimatedtriangles.resize(myMDecimator.GetNTriangles());

        // retreive decimated mesh
        myMDecimator.GetMeshData(&decimatedPoints[0], &decimatedtriangles[0]);

        SaveOBJ(outputFilename, decimatedPoints, decimatedtriangles);

        */

    }


    
    return EXIT_SUCCESS;
}

