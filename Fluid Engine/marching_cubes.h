#ifndef INCLUDE_JET_MARCHING_CUBES_H_
#define INCLUDE_JET_MARCHING_CUBES_H_

#include "array_accessor3.h"
#include "constants.h"


#include "pch.h"

#include "marching_cubes_table.h"
#include "marching_squares_table.h"

#include "bounding_box2.h"
#include "bounding_box3.h"
#include "level_set_utils.h"

#include <array>
#include <unordered_map>



namespace jet {

    //!
    //! \brief      Computes marching cubes and extract triangle mesh from grid.
    //!
    //! This function computes the marching cube algorithm to extract triangle mesh
    //! from the scalar grid field. The triangle mesh will be the iso surface, and
    //! the iso value can be specified. For the boundaries (the walls), it can be
    //! specified whether to close or open with \p bndClose (default: close all).
    //! Another boundary flag \p bndConnectivity can be used for specifying
    //! topological connectivity of the boundary meshes (default: disconnect all).
    //!
    //! \param[in]  grid            The grid.
    //! \param[in]  gridSize        The grid size.
    //! \param[in]  origin          The origin.
    //! \param[out] mesh            The output triangle mesh.
    //! \param[in]  isoValue        The iso-surface value.
    //! \param[in]  bndClose        The boundary open flag.
    //! \param[in]  bndConnectivity The boundary connectivity flag.
    //!
    void marchingCubes(const ConstArrayAccessor3<double>& grid,
        const Vector3D& gridSize, const Vector3D& origin,
        TriangleMesh3* mesh, const FaceCenteredGrid3Ptr& uvInfo, double isoValue = 0,
        int bndClose = kDirectionAll,
        int bndConnectivity = kDirectionNone);

    typedef size_t MarchingCubeVertexHashKey;
    typedef size_t MarchingCubeVertexId;
    typedef std::unordered_map<MarchingCubeVertexHashKey, MarchingCubeVertexId>
        MarchingCubeVertexMap;

    inline bool queryVertexId(const MarchingCubeVertexMap& vertexMap,
        MarchingCubeVertexHashKey vKey,
        MarchingCubeVertexId* vId) {
        auto vItr = vertexMap.find(vKey);
        if (vItr != vertexMap.end()) {
            *vId = vItr->second;
            return true;
        }
        else {
            return false;
        }
    }

    inline Vector3D grad(const ConstArrayAccessor3<double>& grid, ssize_t i,
        ssize_t j, ssize_t k, const Vector3D& invGridSize) {
        Vector3D ret;
        ssize_t ip = i + 1;
        ssize_t im = i - 1;
        ssize_t jp = j + 1;
        ssize_t jm = j - 1;
        ssize_t kp = k + 1;
        ssize_t km = k - 1;
        Size3 dim = grid.size();
        ssize_t dimx = static_cast<ssize_t>(dim.x);
        ssize_t dimy = static_cast<ssize_t>(dim.y);
        ssize_t dimz = static_cast<ssize_t>(dim.z);
        if (i > dimx - 2) {
            ip = i;
        }
        else if (i == 0) {
            im = 0;
        }
        if (j > dimy - 2) {
            jp = j;
        }
        else if (j == 0) {
            jm = 0;
        }
        if (k > dimz - 2) {
            kp = k;
        }
        else if (k == 0) {
            km = 0;
        }
        ret.x = 0.5f * invGridSize.x * (grid(ip, j, k) - grid(im, j, k));
        ret.y = 0.5f * invGridSize.y * (grid(i, jp, k) - grid(i, jm, k));
        ret.z = 0.5f * invGridSize.z * (grid(i, j, kp) - grid(i, j, km));
        return ret;
    }

    inline Vector3D safeNormalize(const Vector3D& n) {
        if (n.lengthSquared() > 0.0) {
            return n.normalized();
        }
        else {
            return n;
        }
    }

    // To compute unique edge ID, map vertices+edges into
    // doubled virtual vertex indices.
    //
    // v  edge   v
    // |----*----|    -->    |-----|-----|
    // i        i+1         2i   2i+1  2i+2
    //
    inline size_t globalEdgeId(size_t i, size_t j, size_t k, const Size3& dim,
        size_t localEdgeId) {
        // See edgeConnection in marching_cubes_table.h for the edge ordering.
        static const int edgeOffset3D[12][3] = {
            {1, 0, 0}, {2, 0, 1}, {1, 0, 2}, {0, 0, 1}, {1, 2, 0}, {2, 2, 1},
            {1, 2, 2}, {0, 2, 1}, {0, 1, 0}, {2, 1, 0}, {2, 1, 2}, {0, 1, 2} };

        return ((2 * k + edgeOffset3D[localEdgeId][2]) * 2 * dim.y +
            (2 * j + edgeOffset3D[localEdgeId][1])) *
            2 * dim.x +
            (2 * i + edgeOffset3D[localEdgeId][0]);
    }

    // To compute unique edge ID, map vertices+edges into
    // doubled virtual vertex indices.
    //
    // v  edge   v
    // |----*----|    -->    |-----|-----|
    // i        i+1         2i   2i+1  2i+2
    //
    inline size_t globalVertexId(size_t i, size_t j, size_t k, const Size3& dim,
        size_t localVertexId) {
        // See edgeConnection in marching_cubes_table.h for the edge ordering.
        static const int vertexOffset3D[8][3] = { {0, 0, 0}, {2, 0, 0}, {2, 0, 2},
                                                 {0, 0, 2}, {0, 2, 0}, {2, 2, 0},
                                                 {2, 2, 2}, {0, 2, 2} };

        return ((2 * k + vertexOffset3D[localVertexId][2]) * 2 * dim.y +
            (2 * j + vertexOffset3D[localVertexId][1])) *
            2 * dim.x +
            (2 * i + vertexOffset3D[localVertexId][0]);
    }

    static void singleSquare(const std::array<double, 4>& data,
        const std::array<size_t, 8>& vertAndEdgeIds,
        const Vector3D& normal,
        const std::array<Vector3D, 4>& corners,
        MarchingCubeVertexMap* vertexMap, TriangleMesh3* mesh,
        const FaceCenteredGrid3Ptr& uvInfo, double isoValue) {
        int itrVertex, itrEdge, itrTri;
        int idxFlags = 0, idxEdgeFlags = 0;
        int idxVertexOfTheEdge[2];

        double phi0, phi1, alpha;
        Vector3D pos, pos0, pos1;
        Vector3D e[4];

        // Which vertices are inside? If i-th vertex is inside, mark '1' at i-th
        // bit. of 'idxFlags'.
        for (itrVertex = 0; itrVertex < 4; itrVertex++) {
            if (data[itrVertex] <= isoValue) {
                idxFlags |= 1 << itrVertex;
            }
        }

        // If the rect is entirely outside of the surface,
        // there is no job to be done in this marching-cube cell.
        if (idxFlags == 0) {
            return;
        }

        // If there are vertices which is inside the surface...
        // Which edges intersect the surface?
        // If i-th edge intersects the surface, mark '1' at i-th bit of
        // 'idxEdgeFlags'
        idxEdgeFlags = squareEdgeFlags[idxFlags];

        // Find the point of intersection of the surface with each edge
        for (itrEdge = 0; itrEdge < 4; itrEdge++) {
            // If there is an intersection on this edge
            if (idxEdgeFlags & (1 << itrEdge)) {
                idxVertexOfTheEdge[0] = edgeConnection2D[itrEdge][0];
                idxVertexOfTheEdge[1] = edgeConnection2D[itrEdge][1];

                // Find the phi = 0 position by iteration
                pos0 = corners[idxVertexOfTheEdge[0]];
                pos1 = corners[idxVertexOfTheEdge[1]];

                phi0 = data[idxVertexOfTheEdge[0]] - isoValue;
                phi1 = data[idxVertexOfTheEdge[1]] - isoValue;

                // I think it needs perturbation a little bit.
                if (std::fabs(phi0) + std::fabs(phi1) > 1e-12) {
                    alpha = std::fabs(phi0) / (std::fabs(phi0) + std::fabs(phi1));
                }
                else {
                    alpha = 0.5f;
                }

                if (alpha < 0.000001f) {
                    alpha = 0.000001f;
                }
                if (alpha > 0.999999f) {
                    alpha = 0.999999f;
                }

                pos = ((1.f - alpha) * pos0 + alpha * pos1);

                // What is the position of this vertex of the edge?
                e[itrEdge] = pos;
            }
        }

        // Make triangular patches.
        for (itrTri = 0; itrTri < 4; ++itrTri) {
            // If there isn't any triangle to be built, escape this loop.
            if (triangleConnectionTable2D[idxFlags][3 * itrTri] < 0) {
                break;
            }

            Point3UI face;

            for (int j = 0; j < 3; ++j) {
                int idxVertex = triangleConnectionTable2D[idxFlags][3 * itrTri + j];

                MarchingCubeVertexHashKey vKey = vertAndEdgeIds[idxVertex];
                MarchingCubeVertexId vId;
                if (queryVertexId(*vertexMap, vKey, &vId)) {
                    face[j] = vId;
                }
                else {
                    // if vertex does not exist...
                    face[j] = mesh->numberOfPoints();
                    mesh->addNormal(normal);
                    if (idxVertex < 4) {
                        mesh->addPoint(corners[idxVertex]);
                    }
                    else {
                        mesh->addPoint(e[idxVertex - 4]);
                    }

                    float dataUv = 0.0f;

                    if (uvInfo) {

                        Size3 position = Size3(std::floor(corners[idxVertexOfTheEdge[0]].x * uvInfo->uSize().x), std::floor(corners[idxVertexOfTheEdge[0]].y * uvInfo->uSize().y), std::floor(corners[idxVertexOfTheEdge[0]].z * uvInfo->uSize().z));



                        dataUv =
                            pow(pow(uvInfo->u(position.x, position.y, position.z), 2) +
                                pow(uvInfo->v(position.x, position.y, position.z), 2) +
                                pow(uvInfo->w(position.x, position.y, position.z), 2), 0.5) / 3;



                        if (dataUv > 1.0f) {

                            dataUv = 1.0f;
                        }
                    }

                



                    mesh->addUv(Vector2D(0.0, dataUv));
                    vertexMap->insert(std::make_pair(vKey, face[j]));
                }
            }

            mesh->addPointUvNormalTriangle(face, face, face);
        }
    }

    static void singleCube(const std::array<double, 8>& data,
        const std::array<size_t, 12>& edgeIds,
        const std::array<Vector3D, 8>& normals,
        const BoundingBox3D& bound,
        MarchingCubeVertexMap* vertexMap, TriangleMesh3* mesh,
        const FaceCenteredGrid3Ptr& uvInfo, double isoValue) {
        int itrVertex, itrEdge, itrTri;
        int idxFlagSize = 0, idxEdgeFlags = 0;
        int idxVertexOfTheEdge[2];

        Vector3D pos, pos0, pos1, normal, normal0, normal1;
        double phi0, phi1;
        double alpha;
        Vector3D e[12], n[12];

        // Which vertices are inside? If i-th vertex is inside, mark '1' at i-th
        // bit. of 'idxFlagSize'.
        for (itrVertex = 0; itrVertex < 8; itrVertex++) {
            if (data[itrVertex] <= isoValue) {
                idxFlagSize |= 1 << itrVertex;
            }
        }

        // If the cube is entirely inside or outside of the surface, there is no job
        // to be done in t1his marching-cube cell.
        if (idxFlagSize == 0 || idxFlagSize == 255) {
            return;
        }

        // If there are vertices which is inside the surface...
        // Which edges intersect the surface? If i-th edge intersects the surface,
        // mark '1' at i-th bit of 'itrEdgeFlags'
        idxEdgeFlags = cubeEdgeFlags[idxFlagSize];

        // Find the point of intersection of the surface with each edge
        for (itrEdge = 0; itrEdge < 12; itrEdge++) {
            // If there is an intersection on this edge
            if (idxEdgeFlags & (1 << itrEdge)) {
                idxVertexOfTheEdge[0] = edgeConnection[itrEdge][0];
                idxVertexOfTheEdge[1] = edgeConnection[itrEdge][1];

                // cube vertex ordering to x-major ordering
                static int indexMap[8] = { 0, 1, 5, 4, 2, 3, 7, 6 };

                // Find the phi = 0 position
                pos0 = bound.corner(indexMap[idxVertexOfTheEdge[0]]);
                pos1 = bound.corner(indexMap[idxVertexOfTheEdge[1]]);

                normal0 = normals[idxVertexOfTheEdge[0]];
                normal1 = normals[idxVertexOfTheEdge[1]];

                phi0 = data[idxVertexOfTheEdge[0]] - isoValue;
                phi1 = data[idxVertexOfTheEdge[1]] - isoValue;

                alpha = distanceToZeroLevelSet(phi0, phi1);

                if (alpha < 0.000001) {
                    alpha = 0.000001;
                }
                if (alpha > 0.999999) {
                    alpha = 0.999999;
                }

                pos = (1.0 - alpha) * pos0 + alpha * pos1;
                normal = (1.0 - alpha) * normal0 + alpha * normal1;

                e[itrEdge] = pos;
                n[itrEdge] = normal;
            }
        }

        // Make triangles
        for (itrTri = 0; itrTri < 5; ++itrTri) {
            // If there isn't any triangle to be made, escape this loop.
            if (triangleConnectionTable3D[idxFlagSize][3 * itrTri] < 0) {
                break;
            }

            Point3UI face;

            for (int j = 0; j < 3; j++) {
                int k = 3 * itrTri + j;
                MarchingCubeVertexHashKey vKey =
                    edgeIds[triangleConnectionTable3D[idxFlagSize][k]];
                MarchingCubeVertexId vId;
                if (queryVertexId(*vertexMap, vKey, &vId)) {
                    face[j] = vId;
                }
                else {
                    // If vertex does not exist from the map
                    face[j] = mesh->numberOfPoints();
                    mesh->addNormal(safeNormalize(
                        n[triangleConnectionTable3D[idxFlagSize][k]]));
                    mesh->addPoint(e[triangleConnectionTable3D[idxFlagSize][k]]);

                    float dataUv = 0.0f;
                    
                    
                    if (uvInfo) {
                        //std::cout << dataAccessor.size().x << std::endl;
                        //std::cout << (int)round(bound.lowerCorner.x * dataAccessor.size().x) << " " << (int)round(bound.lowerCorner.y * dataAccessor.size().y) << " " << (int)round(bound.lowerCorner.z * dataAccessor.size().z) << std::endl;
                        Size3 position = Size3((int)round(bound.lowerCorner.x * uvInfo->uSize().x), (int)round(bound.lowerCorner.y * uvInfo->uSize().y), (int)round(bound.lowerCorner.z * uvInfo->uSize().z));


                        dataUv =
                            pow(pow(uvInfo->u(position.x, position.y, position.z), 2) +
                                pow(uvInfo->v(position.x, position.y, position.z), 2) +
                                pow(uvInfo->w(position.x, position.y, position.z), 2), 0.5) / 3;



                        if (dataUv > 1.0f) {

                            dataUv = 1.0f;
                        }

                    }


               
                    
                    mesh->addUv(Vector2D(0.0, dataUv));
                    vertexMap->insert(std::make_pair(vKey, face[j]));
                }
            }
            mesh->addPointUvNormalTriangle(face, face, face);
        }
    }

    void marchingCubes(const ConstArrayAccessor3<double>& grid,
        const Vector3D& gridSize, const Vector3D& origin,
        TriangleMesh3* mesh, const FaceCenteredGrid3Ptr& uvInfo, double isoValue, int bndClose,
        int bndConnectivity) {
        MarchingCubeVertexMap vertexMap;

        const Size3 dim = grid.size();
        const Vector3D invGridSize = 1.0 / gridSize;

        auto pos = [origin, gridSize](ssize_t i, ssize_t j, ssize_t k) {
            return origin + gridSize * Vector3D({ i, j, k });
        };

        ssize_t dimx = static_cast<ssize_t>(dim.x);
        ssize_t dimy = static_cast<ssize_t>(dim.y);
        ssize_t dimz = static_cast<ssize_t>(dim.z);

        for (ssize_t k = 0; k < dimz - 1; ++k) {
            for (ssize_t j = 0; j < dimy - 1; ++j) {
                for (ssize_t i = 0; i < dimx - 1; ++i) {
                    std::array<double, 8> data;
                    std::array<size_t, 12> edgeIds;
                    std::array<Vector3D, 8> normals;
                    BoundingBox3D bound;

                    data[0] = grid(i, j, k);
                    data[1] = grid(i + 1, j, k);
                    data[4] = grid(i, j + 1, k);
                    data[5] = grid(i + 1, j + 1, k);
                    data[3] = grid(i, j, k + 1);
                    data[2] = grid(i + 1, j, k + 1);
                    data[7] = grid(i, j + 1, k + 1);
                    data[6] = grid(i + 1, j + 1, k + 1);

                    normals[0] = grad(grid, i, j, k, invGridSize);
                    normals[1] = grad(grid, i + 1, j, k, invGridSize);
                    normals[4] = grad(grid, i, j + 1, k, invGridSize);
                    normals[5] = grad(grid, i + 1, j + 1, k, invGridSize);
                    normals[3] = grad(grid, i, j, k + 1, invGridSize);
                    normals[2] = grad(grid, i + 1, j, k + 1, invGridSize);
                    normals[7] = grad(grid, i, j + 1, k + 1, invGridSize);
                    normals[6] = grad(grid, i + 1, j + 1, k + 1, invGridSize);

                    for (int e = 0; e < 12; e++) {
                        edgeIds[e] = globalEdgeId(i, j, k, dim, e);
                    }

                    bound.lowerCorner = pos(i, j, k);
                    bound.upperCorner = pos(i + 1, j + 1, k + 1);

                    singleCube(data, edgeIds, normals, bound, &vertexMap, mesh,
                        uvInfo, isoValue);
                }  // i
            }      // j
        }          // k

        // Construct boundaries parallel to x-y plane
        if (bndClose & (kDirectionBack | kDirectionFront)) {
            MarchingCubeVertexMap vertexMapBack;
            MarchingCubeVertexMap vertexMapFront;

            for (ssize_t j = 0; j < dimy - 1; ++j) {
                for (ssize_t i = 0; i < dimx - 1; ++i) {
                    ssize_t k = 0;

                    std::array<double, 4> data;
                    std::array<size_t, 8> vertexAndEdgeIds;
                    std::array<Vector3D, 4> corners;
                    Vector3D normal;
                    BoundingBox2D bound;

                    data[0] = grid(i + 1, j, k);
                    data[1] = grid(i, j, k);
                    data[2] = grid(i, j + 1, k);
                    data[3] = grid(i + 1, j + 1, k);

                    if (bndClose & kDirectionBack) {
                        normal = Vector3D(0, 0, -1);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 1);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 0);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 4);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 5);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 0);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 8);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 4);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 9);

                        corners[0] = pos(i + 1, j, k);
                        corners[1] = pos(i, j, k);
                        corners[2] = pos(i, j + 1, k);
                        corners[3] = pos(i + 1, j + 1, k);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionBack)
                            ? &vertexMap
                            : &vertexMapBack,
                            mesh, uvInfo, isoValue);
                    }

                    k = dimz - 2;
                    data[0] = grid(i, j, k + 1);
                    data[1] = grid(i + 1, j, k + 1);
                    data[2] = grid(i + 1, j + 1, k + 1);
                    data[3] = grid(i, j + 1, k + 1);

                    if (bndClose & kDirectionFront) {
                        normal = Vector3D(0, 0, 1);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 3);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 2);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 6);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 7);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 2);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 10);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 6);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 11);

                        corners[0] = pos(i, j, k + 1);
                        corners[1] = pos(i + 1, j, k + 1);
                        corners[2] = pos(i + 1, j + 1, k + 1);
                        corners[3] = pos(i, j + 1, k + 1);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionFront)
                            ? &vertexMap
                            : &vertexMapFront,
                            mesh, uvInfo, isoValue);
                    }
                }  // i
            }      // j
        }

        // Construct boundaries parallel to y-z plane
        if (bndClose & (kDirectionLeft | kDirectionRight)) {
            MarchingCubeVertexMap vertexMapLeft;
            MarchingCubeVertexMap vertexMapRight;

            for (ssize_t k = 0; k < dimz - 1; ++k) {
                for (ssize_t j = 0; j < dimy - 1; ++j) {
                    ssize_t i = 0;

                    std::array<double, 4> data;
                    std::array<size_t, 8> vertexAndEdgeIds;
                    std::array<Vector3D, 4> corners;
                    Vector3D normal;
                    BoundingBox2D bound;

                    data[0] = grid(i, j, k);
                    data[1] = grid(i, j, k + 1);
                    data[2] = grid(i, j + 1, k + 1);
                    data[3] = grid(i, j + 1, k);

                    if (bndClose & kDirectionLeft) {
                        normal = Vector3D(-1, 0, 0);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 0);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 3);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 7);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 4);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 3);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 11);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 7);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 8);

                        corners[0] = pos(i, j, k);
                        corners[1] = pos(i, j, k + 1);
                        corners[2] = pos(i, j + 1, k + 1);
                        corners[3] = pos(i, j + 1, k);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionLeft)
                            ? &vertexMap
                            : &vertexMapLeft,
                            mesh, uvInfo,  isoValue);
                    }

                    i = dimx - 2;
                    data[0] = grid(i + 1, j, k + 1);
                    data[1] = grid(i + 1, j, k);
                    data[2] = grid(i + 1, j + 1, k);
                    data[3] = grid(i + 1, j + 1, k + 1);

                    if (bndClose & kDirectionRight) {
                        normal = Vector3D(1, 0, 0);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 2);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 1);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 5);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 6);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 1);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 9);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 5);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 10);

                        corners[0] = pos(i + 1, j, k + 1);
                        corners[1] = pos(i + 1, j, k);
                        corners[2] = pos(i + 1, j + 1, k);
                        corners[3] = pos(i + 1, j + 1, k + 1);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionRight)
                            ? &vertexMap
                            : &vertexMapRight,
                            mesh, uvInfo, isoValue);
                    }
                }  // j
            }      // k
        }

        // Construct boundaries parallel to x-z plane
        if (bndClose & (kDirectionDown | kDirectionUp)) {
            MarchingCubeVertexMap vertexMapDown;
            MarchingCubeVertexMap vertexMapUp;

            for (ssize_t k = 0; k < dimz - 1; ++k) {
                for (ssize_t i = 0; i < dimx - 1; ++i) {
                    ssize_t j = 0;

                    std::array<double, 4> data;
                    std::array<size_t, 8> vertexAndEdgeIds;
                    std::array<Vector3D, 4> corners;
                    Vector3D normal;
                    BoundingBox2D bound;

                    data[0] = grid(i, j, k);
                    data[1] = grid(i + 1, j, k);
                    data[2] = grid(i + 1, j, k + 1);
                    data[3] = grid(i, j, k + 1);

                    if (bndClose & kDirectionDown) {
                        normal = Vector3D(0, -1, 0);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 0);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 1);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 2);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 3);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 0);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 1);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 2);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 3);

                        corners[0] = pos(i, j, k);
                        corners[1] = pos(i + 1, j, k);
                        corners[2] = pos(i + 1, j, k + 1);
                        corners[3] = pos(i, j, k + 1);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionDown)
                            ? &vertexMap
                            : &vertexMapDown,
                            mesh, uvInfo, isoValue);
                    }

                    j = dimy - 2;
                    data[0] = grid(i + 1, j + 1, k);
                    data[1] = grid(i, j + 1, k);
                    data[2] = grid(i, j + 1, k + 1);
                    data[3] = grid(i + 1, j + 1, k + 1);

                    if (bndClose & kDirectionUp) {
                        normal = Vector3D(0, 1, 0);

                        vertexAndEdgeIds[0] = globalVertexId(i, j, k, dim, 5);
                        vertexAndEdgeIds[1] = globalVertexId(i, j, k, dim, 4);
                        vertexAndEdgeIds[2] = globalVertexId(i, j, k, dim, 7);
                        vertexAndEdgeIds[3] = globalVertexId(i, j, k, dim, 6);
                        vertexAndEdgeIds[4] = globalEdgeId(i, j, k, dim, 4);
                        vertexAndEdgeIds[5] = globalEdgeId(i, j, k, dim, 7);
                        vertexAndEdgeIds[6] = globalEdgeId(i, j, k, dim, 6);
                        vertexAndEdgeIds[7] = globalEdgeId(i, j, k, dim, 5);

                        corners[0] = pos(i + 1, j + 1, k);
                        corners[1] = pos(i, j + 1, k);
                        corners[2] = pos(i, j + 1, k + 1);
                        corners[3] = pos(i + 1, j + 1, k + 1);

                        singleSquare(data, vertexAndEdgeIds, normal, corners,
                            (bndConnectivity & kDirectionUp)
                            ? &vertexMap
                            : &vertexMapUp,
                            mesh, uvInfo, isoValue);
                    }
                }  // i
            }      // k
        }
    }



}  // namespace jet

#endif  // INCLUDE_JET_MARCHING_CUBES_H_