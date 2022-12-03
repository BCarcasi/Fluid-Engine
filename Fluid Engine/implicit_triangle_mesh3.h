#ifndef INCLUDE_JET_IMPLICIT_TRIANGLE_MESH3_H_
#define INCLUDE_JET_IMPLICIT_TRIANGLE_MESH3_H_

#include "custom_implicit_surface3.h"
#include "implicit_surface3.h"
#include "triangle_mesh3.h"
#include "vertex_centered_scalar_grid3.h"

#include "pch.h"

#include "triangle_mesh_to_sdf.h"

namespace jet {

    //!
    //! \brief  TriangleMesh3 to ImplicitSurface3 converter.
    //!
    //! This class builds signed-distance field for given TriangleMesh3 instance so
    //! that it can be used as an ImplicitSurface3 instance. The mesh is discretized
    //! into a regular grid and the signed-distance is measured at each grid point.
    //! Thus, there is a sampling error and its magnitude depends on the grid
    //! resolution.
    //!
    class ImplicitTriangleMesh3 final : public ImplicitSurface3 {
    public:
        class Builder;

        //! Constructs an ImplicitSurface3 with mesh and other grid parameters.
        ImplicitTriangleMesh3(const TriangleMesh3Ptr& mesh, size_t resolutionX = 32,
            double margin = 0.2,
            const Transform3& transform = Transform3(),
            bool isNormalFlipped = false);

        virtual ~ImplicitTriangleMesh3();

        //! Returns builder fox ImplicitTriangleMesh3.
        static Builder builder();

        //! Returns grid data.
        const VertexCenteredScalarGrid3Ptr& grid() const;

    private:
        TriangleMesh3Ptr _mesh;
        VertexCenteredScalarGrid3Ptr _grid;
        CustomImplicitSurface3Ptr _customImplicitSurface;

        Vector3D closestPointLocal(const Vector3D& otherPoint) const override;

        double closestDistanceLocal(const Vector3D& otherPoint) const override;

        bool intersectsLocal(const Ray3D& ray) const override;

        BoundingBox3D boundingBoxLocal() const override;

        Vector3D closestNormalLocal(const Vector3D& otherPoint) const override;

        double signedDistanceLocal(const Vector3D& otherPoint) const override;

        SurfaceRayIntersection3 closestIntersectionLocal(
            const Ray3D& ray) const override;
    };

    //! Shared pointer for the ImplicitTriangleMesh3 type.
    typedef std::shared_ptr<ImplicitTriangleMesh3> ImplicitTriangleMesh3Ptr;

    //!
    //! \brief Front-end to create ImplicitTriangleMesh3 objects step by step.
    //!
    class ImplicitTriangleMesh3::Builder final
        : public SurfaceBuilderBase3<ImplicitTriangleMesh3::Builder> {
    public:
        //! Returns builder with triangle mesh.
        Builder& withTriangleMesh(const TriangleMesh3Ptr& mesh);

        //! Returns builder with resolution in x axis.
        Builder& withResolutionX(size_t resolutionX);

        //! Returns builder with margin around the mesh.
        Builder& withMargin(double margin);

        //! Builds ImplicitTriangleMesh3.
        ImplicitTriangleMesh3 build() const;

        //! Builds shared pointer of ImplicitTriangleMesh3 instance.
        ImplicitTriangleMesh3Ptr makeShared() const;

    private:
        TriangleMesh3Ptr _mesh;
        size_t _resolutionX = 32;
        double _margin = 0.2;
    };

    ImplicitTriangleMesh3::ImplicitTriangleMesh3(const TriangleMesh3Ptr& mesh,
        size_t resolutionX, double margin,
        const Transform3& transform,
        bool isNormalFlipped)
        : ImplicitSurface3(transform, isNormalFlipped), _mesh(mesh) {
        if (mesh->numberOfTriangles() > 0 && mesh->numberOfPoints() > 0) {
            BoundingBox3D box = _mesh->boundingBox();
            Vector3D scale(box.width(), box.height(), box.depth());
            box.lowerCorner -= margin * scale;
            box.upperCorner += margin * scale;
            size_t resolutionY = static_cast<size_t>(
                std::ceil(resolutionX * box.height() / box.width()));
            size_t resolutionZ = static_cast<size_t>(
                std::ceil(resolutionX * box.depth() / box.width()));

            double dx = box.width() / resolutionX;

            _grid = std::make_shared<VertexCenteredScalarGrid3>();
            _grid->resize(resolutionX, resolutionY, resolutionZ, dx, dx, dx,
                box.lowerCorner.x, box.lowerCorner.y, box.lowerCorner.z);

            triangleMeshToSdf(*_mesh, _grid.get());

            _customImplicitSurface =
                CustomImplicitSurface3::builder()
                .withSignedDistanceFunction([&](const Vector3D& pt) -> double {
                return _grid->sample(pt);
                    })
                .withDomain(_grid->boundingBox())
                        .withResolution(dx)
                        .makeShared();
        }
        else {
            // Empty mesh -- return big/uniform number
            _customImplicitSurface =
                CustomImplicitSurface3::builder()
                .withSignedDistanceFunction(
                    [&](const Vector3D&) -> double { return kMaxD; })
                .makeShared();
        }
    }

    ImplicitTriangleMesh3::~ImplicitTriangleMesh3() {}

    Vector3D ImplicitTriangleMesh3::closestPointLocal(
        const Vector3D& otherPoint) const {
        return _customImplicitSurface->closestPoint(otherPoint);
    }

    double ImplicitTriangleMesh3::closestDistanceLocal(
        const Vector3D& otherPoint) const {
        return _customImplicitSurface->closestDistance(otherPoint);
    }

    bool ImplicitTriangleMesh3::intersectsLocal(const Ray3D& ray) const {
        return _customImplicitSurface->intersects(ray);
    }

    BoundingBox3D ImplicitTriangleMesh3::boundingBoxLocal() const {
        return _mesh->boundingBox();
    }

    Vector3D ImplicitTriangleMesh3::closestNormalLocal(
        const Vector3D& otherPoint) const {
        return _customImplicitSurface->closestNormal(otherPoint);
    }

    double ImplicitTriangleMesh3::signedDistanceLocal(
        const Vector3D& otherPoint) const {
        return _customImplicitSurface->signedDistance(otherPoint);
    }

    SurfaceRayIntersection3 ImplicitTriangleMesh3::closestIntersectionLocal(
        const Ray3D& ray) const {
        return _customImplicitSurface->closestIntersection(ray);
    }

    ImplicitTriangleMesh3::Builder ImplicitTriangleMesh3::builder() {
        return ImplicitTriangleMesh3::Builder();
    }

    const VertexCenteredScalarGrid3Ptr& ImplicitTriangleMesh3::grid() const {
        return _grid;
    }

    ImplicitTriangleMesh3::Builder&
        ImplicitTriangleMesh3::Builder::withTriangleMesh(const TriangleMesh3Ptr& mesh) {
        _mesh = mesh;
        return *this;
    }

    ImplicitTriangleMesh3::Builder& ImplicitTriangleMesh3::Builder::withResolutionX(
        size_t resolutionX) {
        _resolutionX = resolutionX;
        return *this;
    }

    ImplicitTriangleMesh3::Builder& ImplicitTriangleMesh3::Builder::withMargin(
        double margin) {
        _margin = margin;
        return *this;
    }

    ImplicitTriangleMesh3 ImplicitTriangleMesh3::Builder::build() const {
        return ImplicitTriangleMesh3(_mesh, _resolutionX, _margin, _transform,
            _isNormalFlipped);
    }

    ImplicitTriangleMesh3Ptr ImplicitTriangleMesh3::Builder::makeShared() const {
        return std::shared_ptr<ImplicitTriangleMesh3>(
            new ImplicitTriangleMesh3(_mesh, _resolutionX, _margin, _transform,
                _isNormalFlipped),
            [](ImplicitTriangleMesh3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_IMPLICIT_TRIANGLE_MESH3_H_