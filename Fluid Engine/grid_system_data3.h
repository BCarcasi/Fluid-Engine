#ifndef INCLUDE_JET_GRID_SYSTEM_DATA3_H_
#define INCLUDE_JET_GRID_SYSTEM_DATA3_H_

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "face_centered_grid3.h"
#include "scalar_grid3.h"
#include "serialization.h"
#include <memory>
#include <vector>



#include "pch.h"

#include "factory.h"
#include "fbs_helpers.h"
#include "grid_system_data3_generated.h"

#include "grid_system_data3.h"

#include "flatbuffers.h"

#include <algorithm>
#include <vector>

namespace jet {

    //!
    //! \brief      3-D grid system data.
    //!
    //! This class is the key data structure for storing grid system data. To
    //! represent a grid system for fluid simulation, velocity field is defined as a
    //! face-centered (MAC) grid by default. It can also have additional scalar or
    //! vector attributes by adding extra data layer.
    //!
    class GridSystemData3 : public Serializable {
    public:
        //! Constructs empty grid system.
        GridSystemData3();

        //!
        //! \brief      Constructs a grid system with given resolution, grid spacing
        //!             and origin.
        //!
        //! This constructor builds the entire grid layers within the system. Note,
        //! the resolution is the grid resolution, not the data size of each grid.
        //! Depending on the layout of the grid, the data point may lie on different
        //! part of the grid (vertex, cell-center, or face-center), thus can have
        //! different array size internally. The resolution of the grid means the
        //! grid cell resolution.
        //!
        //! \param[in]  resolution  The resolution.
        //! \param[in]  gridSpacing The grid spacing.
        //! \param[in]  origin      The origin.
        //!
        GridSystemData3(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& origin);

        //! Copy constructor.
        GridSystemData3(const GridSystemData3& other);

        //! Destructor.
        virtual ~GridSystemData3();

        //!
        //! \brief      Resizes the whole system with given resolution, grid
        //!             spacing, and origin.
        //!
        //! This function resizes the entire grid layers within the system. Note,
        //! the resolution is the grid resolution, not the data size of each grid.
        //! Depending on the layout of the grid, the data point may lie on different
        //! part of the grid (vertex, cell-center, or face-center), thus can have
        //! different array size internally. The resolution of the grid means the
        //! grid cell resolution.
        //!
        //! \param[in]  resolution  The resolution.
        //! \param[in]  gridSpacing The grid spacing.
        //! \param[in]  origin      The origin.
        //!
        void resize(
            const Size3& resolution,
            const Vector3D& gridSpacing,
            const Vector3D& origin);

        //!
        //! \brief      Returns the resolution of the grid.
        //!
        //! This function resizes the entire grid layers within the system. Note,
        //! the resolution is the grid resolution, not the data size of each grid.
        //! Depending on the layout of the grid, the data point may lie on different
        //! part of the grid (vertex, cell-center, or face-center), thus can have
        //! different array size internally. The resolution of the grid means the
        //! grid cell resolution.
        //!
        //! \return     Grid cell resolution.
        //!
        Size3 resolution() const;

        //! Return the grid spacing.
        Vector3D gridSpacing() const;

        //! Returns the origin of the grid.
        Vector3D origin() const;

        //! Returns the bounding box of the grid.
        BoundingBox3D boundingBox() const;

        //!
        //! \brief      Adds a non-advectable scalar data grid by passing its
        //!     builder and initial value.
        //!
        //! This function adds a new scalar data grid. This layer is not advectable,
        //! meaning that during the computation of fluid flow, this layer won't
        //! follow the flow. For the future access of this layer, its index is
        //! returned.
        //!
        //! \param[in]  builder    The grid builder.
        //! \param[in]  initialVal The initial value.
        //!
        //! \return     Index of the data.
        //!
        size_t addScalarData(
            const ScalarGridBuilder3Ptr& builder,
            double initialVal = 0.0);

        //!
        //! \brief      Adds a non-advectable vector data grid by passing its
        //!     builder and initial value.
        //!
        //! This function adds a new vector data grid. This layer is not advectable,
        //! meaning that during the computation of fluid flow, this layer won't
        //! follow the flow. For the future access of this layer, its index is
        //! returned.
        //!
        //! \param[in]  builder    The grid builder.
        //! \param[in]  initialVal The initial value.
        //!
        //! \return     Index of the data.
        //!
        size_t addVectorData(
            const VectorGridBuilder3Ptr& builder,
            const Vector3D& initialVal = Vector3D());

        //!
        //! \brief      Adds an advectable scalar data grid by passing its builder
        //!     and initial value.
        //!
        //! This function adds a new scalar data grid. This layer is advectable,
        //! meaning that during the computation of fluid flow, this layer will
        //! follow the flow. For the future access of this layer, its index is
        //! returned.
        //!
        //! \param[in]  builder    The grid builder.
        //! \param[in]  initialVal The initial value.
        //!
        //! \return     Index of the data.
        //!
        size_t addAdvectableScalarData(
            const ScalarGridBuilder3Ptr& builder,
            double initialVal = 0.0);

        //!
        //! \brief      Adds an advectable vector data grid by passing its builder
        //!     and initial value.
        //!
        //! This function adds a new vector data grid. This layer is advectable,
        //! meaning that during the computation of fluid flow, this layer will
        //! follow the flow. For the future access of this layer, its index is
        //! returned.
        //!
        //! \param[in]  builder    The grid builder.
        //! \param[in]  initialVal The initial value.
        //!
        //! \return     Index of the data.
        //!
        size_t addAdvectableVectorData(
            const VectorGridBuilder3Ptr& builder,
            const Vector3D& initialVal = Vector3D());

        //!
        //! \brief      Returns the velocity field.
        //!
        //! This class has velocify field by default, and it is part of the
        //! advectable vector data list.
        //!
        //! \return     Pointer to the velocity field.
        //!
        const FaceCenteredGrid3Ptr& velocity() const;

        //!
        //! \brief      Returns the index of the velocity field.
        //!
        //! This class has velocify field by default, and it is part of the
        //! advectable vector data list. This function returns the index of the
        //! velocity field from the list.
        //!
        //! \return     Index of the velocity field.
        //!
        size_t velocityIndex() const;

        //! Returns the non-advectable scalar data at given index.
        const ScalarGrid3Ptr& scalarDataAt(size_t idx) const;

        //! Returns the non-advectable vector data at given index.
        const VectorGrid3Ptr& vectorDataAt(size_t idx) const;

        //! Returns the advectable scalar data at given index.
        const ScalarGrid3Ptr& advectableScalarDataAt(size_t idx) const;

        //! Returns the advectable vector data at given index.
        const VectorGrid3Ptr& advectableVectorDataAt(size_t idx) const;

        //! Returns the number of non-advectable scalar data.
        size_t numberOfScalarData() const;

        //! Returns the number of non-advectable vector data.
        size_t numberOfVectorData() const;

        //! Returns the number of advectable scalar data.
        size_t numberOfAdvectableScalarData() const;

        //! Returns the number of advectable vector data.
        size_t numberOfAdvectableVectorData() const;

        //! Serialize the data to the given buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Serialize the data from the given buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

    private:
        Size3 _resolution;
        Vector3D _gridSpacing;
        Vector3D _origin;

        FaceCenteredGrid3Ptr _velocity;
        size_t _velocityIdx;
        std::vector<ScalarGrid3Ptr> _scalarDataList;
        std::vector<VectorGrid3Ptr> _vectorDataList;
        std::vector<ScalarGrid3Ptr> _advectableScalarDataList;
        std::vector<VectorGrid3Ptr> _advectableVectorDataList;
    };

    //! Shared pointer type of GridSystemData3.
    typedef std::shared_ptr<GridSystemData3> GridSystemData3Ptr;

    GridSystemData3::GridSystemData3()
        : GridSystemData3({ 0, 0, 0 }, { 1, 1, 1 }, { 0, 0, 0 }) {
    }

    GridSystemData3::GridSystemData3(
        const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin) {
        _velocity = std::make_shared<FaceCenteredGrid3>();
        _advectableVectorDataList.push_back(_velocity);
        _velocityIdx = 0;
        resize(resolution, gridSpacing, origin);
    }

    GridSystemData3::GridSystemData3(const GridSystemData3& other) {
        resize(other._resolution, other._gridSpacing, other._origin);

        for (auto& data : other._scalarDataList) {
            _scalarDataList.push_back(data->clone());
        }
        for (auto& data : other._vectorDataList) {
            _vectorDataList.push_back(data->clone());
        }
        for (auto& data : other._advectableScalarDataList) {
            _advectableScalarDataList.push_back(data->clone());
        }
        for (auto& data : other._advectableVectorDataList) {
            _advectableVectorDataList.push_back(data->clone());
        }

        JET_ASSERT(_advectableVectorDataList.size() > 0);

        _velocity = std::dynamic_pointer_cast<FaceCenteredGrid3>(
            _advectableVectorDataList[0]);

        JET_ASSERT(_velocity != nullptr);

        _velocityIdx = 0;
    }

    GridSystemData3::~GridSystemData3() {
    }

    void GridSystemData3::resize(
        const Size3& resolution,
        const Vector3D& gridSpacing,
        const Vector3D& origin) {
        _resolution = resolution;
        _gridSpacing = gridSpacing;
        _origin = origin;

        for (auto& data : _scalarDataList) {
            data->resize(resolution, gridSpacing, origin);
        }
        for (auto& data : _vectorDataList) {
            data->resize(resolution, gridSpacing, origin);
        }
        for (auto& data : _advectableScalarDataList) {
            data->resize(resolution, gridSpacing, origin);
        }
        for (auto& data : _advectableVectorDataList) {
            data->resize(resolution, gridSpacing, origin);
        }
    }

    Size3 GridSystemData3::resolution() const {
        return _resolution;
    }

    Vector3D GridSystemData3::gridSpacing() const {
        return _gridSpacing;
    }

    Vector3D GridSystemData3::origin() const {
        return _origin;
    }

    BoundingBox3D GridSystemData3::boundingBox() const {
        return _velocity->boundingBox();
    }

    size_t GridSystemData3::addScalarData(
        const ScalarGridBuilder3Ptr& builder,
        double initialVal) {
        size_t attrIdx = _scalarDataList.size();
        _scalarDataList.push_back(
            builder->build(resolution(), gridSpacing(), origin(), initialVal));
        return attrIdx;
    }

    size_t GridSystemData3::addVectorData(
        const VectorGridBuilder3Ptr& builder,
        const Vector3D& initialVal) {
        size_t attrIdx = _vectorDataList.size();
        _vectorDataList.push_back(
            builder->build(resolution(), gridSpacing(), origin(), initialVal));
        return attrIdx;
    }

    size_t GridSystemData3::addAdvectableScalarData(
        const ScalarGridBuilder3Ptr& builder,
        double initialVal) {
        size_t attrIdx = _advectableScalarDataList.size();
        _advectableScalarDataList.push_back(
            builder->build(resolution(), gridSpacing(), origin(), initialVal));
        return attrIdx;
    }

    size_t GridSystemData3::addAdvectableVectorData(
        const VectorGridBuilder3Ptr& builder,
        const Vector3D& initialVal) {
        size_t attrIdx = _advectableVectorDataList.size();
        _advectableVectorDataList.push_back(
            builder->build(resolution(), gridSpacing(), origin(), initialVal));
        return attrIdx;
    }

    const FaceCenteredGrid3Ptr& GridSystemData3::velocity() const {
        return _velocity;
    }

    size_t GridSystemData3::velocityIndex() const {
        return _velocityIdx;
    }

    const ScalarGrid3Ptr& GridSystemData3::scalarDataAt(size_t idx) const {
        return _scalarDataList[idx];
    }

    const VectorGrid3Ptr& GridSystemData3::vectorDataAt(size_t idx) const {
        return _vectorDataList[idx];
    }

    const ScalarGrid3Ptr&
        GridSystemData3::advectableScalarDataAt(size_t idx) const {
        return _advectableScalarDataList[idx];
    }

    const VectorGrid3Ptr&
        GridSystemData3::advectableVectorDataAt(size_t idx) const {
        return _advectableVectorDataList[idx];
    }

    size_t GridSystemData3::numberOfScalarData() const {
        return _scalarDataList.size();
    }

    size_t GridSystemData3::numberOfVectorData() const {
        return _vectorDataList.size();
    }

    size_t GridSystemData3::numberOfAdvectableScalarData() const {
        return _advectableScalarDataList.size();
    }

    size_t GridSystemData3::numberOfAdvectableVectorData() const {
        return _advectableVectorDataList.size();
    }

    void GridSystemData3::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);

        auto resolution = jetToFbs(_resolution);
        auto gridSpacing = jetToFbs(_gridSpacing);
        auto origin = jetToFbs(_origin);

        std::vector<flatbuffers::Offset<fbs::ScalarGridSerialized3>> scalarDataList;
        std::vector<flatbuffers::Offset<fbs::VectorGridSerialized3>> vectorDataList;
        std::vector<flatbuffers::Offset<fbs::ScalarGridSerialized3>>
            advScalarDataList;
        std::vector<flatbuffers::Offset<fbs::VectorGridSerialized3>>
            advVectorDataList;

        serializeGrid(
            &builder,
            _scalarDataList,
            fbs::CreateScalarGridSerialized3,
            &scalarDataList);
        serializeGrid(
            &builder,
            _vectorDataList,
            fbs::CreateVectorGridSerialized3,
            &vectorDataList);
        serializeGrid(
            &builder,
            _advectableScalarDataList,
            fbs::CreateScalarGridSerialized3,
            &advScalarDataList);
        serializeGrid(
            &builder,
            _advectableVectorDataList,
            fbs::CreateVectorGridSerialized3,
            &advVectorDataList);

        auto gsd = fbs::CreateGridSystemData3(
            builder,
            &resolution,
            &gridSpacing,
            &origin,
            _velocityIdx,
            builder.CreateVector(scalarDataList),
            builder.CreateVector(vectorDataList),
            builder.CreateVector(advScalarDataList),
            builder.CreateVector(advVectorDataList));

        builder.Finish(gsd);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void GridSystemData3::deserialize(const std::vector<uint8_t>& buffer) {
        auto gsd = fbs::GetGridSystemData3(buffer.data());

        resize(
            fbsToJet(*gsd->resolution()),
            fbsToJet(*gsd->gridSpacing()),
            fbsToJet(*gsd->origin()));

        _scalarDataList.clear();
        _vectorDataList.clear();
        _advectableScalarDataList.clear();
        _advectableVectorDataList.clear();

        deserializeGrid(
            gsd->scalarData(),
            Factory::buildScalarGrid3,
            &_scalarDataList);
        deserializeGrid(
            gsd->vectorData(),
            Factory::buildVectorGrid3,
            &_vectorDataList);
        deserializeGrid(
            gsd->advectableScalarData(),
            Factory::buildScalarGrid3,
            &_advectableScalarDataList);
        deserializeGrid(
            gsd->advectableVectorData(),
            Factory::buildVectorGrid3,
            &_advectableVectorDataList);

        _velocityIdx = static_cast<size_t>(gsd->velocityIdx());
        _velocity = std::dynamic_pointer_cast<FaceCenteredGrid3>(
            _advectableVectorDataList[_velocityIdx]);
    }

}  // namespace jet

#endif  // INCLUDE_JET_GRID_SYSTEM_DATA3_H_