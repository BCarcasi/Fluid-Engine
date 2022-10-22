#ifndef INCLUDE_JET_PARTICLE_SYSTEM_DATA2_H_
#define INCLUDE_JET_PARTICLE_SYSTEM_DATA2_H_

#include "array1.h"
#include "point_neighbor_searcher2.h"
#include "serialization.h"

#include <memory>
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable: 4244)
#endif

#include "pch.h"

#include "factory.h"
#include "fbs_helpers.h"
#include "particle_system_data2_generated.h"

#include "parallel.h"
#include "point_parallel_hash_grid_searcher2.h"
#include "timer.h"

#include <algorithm>
#include <vector>

#ifndef JET_DOXYGEN

namespace flatbuffers {

    class FlatBufferBuilder;
    template<typename T> struct Offset;

}

namespace jet {
    namespace fbs {

        struct ParticleSystemData2;

    }
}

#endif  // JET_DOXYGEN

namespace jet {

    //!
    //! \brief      2-D particle system data.
    //!
    //! This class is the key data structure for storing particle system data. A
    //! single particle has position, velocity, and force attributes by default. But
    //! it can also have additional custom scalar or vector attributes.
    //!
    class ParticleSystemData2 : public Serializable {
    public:
        //! Scalar data chunk.
        typedef Array1<double> ScalarData;

        //! Vector data chunk.
        typedef Array1<Vector2D> VectorData;

        //! Default constructor.
        ParticleSystemData2();

        //! Constructs particle system data with given number of particles.
        explicit ParticleSystemData2(size_t numberOfParticles);

        //! Copy constructor.
        ParticleSystemData2(const ParticleSystemData2& other);

        //! Destructor.
        virtual ~ParticleSystemData2();

        //!
        //! \brief      Resizes the number of particles of the container.
        //!
        //! This function will resize internal containers to store newly given
        //! number of particles including custom data layers. However, this will
        //! invalidate neighbor searcher and neighbor lists. It is users
        //! responsibility to call ParticleSystemData2::buildNeighborSearcher and
        //! ParticleSystemData2::buildNeighborLists to refresh those data.
        //!
        //! \param[in]  newNumberOfParticles    New number of particles.
        //!
        void resize(size_t newNumberOfParticles);

        //! Returns the number of particles.
        size_t numberOfParticles() const;

        //!
        //! \brief      Adds a scalar data layer and returns its index.
        //!
        //! This function adds a new scalar data layer to the system. It can be used
        //! for adding a scalar attribute, such as temperature, to the particles.
        //!
        //! \params[in] initialVal  Initial value of the new scalar data.
        //!
        size_t addScalarData(double initialVal = 0.0);

        //!
        //! \brief      Adds a vector data layer and returns its index.
        //!
        //! This function adds a new vector data layer to the system. It can be used
        //! for adding a vector attribute, such as vortex, to the particles.
        //!
        //! \params[in] initialVal  Initial value of the new vector data.
        //!
        size_t addVectorData(const Vector2D& initialVal = Vector2D());

        //! Returns the radius of the particles.
        double radius() const;

        //! Sets the radius of the particles.
        virtual void setRadius(double newRadius);

        //! Returns the mass of the particles.
        double mass() const;

        //! Sets the mass of the particles.
        virtual void setMass(double newMass);

        //! Returns the position array (immutable).
        ConstArrayAccessor1<Vector2D> positions() const;

        //! Returns the position array (mutable).
        ArrayAccessor1<Vector2D> positions();

        //! Returns the velocity array (immutable).
        ConstArrayAccessor1<Vector2D> velocities() const;

        //! Returns the velocity array (mutable).
        ArrayAccessor1<Vector2D> velocities();

        //! Returns the force array (immutable).
        ConstArrayAccessor1<Vector2D> forces() const;

        //! Returns the force array (mutable).
        ArrayAccessor1<Vector2D> forces();

        //! Returns custom scalar data layer at given index (immutable).
        ConstArrayAccessor1<double> scalarDataAt(size_t idx) const;

        //! Returns custom scalar data layer at given index (mutable).
        ArrayAccessor1<double> scalarDataAt(size_t idx);

        //! Returns custom vector data layer at given index (immutable).
        ConstArrayAccessor1<Vector2D> vectorDataAt(size_t idx) const;

        //! Returns custom vector data layer at given index (mutable).
        ArrayAccessor1<Vector2D> vectorDataAt(size_t idx);

        //!
        //! \brief      Adds a particle to the data structure.
        //!
        //! This function will add a single particle to the data structure. For
        //! custom data layers, zeros will be assigned for new particles.
        //! However, this will invalidate neighbor searcher and neighbor lists. It
        //! is users responsibility to call
        //! ParticleSystemData2::buildNeighborSearcher and
        //! ParticleSystemData2::buildNeighborLists to refresh those data.
        //!
        //! \param[in]  newPosition The new position.
        //! \param[in]  newVelocity The new velocity.
        //! \param[in]  newForce    The new force.
        //!
        void addParticle(
            const Vector2D& newPosition,
            const Vector2D& newVelocity = Vector2D(),
            const Vector2D& newForce = Vector2D());

        //!
        //! \brief      Adds particles to the data structure.
        //!
        //! This function will add particles to the data structure. For custom data
        //! layers, zeros will be assigned for new particles. However, this will
        //! invalidate neighbor searcher and neighbor lists. It is users
        //! responsibility to call ParticleSystemData2::buildNeighborSearcher and
        //! ParticleSystemData2::buildNeighborLists to refresh those data.
        //!
        //! \param[in]  newPositions  The new positions.
        //! \param[in]  newVelocities The new velocities.
        //! \param[in]  newForces     The new forces.
        //!
        void addParticles(
            const ConstArrayAccessor1<Vector2D>& newPositions,
            const ConstArrayAccessor1<Vector2D>& newVelocities
            = ConstArrayAccessor1<Vector2D>(),
            const ConstArrayAccessor1<Vector2D>& newForces
            = ConstArrayAccessor1<Vector2D>());

        //!
        //! \brief      Returns neighbor searcher.
        //!
        //! This function returns currently set neighbor searcher object. By
        //! default, PointParallelHashGridSearcher2 is used.
        //!
        //! \return     Current neighbor searcher.
        //!
        const PointNeighborSearcher2Ptr& neighborSearcher() const;

        //! Sets neighbor searcher.
        void setNeighborSearcher(
            const PointNeighborSearcher2Ptr& newNeighborSearcher);

        //!
        //! \brief      Returns neighbor lists.
        //!
        //! This function returns neighbor lists which is available after calling
        //! PointParallelHashGridSearcher2::buildNeighborLists. Each list stores
        //! indices of the neighbors.
        //!
        //! \return     Neighbor lists.
        //!
        const std::vector<std::vector<size_t>>& neighborLists() const;

        //! Builds neighbor searcher with given search radius.
        void buildNeighborSearcher(double maxSearchRadius);

        //! Builds neighbor lists with given search radius.
        void buildNeighborLists(double maxSearchRadius);

        //! Serializes this particle system data to the buffer.
        void serialize(std::vector<uint8_t>* buffer) const override;

        //! Deserializes this particle system data from the buffer.
        void deserialize(const std::vector<uint8_t>& buffer) override;

        //! Copies from other particle system data.
        void set(const ParticleSystemData2& other);

        //! Copies from other particle system data.
        ParticleSystemData2& operator=(const ParticleSystemData2& other);

    protected:
        void serializeParticleSystemData(
            flatbuffers::FlatBufferBuilder* builder,
            flatbuffers::Offset<fbs::ParticleSystemData2>* fbsParticleSystemData)
            const;

        void deserializeParticleSystemData(
            const fbs::ParticleSystemData2* fbsParticleSystemData);

    private:
        double _radius = 1e-3;
        double _mass = 1e-3;
        size_t _numberOfParticles = 0;
        size_t _positionIdx;
        size_t _velocityIdx;
        size_t _forceIdx;

        std::vector<ScalarData> _scalarDataList;
        std::vector<VectorData> _vectorDataList;

        PointNeighborSearcher2Ptr _neighborSearcher;
        std::vector<std::vector<size_t>> _neighborLists;
    };

    //! Shared pointer type of ParticleSystemData2.
    typedef std::shared_ptr<ParticleSystemData2> ParticleSystemData2Ptr;

    static const size_t kDefaultHashGridResolutionData2 = 64;

    ParticleSystemData2::ParticleSystemData2()
        : ParticleSystemData2(0) {
    }

    ParticleSystemData2::ParticleSystemData2(size_t numberOfParticles) {
        _positionIdx = addVectorData();
        _velocityIdx = addVectorData();
        _forceIdx = addVectorData();

        // Use PointParallelHashGridSearcher2 by default
        _neighborSearcher = std::make_shared<PointParallelHashGridSearcher2>(
            kDefaultHashGridResolutionData2,
            kDefaultHashGridResolutionData2,
            2.0 * _radius);

        resize(numberOfParticles);
    }

    ParticleSystemData2::ParticleSystemData2(const ParticleSystemData2& other) {
        set(other);
    }

    ParticleSystemData2::~ParticleSystemData2() {
    }

    void ParticleSystemData2::resize(size_t newNumberOfParticles) {
        _numberOfParticles = newNumberOfParticles;

        for (auto& attr : _scalarDataList) {
            attr.resize(newNumberOfParticles, 0.0);
        }

        for (auto& attr : _vectorDataList) {
            attr.resize(newNumberOfParticles, Vector2D());
        }
    }

    size_t ParticleSystemData2::numberOfParticles() const {
        return _numberOfParticles;
    }

    size_t ParticleSystemData2::addScalarData(double initialVal) {
        size_t attrIdx = _scalarDataList.size();
        _scalarDataList.emplace_back(numberOfParticles(), initialVal);
        return attrIdx;
    }

    size_t ParticleSystemData2::addVectorData(const Vector2D& initialVal) {
        size_t attrIdx = _vectorDataList.size();
        _vectorDataList.emplace_back(numberOfParticles(), initialVal);
        return attrIdx;
    }

    double ParticleSystemData2::radius() const {
        return _radius;
    }

    void ParticleSystemData2::setRadius(double newRadius) {
        _radius = std::max(newRadius, 0.0);
    }

    double ParticleSystemData2::mass() const {
        return _mass;
    }

    void ParticleSystemData2::setMass(double newMass) {
        _mass = std::max(newMass, 0.0);
    }

    ConstArrayAccessor1<Vector2D> ParticleSystemData2::positions() const {
        return vectorDataAt(_positionIdx);
    }

    ArrayAccessor1<Vector2D> ParticleSystemData2::positions() {
        return vectorDataAt(_positionIdx);
    }

    ConstArrayAccessor1<Vector2D> ParticleSystemData2::velocities() const {
        return vectorDataAt(_velocityIdx);
    }

    ArrayAccessor1<Vector2D> ParticleSystemData2::velocities() {
        return vectorDataAt(_velocityIdx);
    }

    ConstArrayAccessor1<Vector2D> ParticleSystemData2::forces() const {
        return vectorDataAt(_forceIdx);
    }

    ArrayAccessor1<Vector2D> ParticleSystemData2::forces() {
        return vectorDataAt(_forceIdx);
    }

    ConstArrayAccessor1<double> ParticleSystemData2::scalarDataAt(
        size_t idx) const {
        return _scalarDataList[idx].constAccessor();
    }

    ArrayAccessor1<double> ParticleSystemData2::scalarDataAt(size_t idx) {
        return _scalarDataList[idx].accessor();
    }

    ConstArrayAccessor1<Vector2D> ParticleSystemData2::vectorDataAt(
        size_t idx) const {
        return _vectorDataList[idx].constAccessor();
    }

    ArrayAccessor1<Vector2D> ParticleSystemData2::vectorDataAt(size_t idx) {
        return _vectorDataList[idx].accessor();
    }

    void ParticleSystemData2::addParticle(
        const Vector2D& newPosition,
        const Vector2D& newVelocity,
        const Vector2D& newForce) {
        Array1<Vector2D> newPositions = { newPosition };
        Array1<Vector2D> newVelocities = { newVelocity };
        Array1<Vector2D> newForces = { newForce };

        addParticles(
            newPositions.constAccessor(),
            newVelocities.constAccessor(),
            newForces.constAccessor());
    }

    void ParticleSystemData2::addParticles(
        const ConstArrayAccessor1<Vector2D>& newPositions,
        const ConstArrayAccessor1<Vector2D>& newVelocities,
        const ConstArrayAccessor1<Vector2D>& newForces) {
        JET_THROW_INVALID_ARG_IF(
            newVelocities.size() > 0
            && newVelocities.size() != newPositions.size());
        JET_THROW_INVALID_ARG_IF(
            newForces.size() > 0 && newForces.size() != newPositions.size());

        size_t oldNumberOfParticles = numberOfParticles();
        size_t newNumberOfParticles = oldNumberOfParticles + newPositions.size();

        resize(newNumberOfParticles);

        auto pos = positions();
        auto vel = velocities();
        auto frc = forces();

        parallelFor(kZeroSize, newPositions.size(),
            [&](size_t i) {
                pos[i + oldNumberOfParticles] = newPositions[i];
            });

        if (newVelocities.size() > 0) {
            parallelFor(kZeroSize, newPositions.size(),
                [&](size_t i) {
                    vel[i + oldNumberOfParticles] = newVelocities[i];
                });
        }

        if (newForces.size() > 0) {
            parallelFor(kZeroSize, newPositions.size(),
                [&](size_t i) {
                    frc[i + oldNumberOfParticles] = newForces[i];
                });
        }
    }

    const PointNeighborSearcher2Ptr& ParticleSystemData2::neighborSearcher() const {
        return _neighborSearcher;
    }

    void ParticleSystemData2::setNeighborSearcher(
        const PointNeighborSearcher2Ptr& newNeighborSearcher) {
        _neighborSearcher = newNeighborSearcher;
    }

    const std::vector<std::vector<size_t>>&
        ParticleSystemData2::neighborLists() const {
        return _neighborLists;
    }

    void ParticleSystemData2::buildNeighborSearcher(double maxSearchRadius) {
        Timer timer;

        // Use PointParallelHashGridSearcher2 by default
        _neighborSearcher = std::make_shared<PointParallelHashGridSearcher2>(
            kDefaultHashGridResolutionData2,
            kDefaultHashGridResolutionData2,
            2.0 * maxSearchRadius);

        _neighborSearcher->build(positions());

        JET_INFO << "Building neighbor searcher took: "
            << timer.durationInSeconds()
            << " seconds";
    }

    void ParticleSystemData2::buildNeighborLists(double maxSearchRadius) {
        Timer timer;

        _neighborLists.resize(numberOfParticles());

        auto points = positions();
        for (size_t i = 0; i < numberOfParticles(); ++i) {
            Vector2D origin = points[i];
            _neighborLists[i].clear();

            _neighborSearcher->forEachNearbyPoint(
                origin,
                maxSearchRadius,
                [&](size_t j, const Vector2D&) {
                    if (i != j) {
                        _neighborLists[i].push_back(j);
                    }
                });
        }

        JET_INFO << "Building neighbor list took: "
            << timer.durationInSeconds()
            << " seconds";
    }

    void ParticleSystemData2::serialize(std::vector<uint8_t>* buffer) const {
        flatbuffers::FlatBufferBuilder builder(1024);
        flatbuffers::Offset<fbs::ParticleSystemData2> fbsParticleSystemData;

        serializeParticleSystemData(&builder, &fbsParticleSystemData);

        builder.Finish(fbsParticleSystemData);

        uint8_t* buf = builder.GetBufferPointer();
        size_t size = builder.GetSize();

        buffer->resize(size);
        memcpy(buffer->data(), buf, size);
    }

    void ParticleSystemData2::deserialize(const std::vector<uint8_t>& buffer) {
        auto fbsParticleSystemData = fbs::GetParticleSystemData2(buffer.data());
        deserializeParticleSystemData(fbsParticleSystemData);
    }

    void ParticleSystemData2::set(const ParticleSystemData2& other) {
        _radius = other._radius;
        _mass = other._mass;
        _positionIdx = other._positionIdx;
        _velocityIdx = other._velocityIdx;
        _forceIdx = other._forceIdx;
        _numberOfParticles = other._numberOfParticles;

        for (auto& attr : other._scalarDataList) {
            _scalarDataList.emplace_back(attr);
        }

        for (auto& attr : other._vectorDataList) {
            _vectorDataList.emplace_back(attr);
        }

        _neighborSearcher = other._neighborSearcher->clone();
        _neighborLists = other._neighborLists;
    }

    ParticleSystemData2& ParticleSystemData2::operator=(
        const ParticleSystemData2& other) {
        set(other);
        return *this;
    }

    void ParticleSystemData2::serializeParticleSystemData(
        flatbuffers::FlatBufferBuilder* builder,
        flatbuffers::Offset<fbs::ParticleSystemData2>* fbsParticleSystemData)
        const {
        // Copy data
        std::vector<flatbuffers::Offset<fbs::ScalarParticleData2>> scalarDataList;
        for (const auto& scalarData : _scalarDataList) {
            auto fbsScalarData = fbs::CreateScalarParticleData2(
                *builder,
                builder->CreateVector(scalarData.data(), scalarData.size()));
            scalarDataList.push_back(fbsScalarData);
        }
        auto fbsScalarDataList = builder->CreateVector(scalarDataList);

        std::vector<flatbuffers::Offset<fbs::VectorParticleData2>> vectorDataList;
        for (const auto& vectorData : _vectorDataList) {
            std::vector<fbs::Vector2D> newVectorData;
            for (const auto& v : vectorData) {
                newVectorData.push_back(jetToFbs(v));
            }

            auto fbsVectorData = fbs::CreateVectorParticleData2(
                *builder,
                builder->CreateVectorOfStructs(
                    newVectorData.data(), newVectorData.size()));
            vectorDataList.push_back(fbsVectorData);
        }
        auto fbsVectorDataList = builder->CreateVector(vectorDataList);

        // Copy neighbor searcher
        auto neighborSearcherType
            = builder->CreateString(_neighborSearcher->typeName());
        std::vector<uint8_t> neighborSearcherSerialized;
        _neighborSearcher->serialize(&neighborSearcherSerialized);
        auto fbsNeighborSearcher = fbs::CreatePointNeighborSearcherSerialized2(
            *builder,
            neighborSearcherType,
            builder->CreateVector(
                neighborSearcherSerialized.data(),
                neighborSearcherSerialized.size()));

        // Copy neighbor lists
        std::vector<flatbuffers::Offset<fbs::ParticleNeighborList2>> neighborLists;
        for (const auto& neighbors : _neighborLists) {
            std::vector<uint64_t> neighbors64(neighbors.begin(), neighbors.end());
            flatbuffers::Offset<fbs::ParticleNeighborList2> fbsNeighborList
                = fbs::CreateParticleNeighborList2(
                    *builder,
                    builder->CreateVector(neighbors64.data(), neighbors64.size()));
            neighborLists.push_back(fbsNeighborList);
        }

        auto fbsNeighborLists = builder->CreateVector(neighborLists);

        // Copy the searcher
        *fbsParticleSystemData = fbs::CreateParticleSystemData2(
            *builder,
            _radius,
            _mass,
            _positionIdx,
            _velocityIdx,
            _forceIdx,
            fbsScalarDataList,
            fbsVectorDataList,
            fbsNeighborSearcher,
            fbsNeighborLists);
    }

    void ParticleSystemData2::deserializeParticleSystemData(
        const fbs::ParticleSystemData2* fbsParticleSystemData) {
        _scalarDataList.clear();
        _vectorDataList.clear();

        // Copy scalars
        _radius = fbsParticleSystemData->radius();
        _mass = fbsParticleSystemData->mass();
        _positionIdx = static_cast<size_t>(fbsParticleSystemData->positionIdx());
        _velocityIdx = static_cast<size_t>(fbsParticleSystemData->velocityIdx());
        _forceIdx = static_cast<size_t>(fbsParticleSystemData->forceIdx());

        // Copy data
        auto fbsScalarDataList = fbsParticleSystemData->scalarDataList();
        for (const auto& fbsScalarData : (*fbsScalarDataList)) {
            auto data = fbsScalarData->data();

            _scalarDataList.push_back(ScalarData(data->size()));

            auto& newData = *(_scalarDataList.rbegin());

            for (uint32_t i = 0; i < data->size(); ++i) {
                newData[i] = data->Get(i);
            }
        }

        auto fbsVectorDataList = fbsParticleSystemData->vectorDataList();
        for (const auto& fbsVectorData : (*fbsVectorDataList)) {
            auto data = fbsVectorData->data();

            _vectorDataList.push_back(VectorData(data->size()));
            auto& newData = *(_vectorDataList.rbegin());
            for (uint32_t i = 0; i < data->size(); ++i) {
                newData[i] = fbsToJet(*data->Get(i));
            }
        }

        _numberOfParticles = _vectorDataList[0].size();

        // Copy neighbor searcher
        auto fbsNeighborSearcher = fbsParticleSystemData->neighborSearcher();
        _neighborSearcher
            = Factory::buildPointNeighborSearcher2(
                fbsNeighborSearcher->type()->c_str());
        std::vector<uint8_t> neighborSearcherSerialized(
            fbsNeighborSearcher->data()->begin(),
            fbsNeighborSearcher->data()->end());
        _neighborSearcher->deserialize(neighborSearcherSerialized);

        // Copy neighbor list
        auto fbsNeighborLists = fbsParticleSystemData->neighborLists();
        _neighborLists.resize(fbsNeighborLists->size());
        for (uint32_t i = 0; i < fbsNeighborLists->size(); ++i) {
            auto fbsNeighborList = fbsNeighborLists->Get(i);
            _neighborLists[i].resize(fbsNeighborList->data()->size());
            std::transform(
                fbsNeighborList->data()->begin(),
                fbsNeighborList->data()->end(),
                _neighborLists[i].begin(),
                [](uint64_t val) {
                    return static_cast<size_t>(val);
                });
        }
    }

}  // namespace jet

#endif  // INCLUDE_JET_PARTICLE_SYSTEM_DATA2_H_