#ifndef INCLUDE_JET_PARTICLE_EMITTER_SET3_H_
#define INCLUDE_JET_PARTICLE_EMITTER_SET3_H_


#include "pch.h"
#include "particle_emitter3.h"
#include <tuple>
#include <vector>

namespace jet {

    //!
    //! \brief 3-D particle-based emitter set.
    //!
    class ParticleEmitterSet3 final : public ParticleEmitter3 {
    public:
        class Builder;

        //! Constructs an emitter.
        ParticleEmitterSet3();

        //! Constructs an emitter with sub-emitters.
        explicit ParticleEmitterSet3(
            const std::vector<ParticleEmitter3Ptr>& emitters);

        //! Destructor.
        virtual ~ParticleEmitterSet3();

        //! Adds sub-emitter.
        void addEmitter(const ParticleEmitter3Ptr& emitter);

        //! Returns builder fox ParticleEmitterSet3.
        static Builder builder();

    private:
        std::vector<ParticleEmitter3Ptr> _emitters;

        void onSetTarget(const ParticleSystemData3Ptr& particles) override;

        void onUpdate(
            double currentTimeInSeconds,
            double timeIntervalInSecond) override;
    };

    //! Shared pointer type for the ParticleEmitterSet3.
    typedef std::shared_ptr<ParticleEmitterSet3> ParticleEmitterSet3Ptr;


    //!
    //! \brief Front-end to create ParticleEmitterSet3 objects step by step.
    //!
    class ParticleEmitterSet3::Builder final {
    public:
        //! Returns builder with list of sub-emitters.
        Builder& withEmitters(const std::vector<ParticleEmitter3Ptr>& emitters);

        //! Builds ParticleEmitterSet3.
        ParticleEmitterSet3 build() const;

        //! Builds shared pointer of ParticleEmitterSet3 instance.
        ParticleEmitterSet3Ptr makeShared() const;

    private:
        std::vector<ParticleEmitter3Ptr> _emitters;
    };

    ParticleEmitterSet3::ParticleEmitterSet3() {}

    ParticleEmitterSet3::ParticleEmitterSet3(
        const std::vector<ParticleEmitter3Ptr>& emitters)
        : _emitters(emitters) {}

    ParticleEmitterSet3::~ParticleEmitterSet3() {}

    void ParticleEmitterSet3::addEmitter(const ParticleEmitter3Ptr& emitter) {
        _emitters.push_back(emitter);
    }

    void ParticleEmitterSet3::onSetTarget(const ParticleSystemData3Ptr& particles) {
        for (auto& emitter : _emitters) {
            emitter->setTarget(particles);
        }
    }

    void ParticleEmitterSet3::onUpdate(double currentTimeInSeconds,
        double timeIntervalInSeconds) {
        if (!isEnabled()) {
            return;
        }

        for (auto& emitter : _emitters) {
            emitter->update(currentTimeInSeconds, timeIntervalInSeconds);
        }
    }

    ParticleEmitterSet3::Builder ParticleEmitterSet3::builder() {
        return Builder();
    }

    ParticleEmitterSet3::Builder& ParticleEmitterSet3::Builder::withEmitters(
        const std::vector<ParticleEmitter3Ptr>& emitters) {
        _emitters = emitters;
        return *this;
    }

    ParticleEmitterSet3 ParticleEmitterSet3::Builder::build() const {
        return ParticleEmitterSet3(_emitters);
    }

    ParticleEmitterSet3Ptr ParticleEmitterSet3::Builder::makeShared() const {
        return std::shared_ptr<ParticleEmitterSet3>(
            new ParticleEmitterSet3(_emitters),
            [](ParticleEmitterSet3* obj) { delete obj; });
    }

}  // namespace jet

#endif  // INCLUDE_JET_PARTICLE_EMITTER_SET3_H_