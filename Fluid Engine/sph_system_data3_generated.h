#ifndef FLATBUFFERS_GENERATED_SPHSYSTEMDATA3_JET_FBS_H_
#define FLATBUFFERS_GENERATED_SPHSYSTEMDATA3_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"
#include "particle_system_data3_generated.h"

namespace jet {
    namespace fbs {

        struct SphSystemData3;

        struct SphSystemData3 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_BASE = 4,
                VT_TARGETDENSITY = 6,
                VT_TARGETSPACING = 8,
                VT_KERNELRADIUSOVERTARGETSPACING = 10,
                VT_KERNELRADIUS = 12,
                VT_PRESSUREIDX = 14,
                VT_DENSITYIDX = 16
            };
            const jet::fbs::ParticleSystemData3* base() const {
                return GetPointer<const jet::fbs::ParticleSystemData3*>(VT_BASE);
            }
            double targetDensity() const {
                return GetField<double>(VT_TARGETDENSITY, 0.0);
            }
            double targetSpacing() const {
                return GetField<double>(VT_TARGETSPACING, 0.0);
            }
            double kernelRadiusOverTargetSpacing() const {
                return GetField<double>(VT_KERNELRADIUSOVERTARGETSPACING, 0.0);
            }
            double kernelRadius() const {
                return GetField<double>(VT_KERNELRADIUS, 0.0);
            }
            uint64_t pressureIdx() const {
                return GetField<uint64_t>(VT_PRESSUREIDX, 0);
            }
            uint64_t densityIdx() const {
                return GetField<uint64_t>(VT_DENSITYIDX, 0);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyOffset(verifier, VT_BASE) &&
                    verifier.VerifyTable(base()) &&
                    VerifyField<double>(verifier, VT_TARGETDENSITY) &&
                    VerifyField<double>(verifier, VT_TARGETSPACING) &&
                    VerifyField<double>(verifier, VT_KERNELRADIUSOVERTARGETSPACING) &&
                    VerifyField<double>(verifier, VT_KERNELRADIUS) &&
                    VerifyField<uint64_t>(verifier, VT_PRESSUREIDX) &&
                    VerifyField<uint64_t>(verifier, VT_DENSITYIDX) &&
                    verifier.EndTable();
            }
        };

        struct SphSystemData3Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_base(flatbuffers::Offset<jet::fbs::ParticleSystemData3> base) {
                fbb_.AddOffset(SphSystemData3::VT_BASE, base);
            }
            void add_targetDensity(double targetDensity) {
                fbb_.AddElement<double>(SphSystemData3::VT_TARGETDENSITY, targetDensity, 0.0);
            }
            void add_targetSpacing(double targetSpacing) {
                fbb_.AddElement<double>(SphSystemData3::VT_TARGETSPACING, targetSpacing, 0.0);
            }
            void add_kernelRadiusOverTargetSpacing(double kernelRadiusOverTargetSpacing) {
                fbb_.AddElement<double>(SphSystemData3::VT_KERNELRADIUSOVERTARGETSPACING, kernelRadiusOverTargetSpacing, 0.0);
            }
            void add_kernelRadius(double kernelRadius) {
                fbb_.AddElement<double>(SphSystemData3::VT_KERNELRADIUS, kernelRadius, 0.0);
            }
            void add_pressureIdx(uint64_t pressureIdx) {
                fbb_.AddElement<uint64_t>(SphSystemData3::VT_PRESSUREIDX, pressureIdx, 0);
            }
            void add_densityIdx(uint64_t densityIdx) {
                fbb_.AddElement<uint64_t>(SphSystemData3::VT_DENSITYIDX, densityIdx, 0);
            }
            SphSystemData3Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            SphSystemData3Builder& operator=(const SphSystemData3Builder&);
            flatbuffers::Offset<SphSystemData3> Finish() {
                const auto end = fbb_.EndTable(start_, 7);
                auto o = flatbuffers::Offset<SphSystemData3>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<SphSystemData3> CreateSphSystemData3(
            flatbuffers::FlatBufferBuilder& _fbb,
            flatbuffers::Offset<jet::fbs::ParticleSystemData3> base = 0,
            double targetDensity = 0.0,
            double targetSpacing = 0.0,
            double kernelRadiusOverTargetSpacing = 0.0,
            double kernelRadius = 0.0,
            uint64_t pressureIdx = 0,
            uint64_t densityIdx = 0) {
            SphSystemData3Builder builder_(_fbb);
            builder_.add_densityIdx(densityIdx);
            builder_.add_pressureIdx(pressureIdx);
            builder_.add_kernelRadius(kernelRadius);
            builder_.add_kernelRadiusOverTargetSpacing(kernelRadiusOverTargetSpacing);
            builder_.add_targetSpacing(targetSpacing);
            builder_.add_targetDensity(targetDensity);
            builder_.add_base(base);
            return builder_.Finish();
        }

        inline const jet::fbs::SphSystemData3* GetSphSystemData3(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::SphSystemData3>(buf);
        }

        inline bool VerifySphSystemData3Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::SphSystemData3>(nullptr);
        }

        inline void FinishSphSystemData3Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::SphSystemData3> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif  // FLATBUFFERS_GENERATED_SPHSYSTEMDATA3_JET_FBS_H_