#ifndef FLATBUFFERS_GENERATED_SCALARGRID3_JET_FBS_H_
#define FLATBUFFERS_GENERATED_SCALARGRID3_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"

namespace jet {
    namespace fbs {

        struct ScalarGrid3;

        struct ScalarGrid3 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_RESOLUTION = 4,
                VT_GRIDSPACING = 6,
                VT_ORIGIN = 8,
                VT_DATA = 10
            };
            const jet::fbs::Size3* resolution() const {
                return GetStruct<const jet::fbs::Size3*>(VT_RESOLUTION);
            }
            const jet::fbs::Vector3D* gridSpacing() const {
                return GetStruct<const jet::fbs::Vector3D*>(VT_GRIDSPACING);
            }
            const jet::fbs::Vector3D* origin() const {
                return GetStruct<const jet::fbs::Vector3D*>(VT_ORIGIN);
            }
            const flatbuffers::Vector<double>* data() const {
                return GetPointer<const flatbuffers::Vector<double>*>(VT_DATA);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyField<jet::fbs::Size3>(verifier, VT_RESOLUTION) &&
                    VerifyField<jet::fbs::Vector3D>(verifier, VT_GRIDSPACING) &&
                    VerifyField<jet::fbs::Vector3D>(verifier, VT_ORIGIN) &&
                    VerifyOffset(verifier, VT_DATA) &&
                    verifier.Verify(data()) &&
                    verifier.EndTable();
            }
        };

        struct ScalarGrid3Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_resolution(const jet::fbs::Size3* resolution) {
                fbb_.AddStruct(ScalarGrid3::VT_RESOLUTION, resolution);
            }
            void add_gridSpacing(const jet::fbs::Vector3D* gridSpacing) {
                fbb_.AddStruct(ScalarGrid3::VT_GRIDSPACING, gridSpacing);
            }
            void add_origin(const jet::fbs::Vector3D* origin) {
                fbb_.AddStruct(ScalarGrid3::VT_ORIGIN, origin);
            }
            void add_data(flatbuffers::Offset<flatbuffers::Vector<double>> data) {
                fbb_.AddOffset(ScalarGrid3::VT_DATA, data);
            }
            ScalarGrid3Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            ScalarGrid3Builder& operator=(const ScalarGrid3Builder&);
            flatbuffers::Offset<ScalarGrid3> Finish() {
                const auto end = fbb_.EndTable(start_, 4);
                auto o = flatbuffers::Offset<ScalarGrid3>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<ScalarGrid3> CreateScalarGrid3(
            flatbuffers::FlatBufferBuilder& _fbb,
            const jet::fbs::Size3* resolution = 0,
            const jet::fbs::Vector3D* gridSpacing = 0,
            const jet::fbs::Vector3D* origin = 0,
            flatbuffers::Offset<flatbuffers::Vector<double>> data = 0) {
            ScalarGrid3Builder builder_(_fbb);
            builder_.add_data(data);
            builder_.add_origin(origin);
            builder_.add_gridSpacing(gridSpacing);
            builder_.add_resolution(resolution);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<ScalarGrid3> CreateScalarGrid3Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            const jet::fbs::Size3* resolution = 0,
            const jet::fbs::Vector3D* gridSpacing = 0,
            const jet::fbs::Vector3D* origin = 0,
            const std::vector<double>* data = nullptr) {
            return jet::fbs::CreateScalarGrid3(
                _fbb,
                resolution,
                gridSpacing,
                origin,
                data ? _fbb.CreateVector<double>(*data) : 0);
        }

        inline const jet::fbs::ScalarGrid3* GetScalarGrid3(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::ScalarGrid3>(buf);
        }

        inline bool VerifyScalarGrid3Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::ScalarGrid3>(nullptr);
        }

        inline void FinishScalarGrid3Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::ScalarGrid3> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif  // FLATBUFFERS_GENERATED_SCALARGRID3_JET_FBS_H_