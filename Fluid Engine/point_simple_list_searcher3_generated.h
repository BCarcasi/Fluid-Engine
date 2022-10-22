#ifndef FLATBUFFERS_GENERATED_POINTSIMPLELISTSEARCHER3_JET_FBS_H_
#define FLATBUFFERS_GENERATED_POINTSIMPLELISTSEARCHER3_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"

namespace jet {
    namespace fbs {

        struct PointSimpleListSearcher3;

        struct PointSimpleListSearcher3 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_POINTS = 4
            };
            const flatbuffers::Vector<const jet::fbs::Vector3D*>* points() const {
                return GetPointer<const flatbuffers::Vector<const jet::fbs::Vector3D*>*>(VT_POINTS);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyOffset(verifier, VT_POINTS) &&
                    verifier.Verify(points()) &&
                    verifier.EndTable();
            }
        };

        struct PointSimpleListSearcher3Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_points(flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector3D*>> points) {
                fbb_.AddOffset(PointSimpleListSearcher3::VT_POINTS, points);
            }
            PointSimpleListSearcher3Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            PointSimpleListSearcher3Builder& operator=(const PointSimpleListSearcher3Builder&);
            flatbuffers::Offset<PointSimpleListSearcher3> Finish() {
                const auto end = fbb_.EndTable(start_, 1);
                auto o = flatbuffers::Offset<PointSimpleListSearcher3>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<PointSimpleListSearcher3> CreatePointSimpleListSearcher3(
            flatbuffers::FlatBufferBuilder& _fbb,
            flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector3D*>> points = 0) {
            PointSimpleListSearcher3Builder builder_(_fbb);
            builder_.add_points(points);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<PointSimpleListSearcher3> CreatePointSimpleListSearcher3Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            const std::vector<const jet::fbs::Vector3D*>* points = nullptr) {
            return jet::fbs::CreatePointSimpleListSearcher3(
                _fbb,
                points ? _fbb.CreateVector<const jet::fbs::Vector3D*>(*points) : 0);
        }

        inline const jet::fbs::PointSimpleListSearcher3* GetPointSimpleListSearcher3(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::PointSimpleListSearcher3>(buf);
        }

        inline bool VerifyPointSimpleListSearcher3Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::PointSimpleListSearcher3>(nullptr);
        }

        inline void FinishPointSimpleListSearcher3Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::PointSimpleListSearcher3> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif  // FLATBUFFERS_GENERATED_POINTSIMPLELISTSEARCHER3_JET_FBS_H_