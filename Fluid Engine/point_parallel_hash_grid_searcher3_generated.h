#ifndef FLATBUFFERS_GENERATED_POINTPARALLELHASHGRIDSEARCHER3_JET_FBS_H_
#define FLATBUFFERS_GENERATED_POINTPARALLELHASHGRIDSEARCHER3_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"

namespace jet {
    namespace fbs {

        struct PointParallelHashGridSearcher3;

        struct PointParallelHashGridSearcher3 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_GRIDSPACING = 4,
                VT_RESOLUTION = 6,
                VT_POINTS = 8,
                VT_KEYS = 10,
                VT_STARTINDEXTABLE = 12,
                VT_ENDINDEXTABLE = 14,
                VT_SORTEDINDICES = 16
            };
            double gridSpacing() const {
                return GetField<double>(VT_GRIDSPACING, 0.0);
            }
            const jet::fbs::Size3* resolution() const {
                return GetStruct<const jet::fbs::Size3*>(VT_RESOLUTION);
            }
            const flatbuffers::Vector<const jet::fbs::Vector3D*>* points() const {
                return GetPointer<const flatbuffers::Vector<const jet::fbs::Vector3D*>*>(VT_POINTS);
            }
            const flatbuffers::Vector<uint64_t>* keys() const {
                return GetPointer<const flatbuffers::Vector<uint64_t>*>(VT_KEYS);
            }
            const flatbuffers::Vector<uint64_t>* startIndexTable() const {
                return GetPointer<const flatbuffers::Vector<uint64_t>*>(VT_STARTINDEXTABLE);
            }
            const flatbuffers::Vector<uint64_t>* endIndexTable() const {
                return GetPointer<const flatbuffers::Vector<uint64_t>*>(VT_ENDINDEXTABLE);
            }
            const flatbuffers::Vector<uint64_t>* sortedIndices() const {
                return GetPointer<const flatbuffers::Vector<uint64_t>*>(VT_SORTEDINDICES);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyField<double>(verifier, VT_GRIDSPACING) &&
                    VerifyField<jet::fbs::Size3>(verifier, VT_RESOLUTION) &&
                    VerifyOffset(verifier, VT_POINTS) &&
                    verifier.Verify(points()) &&
                    VerifyOffset(verifier, VT_KEYS) &&
                    verifier.Verify(keys()) &&
                    VerifyOffset(verifier, VT_STARTINDEXTABLE) &&
                    verifier.Verify(startIndexTable()) &&
                    VerifyOffset(verifier, VT_ENDINDEXTABLE) &&
                    verifier.Verify(endIndexTable()) &&
                    VerifyOffset(verifier, VT_SORTEDINDICES) &&
                    verifier.Verify(sortedIndices()) &&
                    verifier.EndTable();
            }
        };

        struct PointParallelHashGridSearcher3Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_gridSpacing(double gridSpacing) {
                fbb_.AddElement<double>(PointParallelHashGridSearcher3::VT_GRIDSPACING, gridSpacing, 0.0);
            }
            void add_resolution(const jet::fbs::Size3* resolution) {
                fbb_.AddStruct(PointParallelHashGridSearcher3::VT_RESOLUTION, resolution);
            }
            void add_points(flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector3D*>> points) {
                fbb_.AddOffset(PointParallelHashGridSearcher3::VT_POINTS, points);
            }
            void add_keys(flatbuffers::Offset<flatbuffers::Vector<uint64_t>> keys) {
                fbb_.AddOffset(PointParallelHashGridSearcher3::VT_KEYS, keys);
            }
            void add_startIndexTable(flatbuffers::Offset<flatbuffers::Vector<uint64_t>> startIndexTable) {
                fbb_.AddOffset(PointParallelHashGridSearcher3::VT_STARTINDEXTABLE, startIndexTable);
            }
            void add_endIndexTable(flatbuffers::Offset<flatbuffers::Vector<uint64_t>> endIndexTable) {
                fbb_.AddOffset(PointParallelHashGridSearcher3::VT_ENDINDEXTABLE, endIndexTable);
            }
            void add_sortedIndices(flatbuffers::Offset<flatbuffers::Vector<uint64_t>> sortedIndices) {
                fbb_.AddOffset(PointParallelHashGridSearcher3::VT_SORTEDINDICES, sortedIndices);
            }
            PointParallelHashGridSearcher3Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            PointParallelHashGridSearcher3Builder& operator=(const PointParallelHashGridSearcher3Builder&);
            flatbuffers::Offset<PointParallelHashGridSearcher3> Finish() {
                const auto end = fbb_.EndTable(start_, 7);
                auto o = flatbuffers::Offset<PointParallelHashGridSearcher3>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<PointParallelHashGridSearcher3> CreatePointParallelHashGridSearcher3(
            flatbuffers::FlatBufferBuilder& _fbb,
            double gridSpacing = 0.0,
            const jet::fbs::Size3* resolution = 0,
            flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector3D*>> points = 0,
            flatbuffers::Offset<flatbuffers::Vector<uint64_t>> keys = 0,
            flatbuffers::Offset<flatbuffers::Vector<uint64_t>> startIndexTable = 0,
            flatbuffers::Offset<flatbuffers::Vector<uint64_t>> endIndexTable = 0,
            flatbuffers::Offset<flatbuffers::Vector<uint64_t>> sortedIndices = 0) {
            PointParallelHashGridSearcher3Builder builder_(_fbb);
            builder_.add_gridSpacing(gridSpacing);
            builder_.add_sortedIndices(sortedIndices);
            builder_.add_endIndexTable(endIndexTable);
            builder_.add_startIndexTable(startIndexTable);
            builder_.add_keys(keys);
            builder_.add_points(points);
            builder_.add_resolution(resolution);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<PointParallelHashGridSearcher3> CreatePointParallelHashGridSearcher3Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            double gridSpacing = 0.0,
            const jet::fbs::Size3* resolution = 0,
            const std::vector<const jet::fbs::Vector3D*>* points = nullptr,
            const std::vector<uint64_t>* keys = nullptr,
            const std::vector<uint64_t>* startIndexTable = nullptr,
            const std::vector<uint64_t>* endIndexTable = nullptr,
            const std::vector<uint64_t>* sortedIndices = nullptr) {
            return jet::fbs::CreatePointParallelHashGridSearcher3(
                _fbb,
                gridSpacing,
                resolution,
                points ? _fbb.CreateVector<const jet::fbs::Vector3D*>(*points) : 0,
                keys ? _fbb.CreateVector<uint64_t>(*keys) : 0,
                startIndexTable ? _fbb.CreateVector<uint64_t>(*startIndexTable) : 0,
                endIndexTable ? _fbb.CreateVector<uint64_t>(*endIndexTable) : 0,
                sortedIndices ? _fbb.CreateVector<uint64_t>(*sortedIndices) : 0);
        }

        inline const jet::fbs::PointParallelHashGridSearcher3* GetPointParallelHashGridSearcher3(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::PointParallelHashGridSearcher3>(buf);
        }

        inline bool VerifyPointParallelHashGridSearcher3Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::PointParallelHashGridSearcher3>(nullptr);
        }

        inline void FinishPointParallelHashGridSearcher3Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::PointParallelHashGridSearcher3> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif  // FLATBUFFERS_GENERATED_POINTPARALLELHASHGRIDSEARCHER3_JET_FBS_H_