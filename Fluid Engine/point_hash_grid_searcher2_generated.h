#ifndef FLATBUFFERS_GENERATED_POINTHASHGRIDSEARCHER2_JET_FBS_H_
#define FLATBUFFERS_GENERATED_POINTHASHGRIDSEARCHER2_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"

namespace jet {
    namespace fbs {

        struct PointHashGridSearcherBucket2;

        struct PointHashGridSearcher2;

        struct PointHashGridSearcherBucket2 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_DATA = 4
            };
            const flatbuffers::Vector<uint64_t>* data() const {
                return GetPointer<const flatbuffers::Vector<uint64_t>*>(VT_DATA);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyOffset(verifier, VT_DATA) &&
                    verifier.Verify(data()) &&
                    verifier.EndTable();
            }
        };

        struct PointHashGridSearcherBucket2Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_data(flatbuffers::Offset<flatbuffers::Vector<uint64_t>> data) {
                fbb_.AddOffset(PointHashGridSearcherBucket2::VT_DATA, data);
            }
            PointHashGridSearcherBucket2Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            PointHashGridSearcherBucket2Builder& operator=(const PointHashGridSearcherBucket2Builder&);
            flatbuffers::Offset<PointHashGridSearcherBucket2> Finish() {
                const auto end = fbb_.EndTable(start_, 1);
                auto o = flatbuffers::Offset<PointHashGridSearcherBucket2>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<PointHashGridSearcherBucket2> CreatePointHashGridSearcherBucket2(
            flatbuffers::FlatBufferBuilder& _fbb,
            flatbuffers::Offset<flatbuffers::Vector<uint64_t>> data = 0) {
            PointHashGridSearcherBucket2Builder builder_(_fbb);
            builder_.add_data(data);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<PointHashGridSearcherBucket2> CreatePointHashGridSearcherBucket2Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            const std::vector<uint64_t>* data = nullptr) {
            return jet::fbs::CreatePointHashGridSearcherBucket2(
                _fbb,
                data ? _fbb.CreateVector<uint64_t>(*data) : 0);
        }

        struct PointHashGridSearcher2 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_GRIDSPACING = 4,
                VT_RESOLUTION = 6,
                VT_POINTS = 8,
                VT_BUCKETS = 10
            };
            double gridSpacing() const {
                return GetField<double>(VT_GRIDSPACING, 0.0);
            }
            const jet::fbs::Size2* resolution() const {
                return GetStruct<const jet::fbs::Size2*>(VT_RESOLUTION);
            }
            const flatbuffers::Vector<const jet::fbs::Vector2D*>* points() const {
                return GetPointer<const flatbuffers::Vector<const jet::fbs::Vector2D*>*>(VT_POINTS);
            }
            const flatbuffers::Vector<flatbuffers::Offset<PointHashGridSearcherBucket2>>* buckets() const {
                return GetPointer<const flatbuffers::Vector<flatbuffers::Offset<PointHashGridSearcherBucket2>>*>(VT_BUCKETS);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyField<double>(verifier, VT_GRIDSPACING) &&
                    VerifyField<jet::fbs::Size2>(verifier, VT_RESOLUTION) &&
                    VerifyOffset(verifier, VT_POINTS) &&
                    verifier.Verify(points()) &&
                    VerifyOffset(verifier, VT_BUCKETS) &&
                    verifier.Verify(buckets()) &&
                    verifier.VerifyVectorOfTables(buckets()) &&
                    verifier.EndTable();
            }
        };

        struct PointHashGridSearcher2Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_gridSpacing(double gridSpacing) {
                fbb_.AddElement<double>(PointHashGridSearcher2::VT_GRIDSPACING, gridSpacing, 0.0);
            }
            void add_resolution(const jet::fbs::Size2* resolution) {
                fbb_.AddStruct(PointHashGridSearcher2::VT_RESOLUTION, resolution);
            }
            void add_points(flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector2D*>> points) {
                fbb_.AddOffset(PointHashGridSearcher2::VT_POINTS, points);
            }
            void add_buckets(flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<PointHashGridSearcherBucket2>>> buckets) {
                fbb_.AddOffset(PointHashGridSearcher2::VT_BUCKETS, buckets);
            }
            PointHashGridSearcher2Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            PointHashGridSearcher2Builder& operator=(const PointHashGridSearcher2Builder&);
            flatbuffers::Offset<PointHashGridSearcher2> Finish() {
                const auto end = fbb_.EndTable(start_, 4);
                auto o = flatbuffers::Offset<PointHashGridSearcher2>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<PointHashGridSearcher2> CreatePointHashGridSearcher2(
            flatbuffers::FlatBufferBuilder& _fbb,
            double gridSpacing = 0.0,
            const jet::fbs::Size2* resolution = 0,
            flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector2D*>> points = 0,
            flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<PointHashGridSearcherBucket2>>> buckets = 0) {
            PointHashGridSearcher2Builder builder_(_fbb);
            builder_.add_gridSpacing(gridSpacing);
            builder_.add_buckets(buckets);
            builder_.add_points(points);
            builder_.add_resolution(resolution);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<PointHashGridSearcher2> CreatePointHashGridSearcher2Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            double gridSpacing = 0.0,
            const jet::fbs::Size2* resolution = 0,
            const std::vector<const jet::fbs::Vector2D*>* points = nullptr,
            const std::vector<flatbuffers::Offset<PointHashGridSearcherBucket2>>* buckets = nullptr) {
            return jet::fbs::CreatePointHashGridSearcher2(
                _fbb,
                gridSpacing,
                resolution,
                points ? _fbb.CreateVector<const jet::fbs::Vector2D*>(*points) : 0,
                buckets ? _fbb.CreateVector<flatbuffers::Offset<PointHashGridSearcherBucket2>>(*buckets) : 0);
        }

        inline const jet::fbs::PointHashGridSearcher2* GetPointHashGridSearcher2(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::PointHashGridSearcher2>(buf);
        }

        inline bool VerifyPointHashGridSearcher2Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::PointHashGridSearcher2>(nullptr);
        }

        inline void FinishPointHashGridSearcher2Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::PointHashGridSearcher2> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif // FLATBUFFERS_GENERATED_POINTHASHGRIDSEARCHER2_JET_FBS_H_