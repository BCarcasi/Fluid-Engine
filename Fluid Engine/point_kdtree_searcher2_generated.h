#ifndef FLATBUFFERS_GENERATED_POINTKDTREESEARCHER2_JET_FBS_H_
#define FLATBUFFERS_GENERATED_POINTKDTREESEARCHER2_JET_FBS_H_

#include "flatbuffers.h"

#include "basic_types_generated.h"

namespace jet {
    namespace fbs {

        struct PointKdTreeSearcherNode2;

        struct PointKdTreeSearcher2;

        MANUALLY_ALIGNED_STRUCT(8) PointKdTreeSearcherNode2 FLATBUFFERS_FINAL_CLASS {
        private:
            uint64_t flags_;
            uint64_t child_;
            uint64_t item_;

        public:
            PointKdTreeSearcherNode2() {
                memset(this, 0, sizeof(PointKdTreeSearcherNode2));
            }
            PointKdTreeSearcherNode2(const PointKdTreeSearcherNode2 & _o) {
                memcpy(this, &_o, sizeof(PointKdTreeSearcherNode2));
            }
            PointKdTreeSearcherNode2(uint64_t _flags, uint64_t _child, uint64_t _item)
                : flags_(flatbuffers::EndianScalar(_flags)),
                child_(flatbuffers::EndianScalar(_child)),
                item_(flatbuffers::EndianScalar(_item)) {
            }
            uint64_t flags() const {
                return flatbuffers::EndianScalar(flags_);
            }
            uint64_t child() const {
                return flatbuffers::EndianScalar(child_);
            }
            uint64_t item() const {
                return flatbuffers::EndianScalar(item_);
            }
        };
        STRUCT_END(PointKdTreeSearcherNode2, 24);

        struct PointKdTreeSearcher2 FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
            enum {
                VT_POINTS = 4,
                VT_NODES = 6
            };
            const flatbuffers::Vector<const jet::fbs::Vector2D*>* points() const {
                return GetPointer<const flatbuffers::Vector<const jet::fbs::Vector2D*>*>(VT_POINTS);
            }
            const flatbuffers::Vector<const PointKdTreeSearcherNode2*>* nodes() const {
                return GetPointer<const flatbuffers::Vector<const PointKdTreeSearcherNode2*>*>(VT_NODES);
            }
            bool Verify(flatbuffers::Verifier& verifier) const {
                return VerifyTableStart(verifier) &&
                    VerifyOffset(verifier, VT_POINTS) &&
                    verifier.Verify(points()) &&
                    VerifyOffset(verifier, VT_NODES) &&
                    verifier.Verify(nodes()) &&
                    verifier.EndTable();
            }
        };

        struct PointKdTreeSearcher2Builder {
            flatbuffers::FlatBufferBuilder& fbb_;
            flatbuffers::uoffset_t start_;
            void add_points(flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector2D*>> points) {
                fbb_.AddOffset(PointKdTreeSearcher2::VT_POINTS, points);
            }
            void add_nodes(flatbuffers::Offset<flatbuffers::Vector<const PointKdTreeSearcherNode2*>> nodes) {
                fbb_.AddOffset(PointKdTreeSearcher2::VT_NODES, nodes);
            }
            PointKdTreeSearcher2Builder(flatbuffers::FlatBufferBuilder& _fbb)
                : fbb_(_fbb) {
                start_ = fbb_.StartTable();
            }
            PointKdTreeSearcher2Builder& operator=(const PointKdTreeSearcher2Builder&);
            flatbuffers::Offset<PointKdTreeSearcher2> Finish() {
                const auto end = fbb_.EndTable(start_, 2);
                auto o = flatbuffers::Offset<PointKdTreeSearcher2>(end);
                return o;
            }
        };

        inline flatbuffers::Offset<PointKdTreeSearcher2> CreatePointKdTreeSearcher2(
            flatbuffers::FlatBufferBuilder& _fbb,
            flatbuffers::Offset<flatbuffers::Vector<const jet::fbs::Vector2D*>> points = 0,
            flatbuffers::Offset<flatbuffers::Vector<const PointKdTreeSearcherNode2*>> nodes = 0) {
            PointKdTreeSearcher2Builder builder_(_fbb);
            builder_.add_nodes(nodes);
            builder_.add_points(points);
            return builder_.Finish();
        }

        inline flatbuffers::Offset<PointKdTreeSearcher2> CreatePointKdTreeSearcher2Direct(
            flatbuffers::FlatBufferBuilder& _fbb,
            const std::vector<const jet::fbs::Vector2D*>* points = nullptr,
            const std::vector<const PointKdTreeSearcherNode2*>* nodes = nullptr) {
            return jet::fbs::CreatePointKdTreeSearcher2(
                _fbb,
                points ? _fbb.CreateVector<const jet::fbs::Vector2D*>(*points) : 0,
                nodes ? _fbb.CreateVector<const PointKdTreeSearcherNode2*>(*nodes) : 0);
        }

        inline const jet::fbs::PointKdTreeSearcher2* GetPointKdTreeSearcher2(const void* buf) {
            return flatbuffers::GetRoot<jet::fbs::PointKdTreeSearcher2>(buf);
        }

        inline bool VerifyPointKdTreeSearcher2Buffer(
            flatbuffers::Verifier& verifier) {
            return verifier.VerifyBuffer<jet::fbs::PointKdTreeSearcher2>(nullptr);
        }

        inline void FinishPointKdTreeSearcher2Buffer(
            flatbuffers::FlatBufferBuilder& fbb,
            flatbuffers::Offset<jet::fbs::PointKdTreeSearcher2> root) {
            fbb.Finish(root);
        }

    }  // namespace fbs
}  // namespace jet

#endif  // FLATBUFFERS_GENERATED_POINTKDTREESEARCHER2_JET_FBS_H_