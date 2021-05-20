////
//// Created by mlaupichler on 18.05.21.
////
//
//#ifndef KAHYPAR_I_THREAD_LOCAL_ASYNCH_REFINERS_H
//#define KAHYPAR_I_THREAD_LOCAL_ASYNCH_REFINERS_H
//
//#include <memory>
//
//#include "do_nothing_refiner.h"
//#include "label_propagation/asynch_lp_refiner.h"
//
//namespace mt_kahypar {
//
//class IThreadLocalAsynchRefiners {
//
//public:
//    IThreadLocalAsynchRefiners(const IThreadLocalAsynchRefiners&) = delete;
//    IThreadLocalAsynchRefiners(IThreadLocalAsynchRefiners&&) = delete;
//    IThreadLocalAsynchRefiners & operator= (const IThreadLocalAsynchRefiners &) = delete;
//    IThreadLocalAsynchRefiners & operator= (IThreadLocalAsynchRefiners &&) = delete;
//
//    virtual ~IThreadLocalAsynchRefiners() = default;
//
//    IAsynchRefiner& local() {
//        return localImpl();
//    }
//
//protected:
//    IThreadLocalAsynchRefiners() = default;
//
//private:
//    virtual IAsynchRefiner& localImpl() = 0;
//};
//
//template <template <typename> class LocalGainPolicy> class ThreadLocalAsynchLPRefiners : public IThreadLocalAsynchRefiners {
//
//public:
//    ThreadLocalAsynchLPRefiners(Hypergraph &hypergraph, const Context &context, const TaskGroupID task_group_id,
//                                ds::GroupLockManager *lockManager) :
//                                    _hypergraph(hypergraph),
//                                    _context(context),
//                                    _task_group_id(task_group_id),
//                                    _lock_manager(lockManager) {}
//
//private:
//    IAsynchRefiner &localImpl() override {
//        if (!_local_refiner) {
//            _local_refiner = std::make_unique<RefinerType>(_hypergraph, _context,_task_group_id,_lock_manager);
//        }
//        return *_local_refiner;
//    }
//
//private:
//
//    using RefinerType = AsynchLPRefiner<LocalGainPolicy>;
//
//    Hypergraph& _hypergraph;
//    const Context& _context;
//    const TaskGroupID _task_group_id;
//    ds::GroupLockManager* _lock_manager;
//    static thread_local std::unique_ptr<RefinerType> _local_refiner;
//
//};
//
//using ThreadLocalAsynchLPKm1Refiners = ThreadLocalAsynchLPRefiners<LocalKm1Policy>;
//using ThreadLocalAsynchLPCutRefiners = ThreadLocalAsynchLPRefiners<LocalCutPolicy>;
//
//
//class ThreadLocalDoNothingRefiners : public IThreadLocalAsynchRefiners {
//public:
//    ThreadLocalDoNothingRefiners(Hypergraph &hypergraph, const Context &context, const TaskGroupID task_group_id,
//                                 ds::GroupLockManager *lockManager) {};
//
//private:
//    IAsynchRefiner &localImpl() override {
//        if (!_local_refiner) {
//            _local_refiner = std::make_unique<DoNothingRefiner>();
//        }
//        return *_local_refiner;
//    }
//
//    thread_local static std::unique_ptr<DoNothingRefiner> _local_refiner;
//};
//
//
//} // namespace mt_kahypar
//
//#endif //KAHYPAR_I_THREAD_LOCAL_ASYNCH_REFINERS_H
