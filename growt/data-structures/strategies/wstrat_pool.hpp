/*******************************************************************************
 * data-structures/strategy/wstrat_pool.h
 *
 * see below
 *
 * Part of Project growt - https://github.com/TooBiased/growt.git
 *
 * Copyright (C) 2015-2016 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#pragma once

#include "counting_wait.hpp"
#include <atomic>
#include <string>
#include <thread>


/*******************************************************************************
 *
 * This is a worker strategy for our growtable.
 *
 * Every worker strategy has to implement the following
 *  - subclass: global_data_type      (is stored at the growtable object)
 *  - subclass: local_data_type       (is stored at each handle)
 *     - init(...)
 *     - deinit()
 *     - execute_migration(...)
 *
 * This specific strategy uses a thread-pool for growing.
 * Every thread who creates a handle will generate a growing
 * thread, which will help with each migration. The thread will
 * be stopped once the handle is deleted.
 *
 * NOTE: The migration thread will be pinned to the core, it was created from.
 *       This is good if all hardware threads are used and pinned .
 *
 ******************************************************************************/

namespace growt
{

template <class Parent> class wstrat_pool
{
  public:
    // Globaly we have to store two "waiting objects" (futexes).
    // They are used to sleep until the next grow/nongrow phase (+wake up).
    class global_data_type
    {
      public:
        global_data_type() : _grow_wait(0), _user_wait(0) {}
        global_data_type(const global_data_type&) = delete;
        global_data_type& operator=(const global_data_type&) = delete;

        ~global_data_type() = default;

        counting_wait _grow_wait;
        counting_wait _user_wait;
    };


    // This is the function executed by the growing threads
    // wait for wakeup -> check if destroyed
    //                 -> check if growing -> help grow -> repeat
    template <class EStrat>
    static void grow_thread_loop(EStrat& estrat, global_data_type& global,
                                 std::atomic_size_t& finished, cpu_set_t* aff);


    // On init the growing thread is created, on deinit it is joined.
    // All migrations are reduced to waiting for the new table version
    // which is automatically created by the thread-pool.
    class local_data_type
    {
      public:
        Parent&                             _parent;
        global_data_type&                   _global;
        std::thread                         _grow_thread;
        std::unique_ptr<std::atomic_size_t> _finished;

        local_data_type(Parent& parent);
        local_data_type(const local_data_type& source) = delete;
        local_data_type& operator=(const local_data_type& source) = delete;
        local_data_type(local_data_type&& rhs);
        local_data_type& operator=(local_data_type&& rhs);
        ~local_data_type() {}

        // creates and pins the thread
        template <class EStrat> inline void init(EStrat& estrat);

        // sets a local destroy flag, and wakes up all growing threads
        // only the one local thread will be destroyed though. Since there is no
        // new table, no migration will be executed by the growing threads.
        inline void deinit();

        template <class ESLocal>
        inline void execute_migration(ESLocal&, size_t epoch);
    };

    static std::string name() { return "w_pool"; }
};


// This is the function executed by the growing threads
// wait for wakeup -> check if destroyed
//                 -> check if growing -> help grow -> repeat
template <class P>
template <class EStrat>
void wstrat_pool<P>::grow_thread_loop(EStrat& estrat, global_data_type& global,
                                      std::atomic_size_t& finished,
                                      cpu_set_t*          aff)
{
    uint epoch = 0;
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), aff);

    while (true)
    {
        global._grow_wait.wait_if(epoch);
        if (finished) break;

        auto next = estrat.migrate();

        global._user_wait.inc_if(epoch);
        global._user_wait.wake();
        epoch = next;
    }
    finished.store(2, std::memory_order_release);
}


template <class P>
wstrat_pool<P>::local_data_type::local_data_type(P& parent)
    : _parent(parent), _global(parent._global_worker),
      _finished(new std::atomic_size_t(0))
{
}


template <class P>
wstrat_pool<P>::local_data_type::local_data_type(local_data_type&& rhs)
    : _parent(rhs._parent), _global(rhs._global),
      _grow_thread(std::move(rhs._grow_thread)), _finished(std::move(_finished))
{
}


template <class P>
typename wstrat_pool<P>::local_data_type&
wstrat_pool<P>::local_data_type::operator=(local_data_type&& rhs)
{
    _parent = rhs._parent;
    _global = rhs._global;
    deinit();
    _grow_thread = std::move(rhs._grow_thread);
    _finished    = std::move(rhs._finished);
}


// creates and pins the thread
template <class P>
template <class EStrat>
void wstrat_pool<P>::local_data_type::init(EStrat& estrat)
{
    cpu_set_t cpuset;
    pthread_getaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

    _grow_thread =
        std::thread(grow_thread_loop<EStrat>, std::ref(estrat),
                    std::ref(_global), std::ref(*_finished), &cpuset);
}


// sets a local destroy flag, and wakes up all growing threads
// only the one local thread will be destroyed though. Since there is no
// new table, no migration will be executed by the growing threads.
template <class P> void wstrat_pool<P>::local_data_type::deinit()
{
    if (_grow_thread.joinable())
    {
        _finished->store(1, std::memory_order_release);

        while (_finished->load(std::memory_order_acquire) < 2)
            _global._grow_wait.wake();

        _grow_thread.join();
    }
}


template <class P>
template <class ESLocal>
void wstrat_pool<P>::local_data_type::execute_migration(ESLocal&, size_t epoch)
{
    // lets instead tell somebody else and ...
    // wait lazily until somebody did this zzzzZZZzz
    if (_global._grow_wait.inc_if(epoch)) _global._grow_wait.wake();

    _global._user_wait.wait_if(epoch);
}

} // namespace growt
