/*******************************************************************************
 * utils/counting_wait.h
 *
 * Simple counter that allows threads to wait until the counter is changed.
 * Waiting threads are sleeping using the futex syscall.
 *
 * Part of Project growt - https://github.com/TooBiased/growt.git
 *
 * Copyright (C) 2015-2016 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#ifndef COUNTING_WAIT_H
#define COUNTING_WAIT_H


#include <atomic>
#include <iostream>

namespace growt
{

class alignas(64) counting_wait
{
  public:
    inline counting_wait(int start = 0) : counter(start)
    {
        if (sizeof(std::atomic_int) != 4)
            std::cout << "std::atomic_int has wrong size:"
                      << sizeof(std::atomic_int) << std::endl;
    }

    inline bool inc_if(int exp)
    {
        auto temp = exp;
        return counter.compare_exchange_strong(temp, exp + 1,
                                               std::memory_order_acq_rel,
                                               std::memory_order_acquire);
    }

    inline void wait_if(int exp)
    {
        // while (counter.load(std::memory_order_acquire) < l_epoch) ;
        // //temporary should soon be removed
        counter.wait(exp);
    }

    inline void wake()
    {
        counter.notify_all();
    }

  private:
    std::atomic_int counter;
};

} // namespace growt

#endif // COUNTING_WAIT_H
