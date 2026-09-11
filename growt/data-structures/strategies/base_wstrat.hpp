/*******************************************************************************
 * data-structures/strategy/wstrat_user.h
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

#include <atomic>
#include <cstddef>
#include <string>

namespace growt
{

template <class Parent> class wstrat
{
  public:
    class global_data_type
    {
      public:
        global_data_type() {}
        global_data_type(const global_data_type& source) = delete;
        global_data_type& operator=(const global_data_type& source) = delete;
        ~global_data_type()                                         = default;
    };


    class local_data_type
    {
      public:
        local_data_type(Parent& parent) : _parent(parent) {}
        local_data_type(const local_data_type& source) = delete;
        local_data_type& operator=(const local_data_type& source) = delete;
        local_data_type(local_data_type&& source)                 = default;
        local_data_type& operator=(local_data_type&& source) = default;
        ~local_data_type()                                   = default;

        Parent& _parent;

        template <class EStrat> inline void init(EStrat&);
        inline void                         deinit();

        template <class ESLocal>
        inline void execute_migration(ESLocal& estrat, size_t);
    };

    static std::string name() { return "w_base"; }
};

} // namespace growt
