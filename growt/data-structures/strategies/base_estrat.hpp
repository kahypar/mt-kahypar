/*******************************************************************************
 * data-structures/strategy/base_estrat.hpp
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
#include <memory>
#include <mutex>
#include <string>

namespace growt
{

template <class table_type>
size_t blockwise_migrate(table_type source, table_type target);


template <class Parent> class estrat
{
  private:
    using this_type   = estrat<Parent>;
    using parent_type = Parent;

  public:
    using base_table_type    = typename Parent::base_table_type;
    using hash_ptr_reference = std::shared_ptr<base_table_type>&;
    using hash_ptr           = std::shared_ptr<base_table_type>;

    // static_assert(std::is_same<typename base_table_type::value_intern,
    // markable_element>::value,
    //               "Asynchroneous migration can only be chosen with
    //               markable_element!!!" );

    class local_data_type;

    // STORED AT THE GLOBAL OBJECT
    class global_data_type
    {
      public:
        global_data_type(size_t size_);
        global_data_type(const global_data_type& source) = delete;
        global_data_type& operator=(const global_data_type& source) = delete;

      private:
        hash_ptr _g_table_w;
    };

    // STORED AT EACH HANDLE
    class local_data_type
    {
      private:
        using worker_strat_local =
            typename Parent::worker_strat::local_data_type;

      public:
        local_data_type(Parent& parent, worker_strat_local& wstrat)
            : _parent(parent), _global(parent._global_exclusion),
              _worker_strat(wstrat)
        {
        }

        local_data_type(const local_data_type& source) = delete;
        local_data_type& operator=(const local_data_type& source) = delete;

        local_data_type(local_data_type&& source) = default;
        local_data_type& operator=(local_data_type&& source) = default;
        ~local_data_type()                                   = default;

        inline void init();
        inline void deinit();

      private:
        Parent&             _parent;
        global_data_type&   _global;
        worker_strat_local& _worker_strat;

      public:
        inline hash_ptr_reference get_table();

        inline void rls_table();

        void          grow();
        void          help_grow();
        inline size_t migrate();
    };

    static std::string name() { return "e_base"; }
};

} // namespace growt
