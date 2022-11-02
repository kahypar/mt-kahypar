/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2021 Noah Wahl <noah.wahl@student.kit.edu>
 * Copyright (C) 2021 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2021 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"

namespace mt_kahypar {

  class IUncoarsener {

  public:
    IUncoarsener(const IUncoarsener&) = delete;
    IUncoarsener(IUncoarsener&&) = delete;
    IUncoarsener & operator= (const IUncoarsener &) = delete;
    IUncoarsener & operator= (IUncoarsener &&) = delete;

    PartitionedHypergraph&& uncoarsen(std::unique_ptr<IRefiner>& label_propagation,
                                      std::unique_ptr<IRefiner>& fm) {
      return uncoarsenImpl(label_propagation, fm);
    }

    virtual ~IUncoarsener() = default;

  protected:
    IUncoarsener() = default;

  private:
    virtual PartitionedHypergraph&& uncoarsenImpl(std::unique_ptr<IRefiner>& label_propagation,
                                                  std::unique_ptr<IRefiner>& fm) = 0;
  };
}
