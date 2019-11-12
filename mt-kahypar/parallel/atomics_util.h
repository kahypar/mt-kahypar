/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <atomic>

namespace mt_kahypar {

template<typename T>
void fetch_add(std::atomic<T>& x, T y) {
	T cur_x = x.load();
	while(!x.compare_exchange_weak(cur_x, cur_x + y, std::memory_order_relaxed));
}

template<typename T>
void fetch_sub(std::atomic<T>& x, T y) {
	T cur_x = x.load();
	while(!x.compare_exchange_weak(cur_x, cur_x - y, std::memory_order_relaxed));
}

template<class T>
class AtomicWrapper : public std::atomic<T> {
public:
	void operator+=(T other) {
		fetch_add(*this, other);
	}

	void operator-=(T other) {
		fetch_sub(*this, other);
	}
};

}