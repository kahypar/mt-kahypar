#pragma once

#include "include/libmtkahypartypes.h"

#include "mt-kahypar/utils/reproducible_random.h"
#include "mt-kahypar/utils/utilities.h"

namespace mt_kahypar {
class ClusterTieBreaker {
public:
    virtual ~ClusterTieBreaker() = default;

    size_t select(size_t max, size_t seed) {
        return selectImpl(max, seed);
    }
protected:
    ClusterTieBreaker() = default;
private:
    virtual size_t selectImpl(size_t max, size_t seed) = 0;

};


class SimpleHashUniform final : public ClusterTieBreaker {
private:
    size_t selectImpl(size_t max, size_t seed) {
        hashing::SimpleIntHash<uint32_t> sih;
        hashing::HashRNG<hashing::SimpleIntHash<uint32_t>> hash_prng(sih, seed);
        return std::uniform_int_distribution<uint32_t>(0, max)(hash_prng);
    }
};


class MtUniform final : public ClusterTieBreaker {
private:
    size_t selectImpl(size_t max, size_t seed) {
        std::mt19937 prng(seed);
        return std::uniform_int_distribution<uint32_t>(0, max)(prng);
    }
};

class MtGeometric final : public ClusterTieBreaker {
private:
    size_t selectImpl(size_t max, size_t seed) {
        std::mt19937 prng(seed);
        return  std::geometric_distribution<uint32_t>()(prng) % max;
    }
};

class SimpleHashGeometric final : public ClusterTieBreaker {
private:
    size_t selectImpl(size_t max, size_t seed) {
        hashing::SimpleIntHash<uint32_t> sih;
        hashing::HashRNG<hashing::SimpleIntHash<uint32_t>> hash_prng(sih, seed);
        return  std::geometric_distribution<uint32_t>()(hash_prng) % max;
    }
};

class First final : public ClusterTieBreaker {
private:
    size_t selectImpl(size_t max, size_t seed) {
        unused(seed);
        unused(max);
        return 0;
    }
};

class Last final : public ClusterTieBreaker {

private:
    size_t selectImpl(size_t max, size_t seed) {
        unused(seed);
        unused(max);
        return max;
    }
};
}  // namespace mt_kahypar
