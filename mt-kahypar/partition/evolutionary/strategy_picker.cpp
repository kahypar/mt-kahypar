#include "mt-kahypar/partition/evolutionary/strategy_picker.h"

namespace mt_kahypar::pick {
  EvoMutateStrategy decideNextMutation(const Context& context, std::mt19937* rng) {
    if (context.partition.deterministic) {
      if (rng == nullptr) {
        throw UnsupportedOperationException("Catastrophic Error! Deterministic mode requires passing rng!");
      }
      std::uniform_real_distribution<float> dist(0.0f, 1.0f);
      if (const float rand_val = dist(*rng); rand_val < 0.5f ) {
        return EvoMutateStrategy::vcycle;
      }
      return EvoMutateStrategy::new_initial_partitioning_vcycle;
    }
    else {
      if (utils::Randomize::instance().flipCoin(THREAD_ID)) {
        return EvoMutateStrategy::vcycle;
      }
      return EvoMutateStrategy::new_initial_partitioning_vcycle;
    }
  }

  EvoDecision decideNextMove(const Context& context, std::mt19937* rng) {

    float rand_val;
    if (!context.partition.deterministic)
      rand_val = utils::Randomize::instance().getRandomFloat(0, 1, THREAD_ID);
    else if (rng != nullptr) {
      std::uniform_real_distribution<float> dist(0.0f, 1.0f);
      rand_val = dist(*rng);
    }
    else {
      throw UnsupportedOperationException("Catastrophic Error! Deterministic mode requires passing rng!");
    }

    if (!context.evolutionary.enable_modified_combine) {
      if ( rand_val < context.evolutionary.mutation_chance ) {
        return EvoDecision::mutation;
      }
      return EvoDecision::combine;
    }
    else {
      if ( rand_val < context.evolutionary.mutation_chance) {
        return EvoDecision::mutation;
      }
      else if ( rand_val < context.evolutionary.mutation_chance + context.evolutionary.modified_combine_chance ) {
        return EvoDecision::modified_combine;
      }
      else {
        return EvoDecision::combine;
      }
    }

    return EvoDecision::combine;
  }

}

