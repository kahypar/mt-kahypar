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

  EvoCombineStrategy decideNextCombination(const Context& context, std::mt19937* rng) {
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

    const float basic = context.evolutionary.basic_combine_weight;
    const float edge_frequency = context.evolutionary.edge_frequency_combine_weight;
    const float synthetic_parent = context.evolutionary.synthetic_parent_combine_weight;
    const float total = basic + edge_frequency + synthetic_parent;

    if (total <= 0.0) {
      throw UnsupportedOperationException("No combine operator enabled");
    }

    rand_val = rand_val * total;

    if (rand_val < basic) return EvoCombineStrategy::basic;
    if (rand_val < basic + edge_frequency) return EvoCombineStrategy::edge_frequency;
    return EvoCombineStrategy::synthetic_parent;

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

    if ( rand_val < context.evolutionary.mutation_chance) {
      return EvoDecision::mutation;
    } else {
      return EvoDecision::combine;
    }
  }


   ContextModifierParameters decideContextModificationParameters(const Context& context, std::mt19937* rng) {

    ContextModifierParameters params;
    params.k = context.partition.k;
    params.epsilon = context.partition.epsilon;

    int choice = 0;
    if (context.partition.deterministic) {
      if (!rng) throw UnsupportedOperationException("Deterministic mode requires rng");
      std::uniform_int_distribution<int> dist(0, 4);
      choice = dist(*rng);
    } else {
      choice = utils::Randomize::instance().getRandomInt(0, 4, THREAD_ID);
    }

    switch (choice) {
      case 0: params.epsilon = 3.0 * context.partition.epsilon; break;
      case 1: params.use_random_partitions = true; break;
      case 2: params.use_degree_sorted_partitions = true; break;
      case 3: params.k = 2 * context.partition.k; break;
      case 4: params.recursive_bipartitioning = true; break;
      default: throw InvalidParameterException("Invalid choice for modified combine strategy");
    }
    return params;

  }

}

