/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/io/presets.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/utilities.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/dynamic/dynamic_partitioner.h"

using namespace mt_kahypar;

int main(int argc, char* argv[]) {

  Context context(false);
  processCommandLineInput(context, argc, argv, nullptr);

  if ( context.partition.preset_file.empty() ) {
    if ( context.partition.preset_type != PresetType::UNDEFINED ) {
      // Only a preset type specified => load according preset
      auto preset_option_list = loadPreset(context.partition.preset_type);
      processCommandLineInput(context, argc, argv, &preset_option_list);
    } else {
      throw InvalidInputException("No preset specified");
    }
  }

  // Determine instance (graph or hypergraph) and partition type
  if ( context.partition.instance_type == InstanceType::UNDEFINED ) {
    context.partition.instance_type = to_instance_type(context.partition.file_format);
  }
  context.partition.partition_type = to_partition_c_type(
          context.partition.preset_type, context.partition.instance_type);


  context.utility_id = utils::Utilities::instance().registerNewUtilityObjects();
  if (context.partition.verbose_output) {
    io::printBanner();
  }

  utils::Randomize::instance().setSeed(context.partition.seed);
  if ( context.shared_memory.use_localized_random_shuffle ) {
    utils::Randomize::instance().enableLocalizedParallelShuffle(
            context.shared_memory.shuffle_block_size);
  }

  size_t num_available_cpus = HardwareTopology::instance().num_cpus();
  if ( num_available_cpus < context.shared_memory.num_threads ) {
    WARNING("There are currently only" << num_available_cpus << "cpus available."
                                       << "Setting number of threads from" << context.shared_memory.num_threads
                                       << "to" << num_available_cpus);
    context.shared_memory.num_threads = num_available_cpus;
  }

  // Initialize TBB task arenas on numa nodes
  TBBInitializer::instance(context.shared_memory.num_threads);

  // We set the membind policy to interleaved allocations in order to
  // distribute allocations evenly across NUMA nodes
  hwloc_cpuset_t cpuset = TBBInitializer::instance().used_cpuset();
  parallel::HardwareTopology<>::instance().activate_interleaved_membind_policy(cpuset);
  hwloc_bitmap_free(cpuset);


  dyn::partition(context);

  return 0;
}