#include <gmock/gmock.h>
#include <thread>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/registries/registry.h"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  mt_kahypar::register_algorithms_and_policies();
  mt_kahypar::TBBInitializer::initialize(std::thread::hardware_concurrency());
  const int result = RUN_ALL_TESTS();
  mt_kahypar::TBBInitializer::instance().terminate();

  return result;
}
