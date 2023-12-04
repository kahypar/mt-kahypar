#include "hypergraph_statistics.h"
#include "array"
namespace mt_kahypar{
    namespace utils{
        std::array<double, mt_kahypar::dimension> parallel_stdev_hnode(const std::vector<HypernodeWeight>& data, const std::array<double, mt_kahypar::dimension> avg, const size_t n) {
    std::array<double, mt_kahypar::dimension> init;
    init.fill(0.0);
    std::array<double, mt_kahypar::dimension> res = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(UL(0), data.size()), init,
            [&](tbb::blocked_range<size_t>& range, const std::array<double, mt_kahypar::dimension> init) -> std::array<double, mt_kahypar::dimension> {
            std::array<double,mt_kahypar::dimension> tmp_stdev = init;
            for ( size_t i = range.begin(); i < range.end(); ++i ) {
                for(int j = 0; j < mt_kahypar::dimension; j++){
                    tmp_stdev[j] += (data[i].weights[j] - avg[j]) * (data[i].weights[j] - avg[j]);
                }
            }
            return tmp_stdev;
            },  [](std::array<double, mt_kahypar::dimension> a1, std::array<double, mt_kahypar::dimension> a2) -> std::array<double, mt_kahypar::dimension>{
                std::array<double, mt_kahypar::dimension> res;
                for(int i = 0; i < mt_kahypar::dimension; i++){
                    res[i] = a1[i] + a2[i];
                }
                return res;
            });
    for(int i = 0; i < mt_kahypar::dimension; i++){
        res[i] = std::sqrt(res[i]) / static_cast<double>(n - 1);
    }
    return res;
}
    }
}