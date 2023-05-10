# Guide for Implementing a New Objective Function

*TODO*: currently only minimization possible

## Setup

- ```partition/context_enum_classes.h```: Add a new enum type to the enum classes ```Objective``` and ```GainPolicy``` representing your new objective function.
- ```partition/context_enum_classes.cpp```: Create a mapping between strings and your new enum type in ```operator<< (std::ostream& os, const Objective& objective)```, ```operator<< (std::ostream& os, const GainPolicy& type)``` and ```objectiveFromString(const std::string& obj)```
- ```partition/metrics.cpp```: Create a template specialization of the ```ObjectiveFunction``` struct for your ```Objective``` enum type and override ```operator()(const PartitionedHypergraph& phg, const HyperedgeID he)```. The function takes a hypergraph and a hyperedge ID and computes the contribution of the hyperedge to the objective function. Moreover, add your new objective function to the switch statements in ```quality(...)``` and ```contribution(...)```.
- ```partition/refinement/gains/gain_definitions.h```: Create a gain type struct for your new objective function. You can copy one of the existing structures. This struct contains all relevant implementations for the gain computation in our refinement algorithms. We will later replace them with custom implementations for the new objective function. You also have to add this struct to the type list ```GainTypes``` and to the macro ```INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES```.
- ```partition/refinement/gains/gain_cache_ptr.h```: Add the ```GainPolicy``` type of your new objective function to all switch statements in the ```GainCachePtr``` class. Use the gain cache implementation of your gain type struct defined in ```gain_definitions.h```. We will later replace it by the concrete gain cache implementation for your new objective function. Moreover, add the ```GainPolicy``` type of your new objective function to the switch statement of the ```bipartition_each_block(...)``` function in ```partition/deep_multilevel.cpp```
- ```partition/context.cpp```: Create a mapping between the enum types ```Objective``` and ```GainPolicy``` in the ```sanityCheck(...)``` function.
- ```partition/registries/register_policies.cpp```: Create a mapping between the enum class ```GainPolicy``` and your gain type struct.
- ```partition/refinement/gains/bipartitioning_policy.h```: Add the ```GainPolicy``` type of your new objective function to the switch statements in ```useCutNetSplitting(...)``` and ```nonCutEdgeMultiplier(...)```. You can copy one of the existing parameters of an other objective function for now. An explanation how to configure these functions properly follows later.
- Create a folder for your objective function in ```partition/refinement/gains```. We will later add here all relevant gain computation techniques.

*TODO*: describe how the user can see that the initial setup works and describe what still fails. Moreover, describe how the new objective function can be selected from the command line interface.

## Initial Partitioning

We perform recursive bipartitioning to compute an initial k-way partition. The scheme recursively bipartitions the hypergraph until we reach the desired number of blocks. Each bipartitioning call optimizes the cut metric (weight of all cut nets). However, other objective functions can be optimized implicitly by implementing the two functions defined in ```partition/refinement/gains/bipartitioning_policy.h```. The main invariant of our recursive bipartitioning algorithm is that the cut of all bipartitions sum up to the objective value of the initial k-way partition. There are unit tests that asserts this invariant in ```tests/partition/refinement/bipartitioning_gain_policy_test.cc``` (build via ```make mt_kahypar_tests```).

### Cut Net Splitting and Removal

The recursive bipartitioning algorithm bipartitions the hypergraph and then extracts both blocks as separate hypergraphs on which we then perform the recursive calls. The block extraction algorithm requires information how deal with cut hyperedges of the bipartition. There are two options: *cut net splitting* and *cut net removal*. The former splits a hyperedge containing only the pins of the considered block, while the latter removes cut hyperedges in the extracted hypergraph. For example, if we optimize the cut metric, we can remove all cut hyperedges as they do not contribute to the cut in further bipartitions. If we optimize the connectivity metric (connectivity minus one of each cut hyperedge times their weight), we have to split cut hyperedges as increasing their connectivity can still increase the connectivity of the final k-way partition in further bipartitions.

### Non-Cut Edge Multiplier

We multiply the weight of each hyperedge that was not cut in any of the previous bipartition by this multiplier before bipartitioning a hypergraph. This feature was introduced due to a percularity of the sum-of-external-degree metric (connectivity of each cut hyperedge times their weight, also called *soed*). We can reduce the soed metric by 2 * w(e) if we remove hyperedge e from the cut (w(e) is the weight of hyperedge e). However, the soed metric only reduces by w(e) if we reduce the connectivity of e by one, but e is still cut afterwards. Therefore, we multiply the weight of all hyperedges that were not cut in any of the previous bipartitions by two. This maintains the invariant that the cut of all bipartitions sum up to the value of the soed metric of the initial k-way partition.

## Label Propagation Refinement

Our label propagation algorithm iterates over all nodes in parallel and moves each node to the block with the highest gain. The algorithm requires two gain techniques: A *gain computation* algorithm to compute the highest gain move for a node and the *attributed gain* technique to double-check the gain of a node move at the time performed on the partition. Create for both techniques a separate header file in the folder of your objective function. You can copy an existing implementation from another objective function (rename the classes appropriately). Afterwards, include both files in ```partition/refinement/gains/gain_definitions.h``` and replace the ```GainComputation``` and ```AttributedGains``` member of your gain type struct with the concrete implementations.

### Attributed Gains

The gain of a node move can change between its initial calculation and execution due to concurrent node moves in its neighborhood. We therefore double-check the gain of a node move at the time performed on the partition via synchronized data structure updates. This technique is called *attributed gains*. The label propagation algorithm reverts node moves that worsen the solution quality by checking the attributed gain value. The attributed gain function implements the following interface:
```cpp
static HyperedgeWeight gain(const HyperedgeID he,
                            const HyperedgeWeight edge_weight,
                            const HypernodeID edge_size,
                            const HypernodeID pin_count_in_from_part_after,
                            const HypernodeID pin_count_in_to_part_after);
```
When we move a node from its *source* to a *target* block, we iterate over all hyperedges, perform syncronized data structure updates and call this function for each incident hyperedge of the moved node. The sum of all calls to this function is the attributed gain of the node move. The most important parameters of this function are ```pin_count_in_from_part_after``` and ```pin_count_in_to_part_after```, which are the number of pins contained in the source and target block of hyperedge ```he``` after the move. For example, the node move removes an hyperedge from the cut if ```pin_count_in_to_part_after == edge_size```. If ```pin_count_in_from_part_after == 0```, then the node move reduces the connectivity of the hyperedge by one. Conversely, if ```pin_count_in_to_part_after == 1```, then the node move increases the connectivity of the hyperedge by one.

### Gain Computation

All gain computation techniques inherit from ```GainComputationBase``` (see ```partition/refinement/gains/gain_compute_base.h```). The base class has two template parameters: the derived class and the attributed gain implementation of the objective function (curiously recurring template pattern, avoids vtable lookups). The base class calls the ```precomputeGains(...)``` function of the derived class, which has the following interface:
```cpp
template<typename PartitionedHypergraph>
void precomputeGains(const PartitionedHypergraph& phg,
                     const HypernodeID hn,
                     RatingMap& tmp_scores,
                     Gain& isolated_block_gain);
```
We split the gain computation in two steps: (i) compute the gain of moving the node into an isolated block (stored in ```isolated_block_gain```) and (ii) moving the node from the isolated block to all adjacent blocks (stored in ```tmp_scores```). ```tmp_scores``` can be used similar to an ```std::vector``` and has exactly k entries. The gain of moving a node to a block ```to``` can then be computed by ```isolated_block_gain - tmp_scores[to]``` (a negative value means that moving the node to block ```to``` improves the objective function). However, the derived class still implements a function that computes the gain to a particular block:
```cpp
HyperedgeWeight gain(const Gain to_score, // tmp_scores[to]
                     const Gain isolated_block_gain);
```
The implementation of this function is most likely ```isolated_block_gain - to_score```. However, we plan to integrate objective functions in the future where this is not the case.

At this point, you should be able to run the ```default``` configuration of Mt-KaHyPar in debug mode without failing assertions if you disable the FM algorithm. To test this, add the following command line parameters to the Mt-KaHyPar call: ```--i-r-fm-type=do_nothing``` and ```--r-fm-type=do_nothing```. If you discover failing assertions, please check the implementations of the techniques described in the initial partitioning and label propagation section for bugs.

## TODOs

- make C and Python library extensible with new objective functions