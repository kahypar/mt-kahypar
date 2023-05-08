# Guide for Implementing a New Objective Function

## Setup

- ```partition/context_enum_classes.h```: Add a new enum type to the enum classes ```Objective``` and ```GainPolicy```
- ```partition/context_enum_classes.cpp```: Create a mapping between strings and your new enum type in ```operator<< (std::ostream& os, const Objective& objective)```, ```operator<< (std::ostream& os, const GainPolicy& type)``` and ```objectiveFromString(const std::string& obj)```
- ```partition/metrics.cpp```: Create a template specialization of the ```ObjectiveFunction``` struct for your ```Objective``` enum type and override ```operator()(const PartitionedHypergraph& phg, const HyperedgeID he)```. The function takes a hypergraph and a hyperedge ID and computes the contribution of the hyperedge to the objective function. Moreover, add your new objective function to the switch statements in ```quality(...)``` and ```contribution(...)```.
- ```partition/refinement/gains/gain_definitions.h```: Create a gain type struct for your new objective function. You can copy one of the existing structures. This struct contains all relevant implementations for the gain computation in our refinement algorithms. We will later replace them with custom implementations for the new objective function. You also have to add this struct to the type list ```GainTypes``` and to the macro ```INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES```.
- ```partition/refinement/gains/gain_cache_ptr.h```: Add the ```GainPolicy``` type of your new objective function to all switch statements in the ```GainCachePtr``` class. Use the gain cache implementation of your gain type struct in ```gain_definitions.h```. We will later replace it by the concrete gain cache implementation for your new objective function. Moreover, add the ```GainPolicy``` type of your new objective function to the switch statement of the ```bipartition_each_block(...)``` function in ```partition/deep_multilevel.cpp```
- ```partition/context.cpp```: Create a mapping between the enum types ```Objective``` and ```GainPolicy``` in the ```sanityCheck(...)``` function.
- ```partition/registries/register_policies.cpp```: Create a mapping between the enum class ```GainPolicy``` and your gain type struct.

## TODOs

- make C and Python library extensible with new objective functions