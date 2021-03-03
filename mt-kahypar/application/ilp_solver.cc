#include <vector>

#include "gurobi_c++.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/io/command_line_options.h"
#include "mt-kahypar/partition/partitioner.h"

using namespace mt_kahypar;

class HGP_ILP {

  class ILPCallback : public GRBCallback {
   protected:
     void callback() { }
  };

 public:
  HGP_ILP(Hypergraph& hg,
          PartitionedHypergraph& initial_solution,
          const PartitionID k,
          const HypernodeWeight max_weight) :
    _hg(hg),
    _initial_solution(initial_solution),
    _k(k),
    _max_weight(max_weight),
    _env(true),
    _model(nullptr),
    _callback(),
    _variables(_hg.initialNumNodes() * _k + _hg.initialNumEdges() * _k) {
    _env.start();
    _model = std::make_unique<GRBModel>(_env);
    _model->set(GRB_IntParam_Threads, std::thread::hardware_concurrency());
    _model->setCallback(&_callback);

    // Add Model Variables
    addModelVariables();

    // Add Objective
    addObjective();

    // Ensure that each vertex is only assigned to one block
    addVertexAssignmentConstraint();

    // Ensure that balance constaint is satisfied
    addBalanceConstraint();

    // Ensure connectivity values are correct
    addConnectivityConstraint();
  }

  PartitionedHypergraph solve() {
    _model->optimize();

    PartitionedHypergraph phg(_k, _hg);
    for ( const HypernodeID& hn : _hg.nodes() ) {
      for ( PartitionID i = 0; i < _k; ++i ) {
        if ( _variables[vertex_offset(hn,i)].get(GRB_DoubleAttr_X) > 0 ) {
          phg.setOnlyNodePart(hn, i);
        }
      }
    }
    phg.initializePartition(0);

    return phg;
  }

  void printAssignment() {
    for ( const HypernodeID& hn : _hg.nodes() ) {
      for ( PartitionID i = 0; i < _k; ++i ) {
        if ( _variables[vertex_offset(hn,i)].get(GRB_DoubleAttr_X) > 0 ) {
          std::cout << "Vertex " << hn << " is assigned to block " << i << std::endl;
        }
      }
    }
  }

  void printObjective() {
    std::cout << "Connectivity Metric: " << _model->get(GRB_DoubleAttr_ObjVal) << std::endl;
  }

 private:

  void addModelVariables() {
    for ( const HypernodeID& hn : _hg.nodes() ) {
      for ( PartitionID i = 0; i < _k; ++i ) {
        _variables[vertex_offset(hn, i)] = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, vertex_var_desc(hn, i));
        _variables[vertex_offset(hn, i)].set(GRB_DoubleAttr_Start, (_initial_solution.partID(hn) == i));
      }
    }

    for ( const HypernodeID& he : _hg.edges() ) {
      for ( PartitionID i = 0; i < _k; ++i ) {
        _variables[hyperedge_offset(he, i)] = _model->addVar(0.0, 1.0, 0.0, GRB_BINARY, hyperedge_var_desc(he, i));
        _variables[hyperedge_offset(he, i)].set(GRB_DoubleAttr_Start, (_initial_solution.pinCountInPart(he, i) > 0));
      }
    }
  }

  void addObjective() {
    GRBLinExpr objective = 0;
    for ( const HypernodeID& he : _hg.edges() ) {
      GRBLinExpr connectivity = 0;
      for ( PartitionID i = 0; i < _k; ++i ) {
        connectivity += _variables[hyperedge_offset(he, i)];
      }
      objective += (connectivity - 1) * _hg.edgeWeight(he);
    }
    _model->setObjective(objective, GRB_MINIMIZE);
  }

  void addVertexAssignmentConstraint() {
    for ( const HypernodeID& hn : _hg.nodes() ) {
      GRBLinExpr constraint = 0;
      for ( PartitionID i = 0; i < _k; ++i ) {
        constraint += _variables[vertex_offset(hn, i)];
      }
      _model->addConstr(constraint == 1, "vertex_assignment_" + std::to_string(hn));
    }
  }

  void addBalanceConstraint() {
    for ( PartitionID i = 0; i < _k; ++i ) {
      GRBLinExpr constraint = 0;
      for ( const HypernodeID& hn : _hg.nodes() ) {
        constraint += _variables[vertex_offset(hn, i)] * _hg.nodeWeight(hn);
      }
      _model->addConstr(constraint <= _max_weight, "balance_constraint_" + std::to_string(i));
    }
  }

  void addConnectivityConstraint() {
    for ( const HyperedgeID& he : _hg.edges() ) {
      for ( const HypernodeID& pin : _hg.pins(he) ) {
        for ( PartitionID i = 0; i < _k; ++i ) {
          _model->addConstr(_variables[hyperedge_offset(he,i)] >= _variables[vertex_offset(pin,i)],
            "connectivity_values_" + std::to_string(he) + "_" +
            std::to_string(pin) + "_" + std::to_string(i));
        }
      }
    }
  }

  std::string vertex_var_desc(const HypernodeID hn, const PartitionID k) {
    return "x_{" + std::to_string(hn) + "," + std::to_string(k) + "}";
  }

  std::string hyperedge_var_desc(const HyperedgeID he, const PartitionID k) {
    return "y_{" + std::to_string(he) + "," + std::to_string(k) + "}";
  }

  size_t vertex_offset(const HypernodeID hn, const PartitionID k) {
    return hn * _k + k;
  }

  size_t hyperedge_offset(const HyperedgeID he, const PartitionID k) {
    return _hg.initialNumNodes() * _k + he * _k + k;
  }

  Hypergraph& _hg;
  PartitionedHypergraph& _initial_solution;
  const PartitionID _k;
  const HypernodeWeight _max_weight;
  GRBEnv _env;
  std::unique_ptr<GRBModel> _model;
  ILPCallback _callback;
  std::vector<GRBVar> _variables;
};

Hypergraph generateRandomHypergraph(const HypernodeID num_hypernodes,
                                    const HyperedgeID num_hyperedges,
                                    const HypernodeID max_edge_size) {
  parallel::scalable_vector<parallel::scalable_vector<HypernodeID>> hyperedges;
  utils::Randomize& rand = utils::Randomize::instance();
  for ( size_t i = 0; i < num_hyperedges; ++i ) {
    parallel::scalable_vector<HypernodeID> net;
    const size_t edge_size = rand.getRandomInt(2, max_edge_size, sched_getcpu());
    for ( size_t i = 0; i < edge_size; ++i ) {
      const HypernodeID pin = rand.getRandomInt(0, num_hypernodes - 1, sched_getcpu());
      if ( std::find(net.begin(), net.end(), pin) == net.end() ) {
        net.push_back(pin);
      }
    }
    hyperedges.emplace_back(std::move(net));
  }
  return HypergraphFactory::construct(
    TBBNumaArena::GLOBAL_TASK_GROUP, num_hypernodes, num_hyperedges, hyperedges);
}

int main() {

  const HypernodeID num_nodes = 100;
  const HyperedgeID num_edges = 100;
  const HypernodeID max_edge_size = 10;

  Hypergraph hg = generateRandomHypergraph(num_nodes, num_edges, max_edge_size);
  const double epsilon = 0.03;
  const PartitionID k = 2;
  const HypernodeWeight l_max = (1.0 + epsilon) * std::ceil(static_cast<double>(hg.totalWeight()) / k);

  Context context;
  parseIniToContext(context, "/home/tobias/mt-kahypar/config/fast_preset.ini");
  context.partition.k = k;
  context.partition.epsilon = epsilon;
  context.partition.objective = kahypar::Objective::km1;
  context.partition.verbose_output = false;

  // Compute Initial Solution
  PartitionedHypergraph initial_solution = partition(hg, context);
  std::cout << "Connectivity Metric: " << metrics::km1(initial_solution) << std::endl;

  try {
    HGP_ILP ilp(hg, initial_solution, k, l_max);
    PartitionedHypergraph phg = ilp.solve();
    ilp.printObjective();
    std::cout << "Connectivity Metric: " << metrics::km1(phg) << std::endl;
  } catch(GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  } catch(...) {
    std::cout << "Exception during optimization" << std::endl;
  }

  return 0;
}
