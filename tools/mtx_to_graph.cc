#include <fstream>
#include <iostream>
#include <cstdint>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <charconv>

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage ./MtxToGraph input-file output-file" << std::endl;
    std::abort();
  }
  std::string input = argv[1];
  std::string output = argv[2];

  std::ifstream in(input);
  if (!in) {
    std::cerr << "Input file " << input << " does not exist" << std::endl;
    std::abort();
  }

  uint64_t num_nodes = 0;
  uint64_t nrows = 0, ncols = 0, nnz = 0;
  std::string dummy, line, symmetry_str;
  bool symmetric = false;

  { // header
    std::getline(in, line);
    std::istringstream sstream(line);
    std::string matrix_market, object, matrix_format, data_format;
    sstream >> matrix_market >> object >> matrix_format >> data_format >> symmetry_str;

    if (matrix_format != "coordinate" || object != "matrix") {
      std::cerr << "not supported format" << std::endl;
      std::abort();
    }
    if (symmetry_str == "symmetric") {
      symmetric = true;
    } else {
      if (symmetry_str != "general") {
        std::cerr << "Unsupported symmetry option " << symmetry_str << std::endl;
        std::abort();
      }
      std::cerr << "Warning. The matrix isn't symmetric. This will result in a directed graph" << std::endl;
    }
  }

  { // dimensions
    do {
      std::getline(in, line);
    } while (line[0] == '%');
    std::istringstream sstream(line);
    sstream >> nrows >> ncols >> nnz;
    if (nrows != ncols) {
      std::cerr << "Num Rows != Num Cols --> This doesn't work for graphs." << std::endl;
      std::abort();
    }
    num_nodes = nrows;
  }

  std::cout << "num nodes = " << num_nodes << " num non-zeroes = " << nnz << " symmetric ? " << (symmetric ? "yes" : "no") << std::endl;

  std::vector<std::vector<uint32_t>> adj_list(num_nodes);

  auto t1 = std::chrono::high_resolution_clock::now();
  int row, col;
  for (uint64_t e = 0; e < nnz; ++e) {
    do {
      std::getline(in, line);
    } while (line[0] == '%');

    size_t pos = 0;
    size_t l = 0;
    while (pos < line.size() && line[pos] != ' ') { ++pos; }
    if (pos == line.size()) {  throw std::runtime_error("Line too short"); }
    std::from_chars(line.data() + l, line.data() + pos, row);

    ++pos;
    l = pos;
    while (pos < line.size() && line[pos] != ' ') { ++pos; }
    if (pos == line.size()) {  throw std::runtime_error("Line too short"); }
    std::from_chars(line.data() + l, line.data() + pos, col);

    --row; --col;
    if (row >= num_nodes || col >= num_nodes) {
      std::cerr << "Row or col index higher than number of nodes " << row << " " << col << " " << num_nodes << std::endl;
      std::cerr << line << std::endl;
      std::abort();
    }
    adj_list[row].push_back(col);
    if (symmetric) adj_list[col].push_back(row);
  }

  auto t3 = std::chrono::high_resolution_clock::now();

  std::cout << (t3-t1).count() << " reading time. " << std::endl;

  bool deg_zero = false;
  for (const auto& n : adj_list) {
    if (n.empty()) {
      deg_zero = true;
      break;
    }
  }
  if (deg_zero) {
    std::cerr << "Has zero degree nodes" << std::endl;

#ifdef false
    std::cerr << "Remap node IDs." << std::endl;
    std::vector<int64_t> remapped_node_ids(num_nodes, -1);
    uint64_t new_node_id = 0;
    for (uint64_t u = 0; u < num_nodes; ++u) {
      if (!adj_list[u].empty()) {
        adj_list[new_node_id] = std::move(adj_list[u]);
        remapped_node_ids[u] = new_node_id;
        new_node_id++;
      }
    }

    adj_list.resize(new_node_id);
    num_nodes = new_node_id;

    for (auto& n : adj_list) {
      for (auto& v : n) {
        v = remapped_node_ids[v];
      }
    }
#endif
  }

#ifdef false
  for (auto neigh : adj_list) {
    std::sort(neigh.begin(), neigh.end());
    if (std::unique(neigh.begin(), neigh.end()) != neigh.end()) {
      std::cerr << "duplicate edges..." << std::endl;
    }
  }
#endif

  uint64_t num_edges = 0;
  for (const auto& n : adj_list) num_edges += n.size();
  if (num_edges % 2 != 0) {
    std::cerr << "Num edges not even " << num_edges << std::endl;
  }
  num_edges /= 2;

  std::ofstream out(output);
  out << num_nodes << " " << num_edges << "\n";
  for (uint64_t u = 0; u < num_nodes; ++u) {
    const auto& n = adj_list[u];
    if (n.empty()) continue;
    out << n[0] + 1;
    for (size_t j = 1; j < n.size(); ++j) {
      out << " " << n[j] + 1;
    }
    out << "\n";
  }

  std::cout << "finished writing." << std::endl;
}
