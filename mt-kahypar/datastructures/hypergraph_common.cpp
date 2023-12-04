#include "hypergraph_common.h"
namespace mt_kahypar{
    uint32_t scalar(NodeWeight w1, NodeWeight w2){
  uint32_t res = 0;
  for(int i = 0; i < dimension; i++){
    res += w1.weights[i]*w2.weights[i];
  }
  return res;
}

std::array<double, dimension> operator*(const double d, NodeWeight w){
    std::array<double, dimension> res;
    for(int i = 0; i < dimension; i++){
      res[i] = d * w.weights[i];
    }
    return res;
  }

std::ostream& operator<<(std::ostream& os, NodeWeight nw){
    for(int i = 0; i < dimension; i++){
      os << nw.weights[i] << ' ';
    }
    return os;
  }

  void operator<<(std::ostringstream os, const std::array<double, dimension> arr){
    for(int i = 0; i < dimension; i++){
      os << arr[i] << ' ';
    }
  }

bool equals_in_one_dimension(NodeWeight w1, NodeWeight w2){
  for(int i = 0; i < dimension; i++){
    if( w1.weights[i] == w2.weights[i]){
      return true;
    }
  }
  return false;
}

std::array<double, dimension> operator+(double d, std::array<double, dimension> nw){
  std::array<double, dimension> res;
  for(int i = 0; i < dimension; i++){
    res[i] = d + static_cast<double>(nw[i]);
  }
  return res;
}


std::array<double, dimension> operator*(std::array<double, dimension> d, NodeWeight nw){
  std::array<double, dimension> res;
  for(int i = 0; i < dimension; i++){
    res[i] = d[i] * nw.weights[i];
  }
  return res;
}

bool operator<=(const std::array<double, dimension> d1, const double d2[dimension]){
  for(int i = 0; i < dimension; i++){
    if(d1[i] > d2[i]){
      return false;
    }
  }
  return true;
}

std::array<double, dimension> operator-(std::array<double, dimension> d, NodeWeight nw){
  std::array<double, dimension> res;
  for(int i = 0; i < dimension; i++){
    res[i] = d[i] - nw.weights[i];
  }
  return res;
}

std::array<double, dimension> divide_to_double(NodeWeight nw1, NodeWeight nw2){
  std::array<double, dimension> res;
  for(int i = 0; i < dimension; i++){
    res[i] = static_cast<double>(nw1.weights[i]) / static_cast<double>(nw2.weights[i]);
  }
  return res;
}
}