// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGYVERTEXHPP
#define GRAPOLOGYVERTEXHPP

namespace grapology {
/**
 * @class vertex
 * @brief Represents an atom in graphene-derived structures
 *
 * This class encapsulates the properties and behavior of atoms in the
 * Grapology library. It manages atomic connections, 3D positioning,
 * and simulation-related attributes like forces and velocities.
 */

template <typename T>
struct vertex final {
  using tpvector = T;
  DataTypeVector3D x, f, v;
  small_unordered_set<vertex*> conn_vertex;
  small_unordered_set<tpvector*,10> conn_vector;
  unsigned int label;
  bool whether_located;
  inline void break_bond(vertex* other);
  inline void break_all_bond();
  inline void connect_bond(vertex* other);
  inline auto get_bonds() -> std::vector<vertex<tpvector>*>;
  inline auto is_connected(vertex* other) const -> bool;
  inline void self_replaced_by(vertex* other);
  inline auto third_nei_coordinate() -> DataTypeVector3D;

  explicit vertex() : x(), f(), v(), label(0), whether_located{false} {}
  vertex(vertex&) = delete;
  vertex(vertex&&) = delete;
  ~vertex() {
    for (auto sv : conn_vertex)
      sv->conn_vertex.erase(this);
    for (auto cv : conn_vector) {
      if (cv->front_atom == this)
        cv->front_atom = nullptr;
      else if (cv->back_atom == this)
        cv->back_atom = nullptr;
    }
  }
  static void * operator new(size_t) = delete;
  static void operator delete(void*) = delete;
};  // class vertex

template <typename T>
inline void vertex<T>::break_bond(vertex* other) {
  if (this == other)
    error(np_debug_info, "error in breaking a bond but with self");
  if (nullptr == other)
    error(np_debug_info, "error in breaking a bond but with None");
  if (this->conn_vertex.contains(other)) {
    conn_vertex.erase(other);
    other->conn_vertex.erase(this);
  }
}

template <typename T>
inline void vertex<T>::break_all_bond() {
  for (auto bonded : conn_vertex) {
    bonded->conn_vertex.erase(this);
  }
  conn_vertex.clear();
}

template <typename T>
inline void vertex<T>::connect_bond(vertex* other) {
  if (this == other)
    error(np_debug_info, "error in building a bond but with self");
  if (nullptr == other)
    error(np_debug_info, "error in building a bond but with None");
  if (this->conn_vertex.contains(other) == false) {
    this->conn_vertex.insert(other);
    other->conn_vertex.insert(this);
  } else if (other->conn_vertex.contains(this) == false) {
    this->conn_vertex.insert(other);
    other->conn_vertex.insert(this);
  }
}

template <typename T>
inline auto vertex<T>::get_bonds() -> std::vector<vertex<T>*> {
  return {conn_vertex.begin(), conn_vertex.end()};
}

template <typename T>
inline auto vertex<T>::is_connected(vertex* other) const -> bool {
  return conn_vertex.contains(other);
}

template <typename T>
inline void vertex<T>::self_replaced_by(vertex* other) {
  if (this == other)
    error(np_debug_info, "error in replacing a vertice but with self");
  if (nullptr == other)
    error(np_debug_info, "error in replacing a vertice but with None");
  for (auto bonded : conn_vertex) {
    other->connect_bond(bonded);
  }
  for (auto cv : conn_vector) {
    if (cv->front_atom == this) {
      cv->front_atom = other;
      other->conn_vector.insert(cv);
    } else if (cv->back_atom == this) {
      cv->back_atom = other;
      other->conn_vector.insert(cv);
    } else
      error(np_debug_info, "error in the process of replacing an atom");
  }
  this->break_all_bond();
  conn_vector.clear();
}

template <typename T>
inline auto vertex<T>::third_nei_coordinate() -> DataTypeVector3D {
  std::vector<vertex*>&& bonds = get_bonds();
  if (bonds.size() != 2)
    error(np_debug_info,
          std::format("error in predicting the position for "
          "the third neighboring atom! current atom has {} neighbors", bonds.size()));

  auto&& va = (bonds[0]->x - x).normal();
  auto&& vb = (bonds[1]->x - x).normal();
  return x - 1.42 * (va + vb).normal();
}

}  // namespace grapology

#endif