// Graphene-topology: A topological vector library for generating
// graphene-derived nanostructures
//
// Author   : Zhengrong Guo (Yulin University)
// Email    : zhengrong_guo@yulinu.edu.cn
// Copyright: Zhengrong Guo (2023-2025)

#ifndef GRAPOLOGYLOCATEHPP
#define GRAPOLOGYLOCATEHPP

namespace grapology {
/**
 * @brief A class for positioning and locating atoms in graphene-derived
 * nanostructures
 *
 * This class provides functionality for dynamic bond simulation, atom
 * positioning, and structural adjustments in graphene-based materials. It
 * manages the spatial relationships between atoms (vertices) and applies
 * various forces and constraints.
 */

template <bool whether_optimizing_second_order = true>
struct locate final {
  explicit locate(const atom_container& v)
      : atoms{v},
        bond{1.42},
        next_nearest{1.42 * std::sqrt(3.0)},
        bond_damp{0.183},
        next_nearest_damp{0.183 / std::sqrt(3.0)},
        out_plane_damp{0.12} {}
  locate(const locate&) = delete;
  locate& operator=(const locate&) = delete;
  /**
   * @brief Applies dynamic relaxation to vertex positions
   * @param parameter Number of steps to relax the system
   */
  inline void dynamic_relaxing(size_t stepnum) const;
  inline void process_relaxing() const;
  inline void set_all_located() const;
  inline void spherical_position_folding() const;
  inline void set_toward_z() const;
  inline void set_edge_to_zero() const;

 private:
  const atom_container& atoms;
  const DataTypeFloat bond, next_nearest;
  const DataTypeFloat bond_damp, next_nearest_damp, out_plane_damp;
  /**
   * @brief Applies dynamic bond forces to a specific vertex
   * @param parameter Pointer to the vertex around which to calculate forces
   */
  inline void dynamic_bonds(vertex<topo_vector>* atom) const {
    for (auto nei : atom->conn_vertex) {
      auto&& bond_vector = nei->x - atom->x;
      auto&& bond_diff = bond_vector - bond_vector.normal() * bond;
      auto&& force = bond_damp * bond_diff;
      atom->f += force;
      nei->f -= force;
    }
  }
  /**
   * @brief Applies dynamic next nearest bond forces around a specific atom's
   * neighboring atoms
   * @param parameter Pointer to the atom around which to calculate forces
   */
  inline void dynamic_next_nearest(vertex<topo_vector>* atom) const {
    const auto& neighbors = atom->get_bonds();
    const auto size = neighbors.size();
    for (size_t i = 0; i < size; i++) {
      for (size_t j = i + 1; j < size; j++) {
        auto&& nei_vector = neighbors[j]->x - neighbors[i]->x;
        auto&& nei_diff = nei_vector - nei_vector.normal() * next_nearest;
        auto&& force = next_nearest_damp * nei_diff;
        neighbors[i]->f += force;
        neighbors[j]->f -= force;
      }
    }
  }

  /**
   * @brief Applies dynamic out of plane forces to a specific atom and its
   * three neighboring atoms
   * @param parameter Pointer to the atom that calculate forces toward to
   */
  inline void dynamic_out_plane(vertex<topo_vector>* atom) const {
    const auto& neighbors = atom->get_bonds();
    constexpr DataTypeFloat onethird = static_cast<DataTypeFloat>(1.0 / 3.0);
    if (neighbors.size() == 3) {
      auto&& out_vector =
          onethird * (neighbors[0]->x + neighbors[1]->x + neighbors[2]->x) -
          atom->x;
      auto&& force = out_plane_damp * out_vector;
      atom->f += force;
      force *= onethird;
      neighbors[0]->f -= force;
      neighbors[1]->f -= force;
      neighbors[2]->f -= force;
    }
  }

  void dynamic_forces(vertex<topo_vector>* atom) const {
    const DataTypeFloat bond = this->bond;
    const DataTypeFloat next_nearest = this->next_nearest;
    const DataTypeFloat bond_damp = this->bond_damp;
    const DataTypeFloat next_nearest_damp = this->next_nearest_damp;
    const DataTypeFloat out_plane_damp = this->out_plane_damp;

    const auto& neighbors = atom->get_bonds();
    const auto size = neighbors.size();
    constexpr DataTypeFloat onethird = static_cast<DataTypeFloat>(1.0 / 3.0);
    switch (size) {
      case 2: {
        auto&& bond_vector_0 = neighbors[0]->x - atom->x;
        bond_vector_0 -= bond_vector_0.normal() * bond;
        bond_vector_0 *= bond_damp;
        atom->f += bond_vector_0;
        neighbors[0]->f -= bond_vector_0;

        auto&& bond_vector_1 = neighbors[1]->x - atom->x;
        bond_vector_1 -= bond_vector_1.normal() * bond;
        bond_vector_1 *= bond_damp;
        atom->f += bond_vector_1;
        neighbors[1]->f -= bond_vector_1;

        auto&& nei_vector = neighbors[1]->x - neighbors[0]->x;
        nei_vector -= nei_vector.normal() * next_nearest;
        nei_vector *= next_nearest_damp;
        neighbors[0]->f += nei_vector;
        neighbors[1]->f -= nei_vector;
        break;
      }
      case 3: {
        auto&& bond_vector_0 = neighbors[0]->x - atom->x;
        bond_vector_0 -= bond_vector_0.normal() * bond;
        bond_vector_0 *= bond_damp;
        atom->f += bond_vector_0;
        neighbors[0]->f -= bond_vector_0;

        auto&& bond_vector_1 = neighbors[1]->x - atom->x;
        bond_vector_1 -= bond_vector_1.normal() * bond;
        bond_vector_1 *= bond_damp;
        atom->f += bond_vector_1;
        neighbors[1]->f -= bond_vector_1;

        auto&& bond_vector_2 = neighbors[2]->x - atom->x;
        bond_vector_2 -= bond_vector_2.normal() * bond;
        bond_vector_2 *= bond_damp;
        atom->f += bond_vector_2;
        neighbors[2]->f -= bond_vector_2;

        auto&& nei_vector_10 = neighbors[1]->x - neighbors[0]->x;
        nei_vector_10 -= nei_vector_10.normal() * next_nearest;
        nei_vector_10 *= next_nearest_damp;
        neighbors[0]->f += nei_vector_10;
        neighbors[1]->f -= nei_vector_10;

        auto&& nei_vector_20 = neighbors[2]->x - neighbors[0]->x;
        nei_vector_20 -= nei_vector_20.normal() * next_nearest;
        nei_vector_20 *= next_nearest_damp;
        neighbors[0]->f += nei_vector_20;
        neighbors[2]->f -= nei_vector_20;

        auto&& nei_vector_21 = neighbors[2]->x - neighbors[1]->x;
        nei_vector_21 -= nei_vector_21.normal() * next_nearest;
        nei_vector_21 *= next_nearest_damp;
        neighbors[1]->f += nei_vector_21;
        neighbors[2]->f -= nei_vector_21;

        auto&& out_vector =
            onethird * (neighbors[0]->x + neighbors[1]->x + neighbors[2]->x) -
            atom->x;
        out_vector *= out_plane_damp;
        atom->f += out_vector;
        out_vector *= onethird;
        neighbors[0]->f -= out_vector;
        neighbors[1]->f -= out_vector;
        neighbors[2]->f -= out_vector;
        break;
      }
      default:
        break;
    }
  }
  /**
   * @brief Applies dynamic intergration to a specific atom
   * @param parameter Pointer to the atom to be intergrated
   */
  inline void dynamic_intergration(vertex<topo_vector>* atom) const {
    auto sigmod = [](const DataTypeFloat x) { return 1.0 / (1.0 + x); };
    if constexpr (whether_optimizing_second_order) {
      atom->v *= 0.9;
      atom->v += atom->f;
      const auto length = atom->v.norm();
      atom->v *= sigmod(length);
      atom->x += atom->v;
    } else {
      atom->v *= 0.9;
      atom->v += atom->f;
      atom->x += atom->v;
    }
  }

  inline void dynamic_momentum(vertex<topo_vector>* atom) const {
    for (auto nei : atom->conn_vertex) {
      auto v = nei->v - atom->v;
      atom->v += v * 0.1;
      nei->v -= v * 0.1;
    }
  }
};  // class locate

template <bool W>
inline void locate<W>::dynamic_relaxing(size_t stepnum) const {
  std::vector<vertex<topo_vector>*> linear_atoms(atoms.begin(), atoms.end());
  for (auto atom : linear_atoms)
    atom->v *= 0;
  for (size_t i = 0; i < stepnum; i++) {
    for (auto atom : linear_atoms) {
      dynamic_forces(atom);
    }
    for (auto atom : linear_atoms) {
      dynamic_intergration(atom);
      atom->f *= 0;
    }
  }
}

template <bool W>
inline void locate<W>::process_relaxing() const {
  for (auto atom : atoms)
    atom->v *= 0.9;
  for (size_t i = 0; i < 5; i++) {
    for (auto atom : atoms) {
      dynamic_forces(atom);
    }
    for (auto atom : atoms) {
      //atom->v *= 0.95;
      dynamic_momentum(atom);
      atom->v += atom->f;
      const auto length = atom->v.length();
      atom->v *= 1.0 / (1.0 + length);
      atom->x += atom->v;
      atom->f *= 0;
    }
  }
}

template <bool W>
inline void locate<W>::set_all_located() const {
  for (auto atom : atoms)
    atom->whether_located = true;
}

template <bool W>
inline void locate<W>::spherical_position_folding() const {
  DataTypeFloat xlo, xhi, ylo, yhi;
  xlo = ylo = std::numeric_limits<DataTypeFloat>::max();
  xhi = yhi = std::numeric_limits<DataTypeFloat>::lowest();
  for (auto atom : atoms) {
    xlo = std::min(xlo, atom->x[0]);
    ylo = std::min(ylo, atom->x[1]);
    xhi = std::max(xhi, atom->x[0]);
    yhi = std::max(yhi, atom->x[1]);
  }
  for (auto atom : atoms) {
    atom->x[0] -= (xlo + xhi) / 2;
    atom->x[1] -= (ylo + yhi) / 2;
  }
  DataTypeFloat radial = 0.;
  for (auto atom : atoms)
    radial = radial < atom->x.length() ? atom->x.length() : radial;
  radial *= 1.2;
  const vector3<DataTypeFloat> central = {0, 0, -radial};
  for (auto atom : atoms) {
    atom->x = (atom->x - central).normal() * radial + central;
  }
}

template <bool W>
inline void locate<W>::set_toward_z() const {
  vector3<DataTypeFloat> orienta{0, 0, 0};
  const vector3<DataTypeFloat> z{0, 0, -1};
  std::set<vertex<topo_vector>*> edge_inner;

  for (auto atom : atoms) {
    if (atom->conn_vertex.size() == 2) {
      for (auto v : atom->conn_vertex) {
        if (v->conn_vertex.size() == 3) {
          orienta += atom->x - v->x;
        }
      }
    }
  }

  orienta = orienta.normal();
  auto angle = orienta.angle(z);
  auto axis = orienta ^ z;

  for (auto atom : atoms) {
    atom->x.rotate(axis, angle);
  }
}

template <bool W>
inline void locate<W>::set_edge_to_zero() const {
  vector3<DataTypeFloat> vbias{0, 0, std::numeric_limits<DataTypeFloat>::max()};
  for (auto atom : atoms) {
    if (vbias[2] > atom->x[2])
      vbias[2] = atom->x[2];
    vbias[0] += atom->x[0];
    vbias[1] += atom->x[1];
  }
  vbias[0] /= atoms.size();
  vbias[1] /= atoms.size();
  for (auto atom : atoms) {
    atom->x -= vbias;
  }
}

}  // namespace grapology

#endif